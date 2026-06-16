#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Shared training engine for the WCH tumor detector and WCH-FHPM.

Manifest columns
----------------
Always required: label, slide_id
Cached-patch mode: patch_path
WSI-coordinate mode: wsi_path, x, y, optional level

Training and validation manifests must be split at slide/patient level.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import time
from collections import OrderedDict
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np
import pandas as pd
from PIL import Image, ImageFile
import torch
import torch.distributed as dist
import torch.nn as nn
from sklearn.metrics import accuracy_score, roc_auc_score
from torch.nn.parallel import DistributedDataParallel as DDP
from torch.utils.data import DataLoader, Dataset
from torch.utils.data.distributed import DistributedSampler
from torchvision import transforms
import timm

try:
    import openslide
except Exception:
    openslide = None

ImageFile.LOAD_TRUNCATED_IMAGES = True


def seed_everything(seed: int, rank: int = 0) -> None:
    seed = seed + rank
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def distributed_enabled() -> bool:
    return int(os.environ.get("WORLD_SIZE", "1")) > 1


def setup_distributed() -> Tuple[int, int, int]:
    if not distributed_enabled():
        return 0, 1, 0
    dist.init_process_group(backend="nccl")
    rank = dist.get_rank()
    world_size = dist.get_world_size()
    local_rank = int(os.environ["LOCAL_RANK"])
    torch.cuda.set_device(local_rank)
    return rank, world_size, local_rank


def cleanup_distributed() -> None:
    if dist.is_available() and dist.is_initialized():
        dist.barrier()
        dist.destroy_process_group()


def reduce_sum(value: float, device: torch.device) -> float:
    tensor = torch.tensor([value], dtype=torch.float64, device=device)
    if dist.is_available() and dist.is_initialized():
        dist.all_reduce(tensor, op=dist.ReduceOp.SUM)
    return float(tensor.item())


def gather_numpy(array: np.ndarray) -> np.ndarray:
    if not (dist.is_available() and dist.is_initialized()):
        return array
    gathered = [None for _ in range(dist.get_world_size())]
    dist.all_gather_object(gathered, array)
    return np.concatenate(gathered, axis=0)


def read_manifest(path: str) -> pd.DataFrame:
    suffix = Path(path).suffix.lower()
    if suffix == ".parquet":
        df = pd.read_parquet(path)
    elif suffix == ".csv":
        df = pd.read_csv(path)
    elif suffix == ".tsv":
        df = pd.read_csv(path, sep="\t")
    else:
        raise ValueError(f"Unsupported manifest type: {path}")

    missing = {"label", "slide_id"} - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    has_patch = "patch_path" in df.columns
    has_coords = {"wsi_path", "x", "y"}.issubset(df.columns)
    if not (has_patch or has_coords):
        raise ValueError("Manifest requires patch_path, or wsi_path/x/y.")

    df = df.copy().reset_index(drop=True)
    df["label"] = pd.to_numeric(df["label"], errors="raise").astype(int)
    if not df["label"].isin([0, 1]).all():
        raise ValueError("Labels must be 0 or 1.")
    return df


class RandomRotate90:
    def __call__(self, image: Image.Image) -> Image.Image:
        return image.rotate(90 * random.randint(0, 3))


class RandomColorShift:
    def __init__(self, max_shift: float = 0.04, p: float = 0.5):
        self.max_shift = max_shift
        self.p = p

    def __call__(self, tensor: torch.Tensor) -> torch.Tensor:
        if random.random() >= self.p:
            return tensor
        shift = torch.empty(3, 1, 1).uniform_(-self.max_shift, self.max_shift)
        return torch.clamp(tensor + shift, 0.0, 1.0)


class RandomGaussianNoise:
    def __init__(self, std: float = 0.015, p: float = 0.3):
        self.std = std
        self.p = p

    def __call__(self, tensor: torch.Tensor) -> torch.Tensor:
        if random.random() >= self.p:
            return tensor
        return torch.clamp(tensor + torch.randn_like(tensor) * self.std, 0.0, 1.0)


def build_transforms(img_size: int, mean, std, training: bool):
    if training:
        return transforms.Compose([
            transforms.Resize((img_size, img_size), antialias=True),
            RandomRotate90(),
            transforms.RandomHorizontalFlip(0.5),
            transforms.RandomVerticalFlip(0.5),
            transforms.ColorJitter(brightness=0.15, contrast=0.15,
                                   saturation=0.10, hue=0.03),
            transforms.ToTensor(),
            RandomColorShift(),
            RandomGaussianNoise(),
            transforms.Normalize(mean=mean, std=std),
        ])
    return transforms.Compose([
        transforms.Resize((img_size, img_size), antialias=True),
        transforms.ToTensor(),
        transforms.Normalize(mean=mean, std=std),
    ])


class PatchManifestDataset(Dataset):
    def __init__(self, manifest_path: str, transform, tile_size: int = 512,
                 max_open_slides: int = 8):
        self.df = read_manifest(manifest_path)
        self.transform = transform
        self.tile_size = int(tile_size)
        self.max_open_slides = int(max_open_slides)
        self._slide_cache: OrderedDict[str, Any] = OrderedDict()
        self.labels = self.df["label"].to_numpy(dtype=np.int64)
        self.has_patch = "patch_path" in self.df.columns
        self.has_coords = {"wsi_path", "x", "y"}.issubset(self.df.columns)

        if self.has_patch:
            self.patch_paths = self.df["patch_path"].astype(str).to_numpy()
        if self.has_coords:
            self.wsi_paths = self.df["wsi_path"].astype(str).to_numpy()
            self.xs = pd.to_numeric(self.df["x"], errors="raise").astype(int).to_numpy()
            self.ys = pd.to_numeric(self.df["y"], errors="raise").astype(int).to_numpy()
            if "level" in self.df.columns:
                self.levels = pd.to_numeric(self.df["level"], errors="coerce").fillna(0).astype(int).to_numpy()
            else:
                self.levels = np.zeros(len(self.df), dtype=np.int64)

    def __len__(self):
        return len(self.df)

    def __getstate__(self):
        state = self.__dict__.copy()
        state["_slide_cache"] = OrderedDict()
        return state

    def _get_slide(self, path: str):
        if openslide is None:
            raise RuntimeError("openslide-python is required for WSI reading.")
        if path in self._slide_cache:
            slide = self._slide_cache.pop(path)
            self._slide_cache[path] = slide
            return slide
        slide = openslide.OpenSlide(path)
        self._slide_cache[path] = slide
        while len(self._slide_cache) > self.max_open_slides:
            _, old = self._slide_cache.popitem(last=False)
            old.close()
        return slide

    def _load_image(self, index: int) -> Image.Image:
        if self.has_patch:
            with Image.open(self.patch_paths[index]) as image:
                return image.convert("RGB")
        slide = self._get_slide(self.wsi_paths[index])
        return slide.read_region(
            (int(self.xs[index]), int(self.ys[index])),
            int(self.levels[index]),
            (self.tile_size, self.tile_size),
        ).convert("RGB")

    def __getitem__(self, index: int):
        image = self.transform(self._load_image(index))
        return image, int(self.labels[index])


def create_model(model_name: str, img_size: int, init_mode: str,
                 init_checkpoint: str) -> nn.Module:
    model = timm.create_model(
        model_name,
        pretrained=(init_mode == "imagenet"),
        img_size=img_size,
        num_classes=2,
    )
    if init_mode == "checkpoint":
        if not init_checkpoint:
            raise ValueError("--init_checkpoint is required.")
        ckpt = torch.load(init_checkpoint, map_location="cpu", weights_only=False)
        state = ckpt.get("state_dict") or ckpt.get("model_state_dict") or ckpt.get("model") or ckpt
        state = {k.replace("module.", "", 1) if k.startswith("module.") else k: v
                 for k, v in state.items()}
        missing, unexpected = model.load_state_dict(state, strict=False)
        print("Checkpoint initialization; missing:", missing)
        print("Checkpoint initialization; unexpected:", unexpected)
    return model


def unwrap_model(model: nn.Module) -> nn.Module:
    return model.module if isinstance(model, DDP) else model


def inference_payload(task: str, model: nn.Module, args) -> Dict[str, Any]:
    state = {k: v.detach().cpu() for k, v in unwrap_model(model).state_dict().items()}
    preprocess = {"mean": args.mean, "std": args.std,
                  "resize": [args.img_size, args.img_size]}
    if task == "tumor":
        return {
            "task": "tumor_detector",
            "version": args.checkpoint_version,
            "model_name": args.model_name,
            "img_size": args.img_size,
            "num_classes": 2,
            "label_mapping": {"normal_tissue": 0, "tumor_associated_tissue": 1},
            "preprocess": preprocess,
            "recommended_inference_rule": {
                "normal_probability_threshold": 0.80,
                "equivalent_tumor_probability_threshold": 0.20,
                "rule": "normal if P(normal)>0.80; otherwise tumor-associated",
            },
            "model_state_dict": state,
        }
    return {
        "format_version": "1.0",
        "model_name": args.model_name,
        "img_size": args.img_size,
        "num_classes": 2,
        "state_dict": state,
        "preprocess": preprocess,
        "label_mapping": {"nonFH": 0, "FH": 1},
        "positive_class": "FH",
        "default_cutoff": 0.5,
        "slide_level_aggregation": "mean_patch_probability",
        "input_description": "512x512 tumor-associated H&E image patches",
        "intended_use": "research_use_only",
    }


def save_training_state(path: Path, epoch: int, model: nn.Module,
                        optimizer, scheduler, scaler, best_auc: float, args) -> None:
    torch.save({
        "epoch": epoch,
        "model_state_dict": unwrap_model(model).state_dict(),
        "optimizer_state_dict": optimizer.state_dict(),
        "scheduler_state_dict": scheduler.state_dict(),
        "scaler_state_dict": scaler.state_dict(),
        "best_auc": best_auc,
        "args": vars(args),
    }, path)


def train_one_epoch(model, loader, criterion, optimizer, scaler, device,
                    grad_accum_steps: int, rank: int, log_interval: int):
    model.train()
    optimizer.zero_grad(set_to_none=True)
    total_loss = total_correct = total_items = 0.0
    started = time.time()

    for step, (images, labels) in enumerate(loader, start=1):
        images = images.to(device, non_blocking=True)
        labels = labels.to(device, non_blocking=True)
        with torch.autocast(device_type="cuda", dtype=torch.float16,
                            enabled=(device.type == "cuda")):
            logits = model(images)
            loss = criterion(logits, labels)
        scaler.scale(loss / grad_accum_steps).backward()
        if step % grad_accum_steps == 0 or step == len(loader):
            scaler.step(optimizer)
            scaler.update()
            optimizer.zero_grad(set_to_none=True)

        total_loss += float(loss.item()) * labels.numel()
        total_correct += float((logits.argmax(1) == labels).sum().item())
        total_items += float(labels.numel())

        if rank == 0 and step % log_interval == 0:
            print(f"step {step}/{len(loader)} | "
                  f"loss={total_loss/max(total_items,1):.5f} | "
                  f"acc={total_correct/max(total_items,1):.5f} | "
                  f"elapsed={(time.time()-started)/60:.1f} min", flush=True)

    total_loss = reduce_sum(total_loss, device)
    total_correct = reduce_sum(total_correct, device)
    total_items = reduce_sum(total_items, device)
    return {"loss": total_loss/max(total_items, 1.0),
            "accuracy": total_correct/max(total_items, 1.0)}


@torch.inference_mode()
def validate(model, loader, criterion, device):
    model.eval()
    local_loss = local_items = 0.0
    probs_all, labels_all = [], []
    for images, labels in loader:
        images = images.to(device, non_blocking=True)
        labels = labels.to(device, non_blocking=True)
        with torch.autocast(device_type="cuda", dtype=torch.float16,
                            enabled=(device.type == "cuda")):
            logits = model(images)
            loss = criterion(logits, labels)
        probs_all.append(torch.softmax(logits, 1)[:, 1].float().cpu().numpy())
        labels_all.append(labels.cpu().numpy())
        local_loss += float(loss.item()) * labels.numel()
        local_items += float(labels.numel())

    probs = np.concatenate(probs_all) if probs_all else np.empty(0, dtype=np.float32)
    labels = np.concatenate(labels_all) if labels_all else np.empty(0, dtype=np.int64)
    probs, labels = gather_numpy(probs), gather_numpy(labels)
    total_loss = reduce_sum(local_loss, device)
    total_items = reduce_sum(local_items, device)

    if labels.size == 0:
        return {"loss": float("nan"), "accuracy": float("nan"), "auroc": float("nan")}
    accuracy = accuracy_score(labels, (probs >= 0.5).astype(np.int64))
    try:
        auroc = roc_auc_score(labels, probs)
    except ValueError:
        auroc = float("nan")
    return {"loss": total_loss/max(total_items, 1.0),
            "accuracy": float(accuracy), "auroc": float(auroc)}


def append_history(path: Path, row: Dict[str, Any]) -> None:
    exists = path.exists()
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def build_parser(task: str):
    parser = argparse.ArgumentParser(description=f"Train WCH task: {task}")
    parser.add_argument("--train_manifest", required=True)
    parser.add_argument("--val_manifest", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--model_name", default="vit_base_patch16_224")
    parser.add_argument("--img_size", type=int, default=512)
    parser.add_argument("--tile_size", type=int, default=512)
    parser.add_argument("--init", choices=["imagenet", "random", "checkpoint"], default="imagenet")
    parser.add_argument("--init_checkpoint", default="")
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--batch_size", type=int, default=16, help="Per GPU/process")
    parser.add_argument("--grad_accum_steps", type=int, default=6)
    parser.add_argument("--num_workers", type=int, default=8)
    parser.add_argument("--learning_rate", type=float, default=1e-4)
    parser.add_argument("--weight_decay", type=float, default=0.0)
    parser.add_argument("--min_learning_rate", type=float, default=1e-6)
    parser.add_argument("--class_weights", nargs=2, type=float, default=None)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--log_interval", type=int, default=100)
    parser.add_argument("--resume", default="")
    parser.add_argument("--checkpoint_version", default="WCH_reference_training")
    parser.add_argument("--mean", nargs=3, type=float, default=[0.485, 0.456, 0.406])
    parser.add_argument("--std", nargs=3, type=float, default=[0.229, 0.224, 0.225])
    return parser


def main(task: str) -> None:
    args = build_parser(task).parse_args()
    rank, world_size, local_rank = setup_distributed()
    seed_everything(args.seed, rank)
    device = torch.device(f"cuda:{local_rank}" if distributed_enabled() else
                          ("cuda" if torch.cuda.is_available() else "cpu"))

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    if rank == 0:
        (outdir / "training_config.json").write_text(
            json.dumps(vars(args), indent=2, ensure_ascii=False), encoding="utf-8")

    train_ds = PatchManifestDataset(
        args.train_manifest,
        build_transforms(args.img_size, args.mean, args.std, True),
        args.tile_size,
    )
    val_ds = PatchManifestDataset(
        args.val_manifest,
        build_transforms(args.img_size, args.mean, args.std, False),
        args.tile_size,
    )

    train_sampler = DistributedSampler(train_ds, num_replicas=world_size, rank=rank,
                                       shuffle=True, seed=args.seed) if distributed_enabled() else None
    val_sampler = DistributedSampler(val_ds, num_replicas=world_size, rank=rank,
                                     shuffle=False, drop_last=False) if distributed_enabled() else None

    loader_args = dict(batch_size=args.batch_size, num_workers=args.num_workers,
                       pin_memory=torch.cuda.is_available(),
                       persistent_workers=args.num_workers > 0)
    train_loader = DataLoader(train_ds, shuffle=train_sampler is None,
                              sampler=train_sampler, drop_last=False, **loader_args)
    val_loader = DataLoader(val_ds, shuffle=False, sampler=val_sampler,
                            drop_last=False, **loader_args)

    model = create_model(args.model_name, args.img_size, args.init, args.init_checkpoint).to(device)
    if distributed_enabled():
        model = DDP(model, device_ids=[local_rank], output_device=local_rank,
                    broadcast_buffers=False)

    weights = None if args.class_weights is None else torch.tensor(
        args.class_weights, dtype=torch.float32, device=device)
    criterion = nn.CrossEntropyLoss(weight=weights)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.learning_rate,
                                  weight_decay=args.weight_decay)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=args.epochs, eta_min=args.min_learning_rate)
    scaler = torch.cuda.amp.GradScaler(enabled=(device.type == "cuda"))

    start_epoch, best_auc = 1, -math.inf
    if args.resume:
        state = torch.load(args.resume, map_location="cpu", weights_only=False)
        unwrap_model(model).load_state_dict(state["model_state_dict"], strict=True)
        optimizer.load_state_dict(state["optimizer_state_dict"])
        scheduler.load_state_dict(state["scheduler_state_dict"])
        scaler.load_state_dict(state.get("scaler_state_dict", {}))
        start_epoch = int(state["epoch"]) + 1
        best_auc = float(state.get("best_auc", -math.inf))

    if rank == 0:
        effective = args.batch_size * world_size * args.grad_accum_steps
        print(f"Task={task}; train={len(train_ds):,}; val={len(val_ds):,}; "
              f"effective_global_batch={effective}")

    history_path = outdir / "training_history.csv"
    for epoch in range(start_epoch, args.epochs + 1):
        if train_sampler is not None:
            train_sampler.set_epoch(epoch)
        if rank == 0:
            print(f"\nEpoch {epoch}/{args.epochs}", flush=True)

        train_metrics = train_one_epoch(
            model, train_loader, criterion, optimizer, scaler, device,
            args.grad_accum_steps, rank, args.log_interval)
        val_metrics = validate(model, val_loader, criterion, device)
        current_lr = float(optimizer.param_groups[0]["lr"])
        scheduler.step()

        if rank == 0:
            row = {
                "epoch": epoch,
                "learning_rate": current_lr,
                "train_loss": train_metrics["loss"],
                "train_accuracy": train_metrics["accuracy"],
                "val_loss": val_metrics["loss"],
                "val_accuracy": val_metrics["accuracy"],
                "val_auroc": val_metrics["auroc"],
            }
            append_history(history_path, row)
            print(json.dumps(row, indent=2), flush=True)
            save_training_state(outdir / "last_training_state.pt", epoch, model,
                                optimizer, scheduler, scaler, best_auc, args)

            auc = val_metrics["auroc"]
            if not math.isnan(auc) and auc > best_auc:
                best_auc = auc
                filename = ("WCH_TumorDetector_ViTBase512.pth" if task == "tumor"
                            else "WCH_FHPM_ViTBase512.pth")
                torch.save(inference_payload(task, model, args), outdir / filename)
                save_training_state(outdir / "best_training_state.pt", epoch, model,
                                    optimizer, scheduler, scaler, best_auc, args)
                print(f"Saved new best checkpoint: {filename}; AUROC={best_auc:.6f}")

        if dist.is_available() and dist.is_initialized():
            dist.barrier()

    cleanup_distributed()


if __name__ == "__main__":
    raise RuntimeError("Run train_tumor_detector.py or train_wch_fhpm.py")

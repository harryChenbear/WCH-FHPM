#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Complete inference workflow for WCH-FHPM.

This script processes H&E-stained whole-slide images (WSIs) and produces:
1. tissue-filtered patch coordinates,
2. tumor-associated patch detection using WCH_TumorDetector_ViTBase512.pth,
3. FH-dRCC patch-level probabilities using WCH_FHPM_ViTBase512.pth,
4. slide-level FH-dRCC probability by averaging tumor-associated patch probabilities,
5. patch-level CSV outputs and a heatmap PNG.

The default tumor-associated threshold is 0.20, which is equivalent to excluding
patches with predicted normal-tissue probability > 0.80 in a two-class softmax
model. This matches the manuscript description: highly normal patches are
removed, while ambiguous, mixed, stromal, and tumor-adjacent patches are retained.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw

import torch
import torch.nn.functional as F
from torchvision import transforms

try:
    import openslide
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "openslide-python is required. Install it with `pip install openslide-python` "
        "and install the OpenSlide system library, e.g. `apt-get install -y openslide-tools libopenslide0`."
    ) from exc

try:
    import timm
except ImportError as exc:  # pragma: no cover
    raise ImportError("timm is required. Install it with `pip install timm`.") from exc


SUPPORTED_EXTENSIONS = {
    ".svs",
    ".ndpi",
    ".tif",
    ".tiff",
    ".mrxs",
    ".scn",
    ".kfb",
}

IMAGENET_MEAN = (0.485, 0.456, 0.406)
IMAGENET_STD = (0.229, 0.224, 0.225)


# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the complete WCH-FHPM inference workflow on one WSI or a folder of WSIs."
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--slide", type=str, help="Path to one whole-slide image.")
    input_group.add_argument("--slide_dir", type=str, help="Path to a folder containing whole-slide images.")

    parser.add_argument("--out_dir", type=str, required=True, help="Output directory.")
    parser.add_argument("--device", type=str, default="cuda", help="Inference device: cuda or cpu. Default: cuda.")

    parser.add_argument("--tumor_model", type=str, default="models/WCH_TumorDetector_ViTBase512.pth")
    parser.add_argument("--fh_model", type=str, default="models/WCH_FHPM_ViTBase512.pth")

    parser.add_argument("--tile_size", type=int, default=512, help="Patch size in pixels. Default: 512.")
    parser.add_argument("--step_size", type=int, default=512, help="Step size in pixels. Default: 512.")
    parser.add_argument(
        "--min_tissue_fraction",
        type=float,
        default=0.20,
        help="Minimum tissue fraction for retaining a patch. Default: 0.20.",
    )
    parser.add_argument(
        "--tumor_threshold",
        type=float,
        default=0.20,
        help=(
            "Threshold for tumor-associated patch selection. Default: 0.20. "
            "For the two-class tumor detector, this corresponds to excluding patches "
            "with predicted normal-tissue probability > 0.80."
        ),
    )
    parser.add_argument(
        "--fh_cutoff",
        type=float,
        default=0.50,
        help="Cutoff for slide-level FH-dRCC classification. Default: 0.50.",
    )
    parser.add_argument("--read_batch_size", type=int, default=128)
    parser.add_argument("--batch_size", type=int, default=128)
    parser.add_argument("--num_workers", type=int, default=0, help="Reserved for compatibility. Default: 0.")
    parser.add_argument(
        "--config",
        type=str,
        default="model_config.json",
        help="Optional JSON configuration file to copy into the output folder. Default: model_config.json.",
    )
    parser.add_argument(
        "--save_heatmap",
        action="store_true",
        default=True,
        help="Save a patch-level FH probability heatmap. Default: True.",
    )
    parser.add_argument(
        "--no_heatmap",
        action="store_false",
        dest="save_heatmap",
        help="Disable heatmap generation.",
    )
    return parser.parse_args()


def resolve_path(path_str: str, base_dir: Optional[Path] = None) -> Path:
    path = Path(path_str)
    if path.is_absolute():
        return path
    if base_dir is None:
        base_dir = Path.cwd()
    return (base_dir / path).resolve()


def find_slides(slide_dir: Path) -> List[Path]:
    slides: List[Path] = []
    for root, _, files in os.walk(slide_dir):
        for filename in files:
            path = Path(root) / filename
            if path.suffix.lower() in SUPPORTED_EXTENSIONS:
                slides.append(path)
    return sorted(slides)


def safe_slide_name(slide_path: Path) -> str:
    name = slide_path.stem
    safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in name)
    return safe or "slide"


def get_tissue_fraction(tile: Image.Image) -> float:
    """Estimate tissue fraction using a conservative RGB/HSV-based background filter."""
    arr = np.asarray(tile.convert("RGB"), dtype=np.uint8)
    if arr.size == 0:
        return 0.0

    # Convert to HSV using PIL to avoid adding extra dependencies.
    hsv = np.asarray(tile.convert("HSV"), dtype=np.uint8)
    saturation = hsv[:, :, 1]
    value = hsv[:, :, 2]

    # White background usually has high value and low saturation.
    # This rule retains colored tissue and excludes bright white background.
    tissue_mask = (saturation > 20) & (value < 245)

    # Very dark regions may correspond to pen marks or scanning artifacts; keep them
    # only if they also have some saturation.
    tissue_fraction = float(tissue_mask.mean())
    return tissue_fraction


def read_tile(slide: "openslide.OpenSlide", x: int, y: int, tile_size: int) -> Image.Image:
    """Read a level-0 tile and return an RGB image."""
    tile = slide.read_region((int(x), int(y)), 0, (tile_size, tile_size)).convert("RGB")
    return tile


def build_transform(tile_size: int) -> transforms.Compose:
    return transforms.Compose(
        [
            transforms.Resize((tile_size, tile_size)),
            transforms.ToTensor(),
            transforms.Normalize(mean=IMAGENET_MEAN, std=IMAGENET_STD),
        ]
    )


def strip_state_dict_prefix(state_dict: Dict[str, torch.Tensor]) -> Dict[str, torch.Tensor]:
    """Remove common wrappers such as DataParallel's 'module.' prefix."""
    cleaned = {}
    for key, value in state_dict.items():
        new_key = key
        for prefix in ("module.", "model.", "net."):
            if new_key.startswith(prefix):
                new_key = new_key[len(prefix) :]
        cleaned[new_key] = value
    return cleaned


def extract_state_dict(checkpoint: object) -> Dict[str, torch.Tensor]:
    """Extract a PyTorch state_dict from common checkpoint formats."""
    if isinstance(checkpoint, dict):
        for key in ("model_state_dict", "state_dict", "model", "net"):
            value = checkpoint.get(key)
            if isinstance(value, dict):
                return strip_state_dict_prefix(value)
        if all(isinstance(k, str) for k in checkpoint.keys()):
            # The checkpoint itself may already be a state_dict.
            if any(isinstance(v, torch.Tensor) for v in checkpoint.values()):
                return strip_state_dict_prefix(checkpoint)  # type: ignore[arg-type]
    raise ValueError(
        "Could not find a model state_dict in the checkpoint. Expected one of: "
        "model_state_dict, state_dict, model, net, or a raw state_dict."
    )


def load_checkpoint(path: Path, map_location: torch.device) -> object:
    try:
        return torch.load(path, map_location=map_location, weights_only=False)
    except TypeError:
        return torch.load(path, map_location=map_location)


def create_vit_base_model(num_classes: int = 2, img_size: int = 512) -> torch.nn.Module:
    """Create the ViT-Base classifier used by the WCH-FHPM workflow."""
    model = timm.create_model(
        "vit_base_patch16_224",
        pretrained=False,
        num_classes=num_classes,
        img_size=img_size,
    )
    return model


def load_model(model_path: Path, device: torch.device, img_size: int = 512) -> torch.nn.Module:
    if not model_path.exists():
        raise FileNotFoundError(
            f"Model weight file not found: {model_path}\n"
            "Download the required weight files from the release page and place them under the models/ folder."
        )
    checkpoint = load_checkpoint(model_path, map_location=device)
    state_dict = extract_state_dict(checkpoint)
    model = create_vit_base_model(num_classes=2, img_size=img_size)
    missing, unexpected = model.load_state_dict(state_dict, strict=False)
    if unexpected:
        print(f"[Warning] Unexpected checkpoint keys in {model_path.name}: {unexpected[:10]}")
    if missing:
        print(f"[Warning] Missing checkpoint keys in {model_path.name}: {missing[:10]}")
    model.to(device)
    model.eval()
    return model


def infer_probabilities(
    model: torch.nn.Module,
    images: Sequence[Image.Image],
    transform: transforms.Compose,
    device: torch.device,
    batch_size: int,
    positive_index: int = 1,
) -> np.ndarray:
    """Run model inference and return class-1 probabilities."""
    if len(images) == 0:
        return np.asarray([], dtype=np.float32)

    probs: List[np.ndarray] = []
    with torch.no_grad():
        for start in range(0, len(images), batch_size):
            batch_imgs = images[start : start + batch_size]
            tensor = torch.stack([transform(img) for img in batch_imgs], dim=0).to(device)
            logits = model(tensor)
            prob = F.softmax(logits, dim=1)[:, positive_index]
            probs.append(prob.detach().cpu().numpy().astype(np.float32))
    return np.concatenate(probs, axis=0)


def write_json(path: Path, obj: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)


def make_heatmap(
    records: pd.DataFrame,
    slide_width: int,
    slide_height: int,
    tile_size: int,
    step_size: int,
    out_path: Path,
) -> None:
    """Create a simple patch-grid FH probability heatmap using PIL only."""
    if records.empty or "FH_probability" not in records.columns:
        return

    n_cols = max(1, math.ceil(max(slide_width - tile_size, 0) / step_size) + 1)
    n_rows = max(1, math.ceil(max(slide_height - tile_size, 0) / step_size) + 1)

    # Keep heatmap reasonably sized.
    max_side = 2000
    cell = max(1, min(max_side // max(n_cols, 1), max_side // max(n_rows, 1)))
    width = n_cols * cell
    height = n_rows * cell

    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)

    for _, row in records.iterrows():
        try:
            fh_prob = float(row["FH_probability"])
        except (TypeError, ValueError):
            continue
        if not np.isfinite(fh_prob):
            continue

        x = int(row["x"])
        y = int(row["y"])
        col = max(0, min(n_cols - 1, x // step_size))
        rr = max(0, min(n_rows - 1, y // step_size))

        # Blue-white-red style without matplotlib dependency.
        # Low FH: blue; high FH: red.
        p = max(0.0, min(1.0, fh_prob))
        if p < 0.5:
            t = p / 0.5
            r = int(255 * t)
            g = int(255 * t)
            b = 255
        else:
            t = (p - 0.5) / 0.5
            r = 255
            g = int(255 * (1.0 - t))
            b = int(255 * (1.0 - t))
        draw.rectangle(
            [col * cell, rr * cell, (col + 1) * cell - 1, (rr + 1) * cell - 1],
            fill=(r, g, b),
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    img.save(out_path)


# -----------------------------------------------------------------------------
# Main slide-processing workflow
# -----------------------------------------------------------------------------

def generate_grid(width: int, height: int, tile_size: int, step_size: int) -> Iterable[Tuple[int, int]]:
    max_x = max(0, width - tile_size)
    max_y = max(0, height - tile_size)
    ys = list(range(0, max_y + 1, step_size))
    xs = list(range(0, max_x + 1, step_size))
    if ys[-1] != max_y:
        ys.append(max_y)
    if xs[-1] != max_x:
        xs.append(max_x)
    for y in ys:
        for x in xs:
            yield x, y


def process_slide(
    slide_path: Path,
    out_dir: Path,
    tumor_model: torch.nn.Module,
    fh_model: torch.nn.Module,
    transform: transforms.Compose,
    device: torch.device,
    tile_size: int,
    step_size: int,
    min_tissue_fraction: float,
    tumor_threshold: float,
    fh_cutoff: float,
    read_batch_size: int,
    batch_size: int,
    save_heatmap: bool,
) -> Dict[str, object]:
    start_time = time.time()
    slide_name = safe_slide_name(slide_path)
    slide_out_dir = out_dir / slide_name
    slide_out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n[Slide] {slide_path}")
    print(f"[Output] {slide_out_dir}")

    records: List[Dict[str, object]] = []
    n_grid_patches = 0
    n_tissue_patches = 0
    n_tumor_patches = 0

    slide = openslide.OpenSlide(str(slide_path))
    width, height = slide.dimensions

    batch_tiles: List[Image.Image] = []
    batch_meta: List[Tuple[int, int, float]] = []

    def flush_batch() -> None:
        nonlocal n_tumor_patches, records, batch_tiles, batch_meta
        if not batch_tiles:
            return

        tumor_probs = infer_probabilities(
            tumor_model,
            batch_tiles,
            transform,
            device,
            batch_size=batch_size,
            positive_index=1,
        )
        is_tumor = tumor_probs >= tumor_threshold

        tumor_tiles = [img for img, keep in zip(batch_tiles, is_tumor) if bool(keep)]
        fh_probs_for_tumor = infer_probabilities(
            fh_model,
            tumor_tiles,
            transform,
            device,
            batch_size=batch_size,
            positive_index=1,
        )

        tumor_idx = 0
        for (x, y, tissue_fraction), tumor_prob, keep in zip(batch_meta, tumor_probs, is_tumor):
            if bool(keep):
                fh_prob = float(fh_probs_for_tumor[tumor_idx])
                tumor_idx += 1
                n_tumor_patches += 1
            else:
                fh_prob = np.nan

            records.append(
                {
                    "slide_name": slide_path.name,
                    "x": int(x),
                    "y": int(y),
                    "tissue_fraction": float(tissue_fraction),
                    "tumor_probability": float(tumor_prob),
                    "is_tumor_associated": bool(keep),
                    "FH_probability": fh_prob,
                }
            )

        batch_tiles = []
        batch_meta = []

    for x, y in generate_grid(width, height, tile_size, step_size):
        n_grid_patches += 1
        tile = read_tile(slide, x, y, tile_size)
        tissue_fraction = get_tissue_fraction(tile)
        if tissue_fraction < min_tissue_fraction:
            continue

        n_tissue_patches += 1
        batch_tiles.append(tile)
        batch_meta.append((x, y, tissue_fraction))

        if len(batch_tiles) >= read_batch_size:
            flush_batch()

    flush_batch()
    slide.close()

    patch_df = pd.DataFrame.from_records(records)
    patch_csv = slide_out_dir / "patch_predictions.csv"
    patch_df.to_csv(patch_csv, index=False)

    if n_tumor_patches > 0:
        final_prob = float(patch_df.loc[patch_df["is_tumor_associated"], "FH_probability"].mean())
        final_label = int(final_prob >= fh_cutoff)
        final_prediction = "FH-dRCC" if final_label == 1 else "non-FH-dRCC"
        status = "success"
    else:
        final_prob = np.nan
        final_label = np.nan
        final_prediction = "NA"
        status = "no_tumor_associated_patches"

    elapsed_minutes = (time.time() - start_time) / 60.0
    result = {
        "slide_name": slide_path.name,
        "slide_path": str(slide_path),
        "status": status,
        "slide_width": int(width),
        "slide_height": int(height),
        "n_grid_patches": int(n_grid_patches),
        "n_tissue_patches": int(n_tissue_patches),
        "n_tumor_patches": int(n_tumor_patches),
        "tissue_fraction_among_grid": float(n_tissue_patches / n_grid_patches) if n_grid_patches else np.nan,
        "tumor_fraction_among_tissue": float(n_tumor_patches / n_tissue_patches) if n_tissue_patches else np.nan,
        "final_fh_probability": final_prob,
        "final_prediction": final_prediction,
        f"final_pred_label_{fh_cutoff}": final_label,
        "elapsed_minutes": float(elapsed_minutes),
    }

    pd.DataFrame([result]).to_csv(slide_out_dir / "slide_result.csv", index=False)

    if save_heatmap and not patch_df.empty:
        make_heatmap(
            patch_df,
            slide_width=width,
            slide_height=height,
            tile_size=tile_size,
            step_size=step_size,
            out_path=slide_out_dir / "FH_probability_heatmap.png",
        )

    print(
        f"[Done] {slide_path.name}: status={status}, "
        f"tumor_patches={n_tumor_patches}, final_fh_probability={final_prob}"
    )
    return result


def main() -> None:
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Resolve default model paths relative to the script directory, so the command
    # can be launched from either the repository folder or another working directory.
    tumor_model_path = resolve_path(args.tumor_model, base_dir=script_dir)
    fh_model_path = resolve_path(args.fh_model, base_dir=script_dir)

    if args.device == "cuda" and not torch.cuda.is_available():
        print("[Warning] CUDA was requested but is not available. Falling back to CPU.")
        device = torch.device("cpu")
    else:
        device = torch.device(args.device)

    if args.slide:
        slides = [Path(args.slide).resolve()]
    else:
        slides = find_slides(Path(args.slide_dir).resolve())

    if not slides:
        raise FileNotFoundError("No supported whole-slide image files were found.")

    print(f"[Info] Number of slides: {len(slides)}")
    print(f"[Info] Device: {device}")
    print(f"[Info] Tumor detector: {tumor_model_path}")
    print(f"[Info] FH predictor: {fh_model_path}")
    print(
        f"[Info] tumor_threshold={args.tumor_threshold}. "
        f"For a two-class softmax tumor detector, this excludes patches with "
        f"predicted normal-tissue probability > {1.0 - args.tumor_threshold:.2f}."
    )

    transform = build_transform(args.tile_size)
    tumor_model = load_model(tumor_model_path, device=device, img_size=args.tile_size)
    fh_model = load_model(fh_model_path, device=device, img_size=args.tile_size)

    run_config = vars(args).copy()
    run_config.update(
        {
            "resolved_tumor_model": str(tumor_model_path),
            "resolved_fh_model": str(fh_model_path),
            "device_resolved": str(device),
            "n_slides": len(slides),
            "tumor_threshold_note": (
                "A tumor-associated probability threshold of 0.20 corresponds to excluding patches "
                "with predicted normal-tissue probability greater than 0.80."
            ),
        }
    )
    write_json(out_dir / "run_config.json", run_config)

    # Copy optional model_config.json content into the output folder if available.
    config_path = resolve_path(args.config, base_dir=script_dir)
    if config_path.exists():
        try:
            with open(config_path, "r", encoding="utf-8") as f:
                config_obj = json.load(f)
            write_json(out_dir / "model_config_used.json", config_obj)
        except Exception as exc:
            print(f"[Warning] Could not copy model_config.json: {exc}")

    all_results: List[Dict[str, object]] = []
    for slide_path in slides:
        try:
            result = process_slide(
                slide_path=slide_path,
                out_dir=out_dir,
                tumor_model=tumor_model,
                fh_model=fh_model,
                transform=transform,
                device=device,
                tile_size=args.tile_size,
                step_size=args.step_size,
                min_tissue_fraction=args.min_tissue_fraction,
                tumor_threshold=args.tumor_threshold,
                fh_cutoff=args.fh_cutoff,
                read_batch_size=args.read_batch_size,
                batch_size=args.batch_size,
                save_heatmap=args.save_heatmap,
            )
        except Exception as exc:
            print(f"[Error] Failed to process {slide_path}: {exc}")
            result = {
                "slide_name": slide_path.name,
                "slide_path": str(slide_path),
                "status": "error",
                "error_message": str(exc),
            }
        all_results.append(result)
        pd.DataFrame(all_results).to_csv(out_dir / "all_slide_results.csv", index=False)

    print(f"\n[Finished] Summary saved to: {out_dir / 'all_slide_results.csv'}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
import argparse
from pathlib import Path
import sys
ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))
import pandas as pd
import torch
from abmil.model import GatedABMIL
from abmil.io import load_features, torch_load
from abmil.metrics import binary_metrics, bootstrap_ci


def load_model(model_path, device):
    ckpt = torch_load(model_path, map_location=device)
    config = ckpt.get("config", {})
    model = GatedABMIL(
        input_dim=int(config.get("input_dim", 1024)),
        hidden_dim=int(config.get("hidden_dim", 512)),
        attn_dim=int(config.get("attn_dim", 256)),
        dropout=float(config.get("dropout", 0.25)),
    ).to(device)
    model.load_state_dict(ckpt["model_state_dict"])
    model.eval()
    return model, config


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--dataset_name", default="External")
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    device = args.device if (args.device == "cpu" or torch.cuda.is_available()) else "cpu"
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.manifest)
    df["label"] = df["label"].astype(int)
    model, config = load_model(args.model, device)
    input_dim = int(config.get("input_dim", 1024))
    print("Device:", device)
    print(f"{args.dataset_name} slides:", len(df))
    print(df["label_name"].value_counts())
    probs = []
    with torch.no_grad():
        for i, row in df.iterrows():
            x = load_features(row["feature_path"], expected_dim=input_dim).to(device)
            logit, _ = model(x)
            prob = float(torch.sigmoid(logit).item())
            probs.append(prob)
            print(f"[{i+1}/{len(df)}] {row['slide_id']} prob_FH={prob:.6f}", flush=True)
    df["prob_FH"] = probs
    df["pred"] = (df["prob_FH"].values >= args.threshold).astype(int)
    pred_path = outdir / f"{args.dataset_name}_ABMIL_final_predictions.csv"
    metrics_path = outdir / f"{args.dataset_name}_ABMIL_final_metrics.csv"
    ci_path = outdir / f"{args.dataset_name}_ABMIL_final_metrics_bootstrap_ci.csv"
    df.to_csv(pred_path, index=False)
    metrics = binary_metrics(df["label"].values, df["prob_FH"].values, threshold=args.threshold)
    pd.DataFrame([metrics]).to_csv(metrics_path, index=False)
    ci = bootstrap_ci(df["label"].values, df["prob_FH"].values, threshold=args.threshold, n_bootstrap=args.bootstrap, seed=args.seed)
    pd.DataFrame([ci]).to_csv(ci_path, index=False)
    print("Saved predictions:", pred_path)
    print("Saved metrics:", metrics_path)
    print("Saved bootstrap CI:", ci_path)
    print(pd.DataFrame([metrics]))


if __name__ == "__main__":
    main()

#!/usr/bin/env python
import argparse
import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

import pandas as pd
import torch

from abmil.model import GatedABMIL
from abmil.io import load_features, torch_load


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


def collect_feature_files(feature_path=None, feature_dir=None):
    files = []

    if feature_path is not None:
        p = Path(feature_path)
        if not p.exists():
            raise FileNotFoundError(f"Feature file not found: {p}")
        files.append(p)

    if feature_dir is not None:
        d = Path(feature_dir)
        if not d.exists():
            raise FileNotFoundError(f"Feature directory not found: {d}")
        files.extend(sorted(d.glob("*.pt")))

    if not files:
        raise RuntimeError("No feature files were provided. Use --feature_path or --feature_dir.")

    return files


def main():
    parser = argparse.ArgumentParser(description="Run ABMIL prediction on unlabeled .pt feature files.")
    parser.add_argument("--feature_path", default=None, help="Path to one slide-level .pt feature file.")
    parser.add_argument("--feature_dir", default=None, help="Directory containing multiple slide-level .pt feature files.")
    parser.add_argument("--model", required=True)
    parser.add_argument("--out_csv", required=True)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()

    device = args.device if (args.device == "cpu" or torch.cuda.is_available()) else "cpu"
    model, config = load_model(args.model, device)
    input_dim = int(config.get("input_dim", 1024))

    rows = []
    feature_files = collect_feature_files(args.feature_path, args.feature_dir)

    with torch.no_grad():
        for i, f in enumerate(feature_files, start=1):
            x = load_features(f, expected_dim=input_dim).to(device)
            logit, attn = model(x)
            prob = float(torch.sigmoid(logit).item())
            pred = int(prob >= args.threshold)

            rows.append({
                "slide_id": f.stem,
                "feature_path": str(f.resolve()),
                "n_patches": int(x.shape[0]),
                "prob_FH": prob,
                "pred": pred,
                "pred_label_name": "FH" if pred == 1 else "nonFH",
                "threshold": args.threshold,
            })

            print(f"[{i}/{len(feature_files)}] {f.stem} n_patches={x.shape[0]} prob_FH={prob:.6f}")

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print("Saved:", out_csv)


if __name__ == "__main__":
    main()

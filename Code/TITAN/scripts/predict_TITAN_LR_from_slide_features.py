#!/usr/bin/env python3
from pathlib import Path
import argparse
import h5py
import joblib
import numpy as np
import pandas as pd

def load_feature(h5_path):
    with h5py.File(h5_path, "r") as f:
        x = np.asarray(f["features"], dtype=np.float32)
    return x.reshape(-1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True)
    parser.add_argument("--feature_dir", required=True)
    parser.add_argument("--out_csv", required=True)
    parser.add_argument("--label_name", default="unknown")
    args = parser.parse_args()

    model = joblib.load(args.model)
    feature_dir = Path(args.feature_dir)
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    X_list = []

    for p in sorted(feature_dir.glob("*.h5")):
        x = load_feature(p)
        if len(x) != 768:
            raise RuntimeError(f"Unexpected feature dimension {len(x)} in {p}")

        rows.append({
            "sample_id": p.stem,
            "label_name": args.label_name,
            "h5_path": str(p),
            "feature_dim": len(x),
        })
        X_list.append(x)

    if not X_list:
        raise RuntimeError(f"No h5 files found in {feature_dir}")

    X = np.vstack(X_list)
    df = pd.DataFrame(rows)

    prob = model.predict_proba(X)[:, 1]
    pred = (prob >= 0.5).astype(int)

    df["TITAN_LR_probability"] = prob
    df["TITAN_LR_pred_0p5"] = pred
    df["TITAN_LR_pred_label_0p5"] = np.where(pred == 1, "FH", "nonFH")

    df.to_csv(out_csv, index=False)
    print("Saved:", out_csv)

if __name__ == "__main__":
    main()

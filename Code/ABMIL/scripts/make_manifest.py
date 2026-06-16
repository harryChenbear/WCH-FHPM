#!/usr/bin/env python
import argparse
from pathlib import Path
import pandas as pd


def collect(pt_dir, label, label_name):
    pt_dir = Path(pt_dir)
    rows = []
    for f in sorted(pt_dir.glob("*.pt")):
        rows.append({
            "slide_id": f.stem,
            "feature_path": str(f.resolve()),
            "label": int(label),
            "label_name": label_name,
        })
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fh_pt_dir", required=True)
    parser.add_argument("--nonfh_pt_dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    rows = []
    rows.extend(collect(args.fh_pt_dir, 1, "FH"))
    rows.extend(collect(args.nonfh_pt_dir, 0, "nonFH"))
    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError("No .pt files found.")
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print("Saved:", out)
    print("Total:", len(df))
    print(df["label_name"].value_counts())


if __name__ == "__main__":
    main()

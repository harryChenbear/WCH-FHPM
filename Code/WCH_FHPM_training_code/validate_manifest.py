#!/usr/bin/env python3
from pathlib import Path
import argparse
import pandas as pd


def read_manifest(path: str) -> pd.DataFrame:
    suffix = Path(path).suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix == ".tsv":
        return pd.read_csv(path, sep="\t")
    raise ValueError(f"Unsupported file type: {suffix}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest")
    args = parser.parse_args()
    df = read_manifest(args.manifest)

    print("=" * 80)
    print("Manifest:", args.manifest)
    print("Rows:", len(df))
    print("Columns:", list(df.columns))

    missing = {"label", "slide_id"} - set(df.columns)
    if missing:
        raise SystemExit(f"Missing required columns: {sorted(missing)}")

    has_patch = "patch_path" in df.columns
    has_coords = {"wsi_path", "x", "y"}.issubset(df.columns)
    if not (has_patch or has_coords):
        raise SystemExit("Need patch_path, or wsi_path/x/y.")

    print("\nLabel counts:")
    print(pd.to_numeric(df["label"], errors="coerce").value_counts(dropna=False).sort_index())
    print("\nUnique slides:", df["slide_id"].nunique())
    print("Duplicate rows:", int(df.duplicated().sum()))

    if has_patch:
        exists = df["patch_path"].astype(str).map(Path.exists)
        print("Existing patch files:", int(exists.sum()), "/", len(df))
        if not exists.all():
            print(df.loc[~exists, "patch_path"].head(10).to_string(index=False))

    if has_coords:
        exists = df["wsi_path"].astype(str).map(Path.exists)
        print("Existing WSI files:", int(exists.sum()), "/", len(df))
        if not exists.all():
            print(df.loc[~exists, "wsi_path"].head(10).to_string(index=False))


if __name__ == "__main__":
    main()

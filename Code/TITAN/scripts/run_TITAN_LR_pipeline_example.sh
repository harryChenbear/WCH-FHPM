#!/usr/bin/env bash
set -euo pipefail
WSI_DIR="${1:-}"
JOB_DIR="${2:-}"
OUT_CSV="${3:-predictions.csv}"
GPU="${4:-0}"
bash scripts/extract_TITAN_features_with_TRIDENT.sh "$WSI_DIR" "$JOB_DIR" "$GPU"
python scripts/predict_TITAN_LR_from_slide_features.py \
  --model models/TITAN_LR_Dataset1_model.joblib \
  --feature_dir "$JOB_DIR/20x_512px_0px_overlap/slide_features_titan" \
  --out_csv "$OUT_CSV" \
  --label_name unknown

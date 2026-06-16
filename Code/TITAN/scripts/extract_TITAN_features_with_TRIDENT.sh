#!/usr/bin/env bash
set -euo pipefail
WSI_DIR="${1:-}"
JOB_DIR="${2:-}"
GPU="${3:-0}"
TRIDENT_DIR="${TRIDENT_DIR:-/root/software/TRIDENT}"
if [ -z "$WSI_DIR" ] || [ -z "$JOB_DIR" ]; then
  echo "Usage: bash scripts/extract_TITAN_features_with_TRIDENT.sh /path/to/wsi_dir /path/to/job_dir 0"
  exit 1
fi
cd "$TRIDENT_DIR"
python run_batch_of_slides.py \
  --task all \
  --wsi_dir "$WSI_DIR" \
  --job_dir "$JOB_DIR" \
  --slide_encoder titan \
  --mag 20 \
  --patch_size 512 \
  --gpus "$GPU"
echo "TITAN features saved to: $JOB_DIR/20x_512px_0px_overlap/slide_features_titan"

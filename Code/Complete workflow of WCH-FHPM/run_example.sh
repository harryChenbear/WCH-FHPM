#!/usr/bin/env bash
set -e

SLIDE="${1:-/path/to/slide.svs}"
OUT="${2:-./wch_fhpm_result}"

python run_wch_fhpm.py \
  --slide "$SLIDE" \
  --out_dir "$OUT" \
  --device cuda \
  --tile_size 512 \
  --step_size 512 \
  --min_tissue_fraction 0.20 \
  --tumor_threshold 0.20 \
  --fh_cutoff 0.50 \
  --read_batch_size 128 \
  --batch_size 128

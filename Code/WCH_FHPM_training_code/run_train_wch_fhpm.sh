#!/usr/bin/env bash
set -euo pipefail

# Four GPUs; effective global batch size = 16 × 4 × 6 = 384
CUDA_VISIBLE_DEVICES=0,1,2,3 \
torchrun --standalone --nproc_per_node=4 train_wch_fhpm.py \
  --train_manifest /path/to/fh_train_manifest.parquet \
  --val_manifest /path/to/fh_val_manifest.parquet \
  --output_dir ./outputs/wch_fhpm \
  --init imagenet \
  --img_size 512 \
  --tile_size 512 \
  --epochs 100 \
  --batch_size 16 \
  --grad_accum_steps 6 \
  --num_workers 8 \
  --learning_rate 0.0001 \
  --weight_decay 0.0 \
  --min_learning_rate 0.000001 \
  --seed 2026

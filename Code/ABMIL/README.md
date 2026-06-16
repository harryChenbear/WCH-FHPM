# ABMIL model for FH-dRCC slide-level prediction

This repository contains the final ABMIL slide-level classifier for predicting FH-deficient renal cell carcinoma (FH-dRCC) from pre-extracted patch-level features.

## Important note

This repository does not directly process raw WSI files such as `.svs`, `.ndpi`, `.tif`, or `.tiff`.

The ABMIL model requires one pre-extracted `.pt` feature file per slide.

Each `.pt` file must contain a PyTorch tensor with shape:

```text
N_patches x 1024
```

where each row is a 1024-dimensional patch-level feature vector.

The model was trained using CLAM-style ResNet50 patch-level features. To obtain reliable predictions, users should generate compatible 1024-dimensional patch-level features using the same or equivalent feature extraction setting.

## Repository contents

```text
abmil/model.py                              ABMIL model architecture
abmil/io.py                                 Feature loading utilities
abmil/metrics.py                            Metric calculation utilities
scripts/check_features.py                   Check whether .pt feature files are compatible
scripts/make_manifest.py                    Create a manifest from feature folders
scripts/evaluate_external.py                Run ABMIL prediction and calculate metrics
models/WCH_ABMIL_FHPM_Dataset1_final.pt     Final ABMIL model checkpoint
models/final_model_config.json              Model configuration
```

## Model input

For each slide, prepare one `.pt` file:

```text
features_FH/pt_files/slide_1.pt
features_nonFH/pt_files/slide_2.pt
```

Each `.pt` file should be loadable by `torch.load()` and should contain either:

```python
torch.Tensor  # shape: [N_patches, 1024]
```

or a dictionary containing:

```python
{"features": torch.Tensor}
```

## How to generate compatible patch-level features

The original workflow used a CLAM-style feature extraction pipeline.

Typical WSI patching settings were:

```text
patch_size = 512
step_size  = 512
patch_level = 0
feature extraction input size = 224 x 224
normalization = ImageNet mean/std
feature dimension = 1024
```

Example CLAM commands:

```bash
python create_patches_fp.py \
  --source /path/to/slides \
  --save_dir /path/to/patches \
  --patch_size 512 \
  --step_size 512 \
  --seg \
  --patch \
  --stitch

python extract_features_fp.py \
  --data_h5_dir /path/to/patches \
  --data_slide_dir /path/to/slides \
  --csv_path /path/to/patches/process_list_autogen.csv \
  --feat_dir /path/to/features \
  --batch_size 128 \
  --slide_ext .svs
```

After feature extraction, the expected output is:

```text
/path/to/features/pt_files/slide_1.pt
/path/to/features/pt_files/slide_2.pt
```

## Check feature files

Before prediction, check that the `.pt` files have the correct shape:

```bash
python scripts/check_features.py \
  --pt_dir /path/to/features/pt_files
```

Expected output example:

```text
slide_1.pt (31888, 1024) torch.float32
```

## Create manifest

```bash
python scripts/make_manifest.py \
  --fh_pt_dir /path/to/features_FH/pt_files \
  --nonfh_pt_dir /path/to/features_nonFH/pt_files \
  --output /path/to/ABMIL_manifest.csv
```

The manifest contains:

```text
slide_id,feature_path,label,label_name
```

where:

```text
label = 1 for FH-dRCC
label = 0 for non-FH-dRCC
```

## Run prediction

```bash
python scripts/evaluate_external.py \
  --manifest /path/to/ABMIL_manifest.csv \
  --model models/WCH_ABMIL_FHPM_Dataset1_final.pt \
  --outdir /path/to/ABMIL_output \
  --dataset_name External \
  --threshold 0.5 \
  --bootstrap 1000
```

## Output files

```text
External_ABMIL_final_predictions.csv
External_ABMIL_final_metrics.csv
External_ABMIL_final_metrics_bootstrap_ci.csv
```

The most important prediction column is:

```text
prob_FH = predicted probability of FH-dRCC
```


## Predict unlabeled slides

If you only want to predict new slides and do not have known FH/non-FH labels, use:

```bash
python scripts/predict_unlabeled.py \
  --feature_dir /path/to/pt_files \
  --model models/WCH_ABMIL_FHPM_Dataset1_final.pt \
  --out_csv /path/to/ABMIL_predictions.csv \
  --threshold 0.5
```

For one slide:

```bash
python scripts/predict_unlabeled.py \
  --feature_path /path/to/slide_1.pt \
  --model models/WCH_ABMIL_FHPM_Dataset1_final.pt \
  --out_csv /path/to/ABMIL_prediction_one_slide.csv \
  --threshold 0.5
```

This produces per-slide prediction probabilities without requiring labels:

```text
slide_id,feature_path,n_patches,prob_FH,pred,pred_label_name,threshold
```

Use `scripts/evaluate_external.py` only when ground-truth labels are available and you want AUROC, sensitivity, specificity, PPV, and NPV.

## What is not included

This repository intentionally does not include raw WSI files, patch coordinate `.h5` files, patch-level feature `.pt` files, training scripts, validation results, or per-sample prediction tables.

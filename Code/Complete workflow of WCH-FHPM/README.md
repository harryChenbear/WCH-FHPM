# Complete Workflow of WCH-FHPM

This repository provides a complete inference workflow for **WCH-FHPM**, a whole-slide image analysis model for FH-deficient renal cell carcinoma (FH-dRCC) research.

The workflow starts from a whole-slide image (WSI) and generates:

1. tissue-filtered image patches by background removal,
2. tumor-associated patch detection using the tumor detector,
3. FH-dRCC probability prediction using the FH predictor,
4. slide-level FH-dRCC probability,
5. patch-level FH-dRCC probability heatmap.

This repository is intended for **research use only**. It is not intended for clinical diagnosis, treatment selection, or patient management.

---

## Workflow overview

```text
Whole-slide image
        |
        v
Tissue filtering / background removal
        |
        v
Tumor-associated patch detection
        |
        v
FH-dRCC probability prediction on tumor-associated patches
        |
        v
Slide-level aggregation
        |
        v
CSV outputs and FH-dRCC probability heatmap
```

The final slide-level FH-dRCC probability is calculated as the mean FH-dRCC probability across tumor-associated patches.

---

## Repository structure

```text
Complete workflow of WCH-FHPM/
├── README.md
├── requirements.txt
├── run_wch_fhpm.py
├── run_example.sh
├── model_config.json
├── .gitignore
└── models/
    ├── WCH_TumorDetector_ViTBase512.pth
    └── WCH_FHPM_ViTBase512.pth
```

The two model weight files are required for inference:

```text
models/WCH_TumorDetector_ViTBase512.pth
models/WCH_FHPM_ViTBase512.pth
```

The required model weight files, **WCH_FHPM_ViTBase512.pth** and **WCH_TumorDetector_ViTBase512.pth**, can be downloaded from the GitHub release page:

```text
https://github.com/harryChenbear/WCH-FHPM/releases/tag/1.0
```

After downloading, place them under the `models/` folder as follows:

```text
Complete workflow of WCH-FHPM/
└── models/
    ├── WCH_TumorDetector_ViTBase512.pth
    └── WCH_FHPM_ViTBase512.pth
```

---

## Model files

### Tumor detector

```text
models/WCH_TumorDetector_ViTBase512.pth
```

This model identifies tumor-associated tissue patches from tissue-filtered patches.

Main metadata:

```text
Model architecture: vit_base_patch16_224
Input size: 512 x 512
Number of classes: 2
Classes:
  normal_tissue: 0
  tumor_associated_tissue: 1
Recommended tumor-associated threshold: 0.20
```

The default `--tumor_threshold` is **0.20**. In the two-class tumor detector, the tumor-associated probability and normal-tissue probability sum to 1. Therefore, a tumor-associated probability threshold of 0.20 corresponds to excluding patches with a predicted normal-tissue probability greater than 0.80. This conservative threshold removes patches that are highly likely to represent normal renal tissue while retaining ambiguous, mixed, stromal, or tumor-adjacent patches for downstream WCH-FHPM analysis.

### FH predictor

```text
models/WCH_FHPM_ViTBase512.pth
```

This model predicts FH-dRCC probability on tumor-associated H&E image patches.

Main metadata:

```text
Model architecture: vit_base_patch16_224
Input size: 512 x 512
Number of classes: 2
Classes:
  nonFH: 0
  FH: 1
Default slide-level cutoff: 0.50
Slide-level aggregation: mean patch probability
```

---

## Installation

### 1. Install Python dependencies

```bash
pip install -r requirements.txt
```

### 2. Install OpenSlide system libraries

OpenSlide is required for reading whole-slide images.

For Ubuntu:

```bash
apt-get update
apt-get install -y openslide-tools libopenslide0
```

---

## Required Python packages

The main Python dependencies are:

```text
torch
torchvision
timm
openslide-python
pandas
numpy
Pillow
```

These packages are listed in `requirements.txt`.

---

## Supported input formats

The workflow supports common whole-slide image formats, including:

```text
.svs
.ndpi
.tif
.tiff
.mrxs
.scn
.kfb
```

The current implementation is recommended for 40× H&E-stained WSIs scanned at approximately 0.225 μm/pixel. Lower-resolution images can be processed, but performance may be less balanced.

---

## Run inference on one slide

```bash
python run_wch_fhpm.py \
  --slide /path/to/slide.svs \
  --out_dir ./results_one_slide \
  --device cuda
```

To run on CPU:

```bash
python run_wch_fhpm.py \
  --slide /path/to/slide.svs \
  --out_dir ./results_one_slide \
  --device cpu
```

---

## Run inference on a folder of slides

```bash
python run_wch_fhpm.py \
  --slide_dir /path/to/slide_folder \
  --out_dir ./results_slide_folder \
  --device cuda
```

The script will recursively search for supported whole-slide image files in the input folder.

---

## Run with custom parameters

```bash
python run_wch_fhpm.py \
  --slide /path/to/slide.svs \
  --out_dir ./results_one_slide \
  --device cuda \
  --tile_size 512 \
  --step_size 512 \
  --min_tissue_fraction 0.20 \
  --tumor_threshold 0.20 \
  --fh_cutoff 0.50 \
  --read_batch_size 128 \
  --batch_size 128
```

Main parameters:

```text
--tile_size              Patch size. Default: 512
--step_size             Step size for patch extraction. Default: 512
--min_tissue_fraction   Minimum tissue fraction for retaining a patch. Default: 0.20
--tumor_threshold       Threshold for tumor-associated patch selection. Default: 0.20
                         A value of 0.20 corresponds to excluding patches with
                         predicted normal-tissue probability > 0.80.
--fh_cutoff             Cutoff for slide-level FH-dRCC classification. Default: 0.50
--read_batch_size       Number of tissue patches read before model inference. Default: 128
--batch_size            Model inference batch size. Default: 128
--device                cuda or cpu. Default: cuda
```

---

## Example script

You can also use the example shell script:

```bash
bash run_example.sh /path/to/slide.svs ./wch_fhpm_result
```

---

## Output files

For each slide, the workflow creates one output subfolder containing:

```text
patch_predictions.csv
slide_result.csv
FH_probability_heatmap.png
```

A global summary file is also saved in the main output directory:

```text
all_slide_results.csv
```

The workflow also saves the runtime configuration:

```text
run_config.json
model_config_used.json
```

---

## Output description

### `patch_predictions.csv`

Patch-level prediction results. Main columns include:

```text
slide_name
x
y
tissue_fraction
tumor_probability
is_tumor_associated
FH_probability
```

### `slide_result.csv`

Slide-level prediction result for one slide. Main columns include:

```text
slide_name
status
n_grid_patches
n_tissue_patches
n_tumor_patches
tumor_fraction_among_tissue
final_fh_probability
final_prediction
final_pred_label_0.5
elapsed_minutes
```

### `all_slide_results.csv`

Summary of slide-level prediction results for all processed slides.

### `FH_probability_heatmap.png`

Patch-level FH-dRCC probability heatmap generated from tumor-associated patches.

---

## Interpretation

The default slide-level interpretation is:

```text
final_fh_probability >= 0.5  -> FH-dRCC high-risk prediction
final_fh_probability < 0.5   -> non-FH-dRCC low-risk prediction
```

The final probability is calculated as:

```text
mean FH-dRCC probability over tumor-associated patches
```

High-risk predictions should be interpreted together with histomorphology and, when appropriate, confirmatory FH/2SC immunohistochemistry and molecular testing.

---

## Notes on model weights

The model files are inference checkpoints. They contain model parameters and minimal inference-related metadata, such as model architecture, input size, preprocessing information, label mapping, and recommended inference thresholds.

The released inference package does not include training slides, patient data, training manifests, or raw dataset files.

---

## Research use only

This repository is provided for research use only. The outputs of this workflow should not be used as the sole basis for clinical diagnosis, treatment selection, or patient management.

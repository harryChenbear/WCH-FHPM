# WCH-FHPM reference training code

This folder provides a clean reference implementation for training:

1. the WCH tumor detector; and
2. WCH-FHPM, the FH-dRCC patch classifier.

The implementation is aligned with the released inference package and the model description used in the manuscript:

- ViT-Base (`vit_base_patch16_224`)
- input size: 512 × 512 RGB
- two output classes
- ImageNet-pretrained initialization
- AdamW optimizer
- initial learning rate: 1e-4
- weight decay: 0
- cosine-annealing learning-rate schedule
- up to 100 epochs
- effective global batch size: 384 in the provided four-GPU examples

## Scope statement

This is a reference implementation reconstructed to match the reported architecture, preprocessing, optimization settings, label definitions, and released checkpoint format. It is not claimed to reproduce the historical training run bit-for-bit, because the exact historical random state, data ordering, hardware kernels, and all intermediate training files are not included.

## Files

```text
train_common.py
train_tumor_detector.py
train_wch_fhpm.py
validate_manifest.py
run_train_tumor_detector.sh
run_train_wch_fhpm.sh
requirements_train.txt
examples/
```

## Manifest formats

Training and validation manifests must be separated at the slide or patient level. Never place patches from the same slide in both sets.

### Cached patch images

```text
patch_path,label,slide_id
/path/to/patch_0001.jpg,1,slide_001
/path/to/patch_0002.jpg,0,slide_002
```

### Direct WSI coordinate reading

```text
wsi_path,x,y,level,label,slide_id
/path/to/slide.svs,1024,2048,0,1,slide_001
```

The `level` column is optional and defaults to 0. Coordinates are interpreted in the OpenSlide level-0 reference frame.

## Labels

### Tumor detector

```text
0 = normal renal tissue
1 = tumor-associated tissue
```

Tumor-associated tissue may include tumor cells and surrounding stroma.

### WCH-FHPM

```text
0 = non-FH-dRCC
1 = FH-dRCC
```

Only tumor-associated patches should be included in the WCH-FHPM manifests. Patch labels correspond to the molecularly defined diagnosis of the source slide.

## Installation

```bash
pip install -r requirements_train.txt
```

For coordinate-based WSI reading on Ubuntu:

```bash
apt-get update
apt-get install -y openslide-tools libopenslide0
```

## Validate manifests

```bash
python validate_manifest.py /path/to/train_manifest.parquet
python validate_manifest.py /path/to/val_manifest.parquet
```

## Four-GPU training

Edit the paths in the shell scripts and run:

```bash
bash run_train_tumor_detector.sh
bash run_train_wch_fhpm.sh
```

The examples use:

```text
per-GPU batch size = 16
number of GPUs = 4
gradient accumulation steps = 6
effective global batch size = 16 × 4 × 6 = 384
```

Adjust per-GPU batch size and gradient accumulation according to GPU memory while retaining the desired effective global batch size.

## Outputs

The tumor-detector output directory contains:

```text
training_config.json
training_history.csv
last_training_state.pt
best_training_state.pt
WCH_TumorDetector_ViTBase512.pth
```

The WCH-FHPM output directory contains the analogous files and:

```text
WCH_FHPM_ViTBase512.pth
```

The exported `.pth` files use the metadata layout expected by the released `run_wch_fhpm.py` inference workflow.

## Resume training

```bash
torchrun --standalone --nproc_per_node=4 train_wch_fhpm.py \
  ... \
  --resume ./outputs/wch_fhpm/last_training_state.pt
```

## Optional class weights

```bash
--class_weights 1.0 2.0
```

The default is unweighted cross-entropy. Any use of class weighting should be reported explicitly.

## Reproducibility notes

For the closest online/offline agreement:

- use identical model checkpoints;
- use identical preprocessing metadata;
- keep software versions consistent;
- record input WSI hashes;
- use fixed seeds;
- use the same patch manifests and slide-level splits.

# WCH-FHPM
A deep learning pipeline for FH-dRCC diagnosis in whole-slide histopathology images using Vision Transformers (ViT).

## Requirements
- Python 3.8+
- OpenSlide (Windows/Linux)
- CUDA-enabled GPU (recommended)

## Data Preparation
### Directory Structure
+data
 -patient_01.svs
 -patient_02.svs
...
+annotations
 -patient_01.geojson
 -patient_02.geojson
...

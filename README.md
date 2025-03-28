# WCH-FHPM
A deep learning pipeline for FH-dRCC diagnosis in whole-slide histopathology images using Vision Transformers (ViT).

## Requirements
- Python 3.8+
- OpenSlide (Windows/Linux)
- CUDA-enabled GPU (recommended)

## Data Preparation
📁 data/ <i>(Whole Slide Images)</i>

├── 🔬 patient_01.svs

├── 🔬 patient_02.svs

└── ⋮

📁 annotations/ <i>(Pathology annotations)</i>

├── 🏷️ patient_01.geojson

├── 🏷️ patient_02.geojson

└── ⋮

## Data Requirements
- SVS whole-slide images (40x magnification)
- GeoJSON annotations containing normal kidney regions

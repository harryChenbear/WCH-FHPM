# WCH-FHPM
A deep learning pipeline for FH-dRCC diagnosis in whole-slide histopathology images using Vision Transformers (ViT).

## Requirements
- Python 3.8+
- OpenSlide (Windows/Linux)
- CUDA-enabled GPU (recommended)

## Installation
```bash
# Clone repository
git clone https://github.com/yourusername/rcc-detection.git
cd rcc-detection

# Install dependencies
conda create -n rcc python=3.8
conda activate rcc
pip install -r requirements.txt

# Install OpenSlide
## Windows
conda install -c conda-forge openslide
## Linux
sudo apt-get install openslide-tools

## Data Preparation
### Directory Structure
data/
├── patient_01.svs
├── patient_02.svs
...
annotations/
├── patient_01.geojson
├── patient_02.geojson
...

# WCH-FHPM: FH-RCC Deep Learning Prediction Model

A specialized deep learning framework utilizing Vision Transformers (ViT) to identify Fumarate Hydratase-deficient Renal Cell Carcinoma (FH-dRCC) from Whole Slide Images (WSI).

## Overview

This project implements a two-stage classification pipeline based on pathological slices and annotations:
1.  **Tissue Screening**: Distinguishes between normal kidney tissue and tumor regions using a binary classifier.
2.  **FH-RCC Diagnosis**: Specifically identifies FH-RCC within the detected tumor regions.

## Key Results

The model's performance and study findings are summarized in the following figures:

### [Figure 1: The overview of the study workflow](Result/Figure%201.png)
![Figure 1](Result/Figure%201.png)
*Five datasets were included in this study for model development, internal testing, multi-center testing, TCGA validation, and proof-of-concept validation. The workflow covers WSI preprocessing, deep learning predictions, pathologist comparisons, interpretability, and web platform validation.*

### [Figure 2: Performance of the WCH-FHPM in multiple datasets](Result/Figure%202.png)
![Figure 2](Result/Figure%202.png)
*ROC curves and confusion matrices showcasing the precision and diagnostic performance of predicting FH-dRCC across Dataset 2, Dataset 3, and Dataset 4.*

### [Figure 3: Three different human-machine strategies and their impact on enhancing pathologists' ability to predict FH-dRCC](Result/Figure%203.png)
![Figure 3](Result/Figure%203.png)
*Illustrates AND, OR, and Modified strategies, demonstrating how WCH-FHPM assists and improves each pathologist's diagnostic performance for FH-dRCC.*

### [Figure 4: Interpretability of WCH-FHPM based on pathological morphology](Result/Figure%204.png)
![Figure 4](Result/Figure%204.png)
*Representative high/low-score patches and heatmaps, with corresponding microscopic morphological review. Highlighting the incidence rates of key morphological features identified by the AI compared to traditional morphological models.*

### [Figure 5: Development of the online platform based on WCH-FHPM and its performance in proof-of-concept validation](Result/Figure%205.png)
![Figure 5](Result/Figure%205.png)
*Illustration of the user-friendly web platform and its practical testing performance using an independent proof-of-concept validation dataset.*

## Technical Implementation

### Requirements
- Python 3.8+
- OpenSlide (Windows/Linux)
- CUDA-enabled GPU (recommended)
- Libraries: `torch`, `torchvision`, `transformers`, `opencv-python`, `openslide-python`, `scikit-learn`

### Data Preparation
The project expects the following directory structure:
```
.
├── data/                  # Whole Slide Images (.svs)
│   ├── patient_01.svs
│   └── ...
├── annotations/           # Annotations (.geojson)
│   ├── patient_01.geojson
│   └── ...
└── WCH-FHPMproject.py     # Main project script
```

### Configuration Note
In `WCH-FHPMproject.py`, you must configure the `OPENSLIDE_PATH` variable to match your local OpenSlide binary installation path before running:
```python
OPENSLIDE_PATH = r"Path\to\your\openslide-bin-4.0.0.6-windows-x64\bin"
```

## Usage
The script includes functions for:
1.  **Training**: `train_normal_tumor_model` and `train_rcc_model`
2.  **Inference**: `slide_level_prediction` aggregates patch scores to predict patient diagnosis.

```python
# Run the main script to perform prediction on test data
python WCH-FHPMproject.py
```

## Online Web Platform

To facilitate clinical use, we have deployed the WCH-FHPM model as an online web service (https://wchrcc.ai4ss.com/).

### Operational Guidelines
1.  **Upload WSI**: Support for SVS, NDPI, or TIF formats (>100MB) of primary renal tumor H&E slides.
2.  **Analysis**: The system processes the image to distinguish tumor tissue and predicts FH-dRCC probability.
3.  **Interpretation**:
    *   **Heatmap**: Visualizes the probability distribution across tumor regions.
    *   **Score**: An overall probability ≥ 0.5 suggests a high risk of FH-dRCC.

*Disclaimer: This tool allows for privacy-preserving analysis (data is not stored) and is intended for research reference, not as a standalone diagnostic confirmation.*


# TITAN+LR Model for FH-dRCC Prediction

This repository provides a trained TITAN+LR downstream classifier for FH-deficient renal cell carcinoma prediction from whole-slide images.

## Included files
- models/TITAN_LR_Dataset1_model.joblib
- scripts/extract_TITAN_features_with_TRIDENT.sh
- scripts/predict_TITAN_LR_from_slide_features.py
- scripts/run_TITAN_LR_pipeline_example.sh
- requirements.txt

## Workflow
1. Use TRIDENT/TITAN to extract slide-level features from WSIs.
2. Apply TITAN_LR_Dataset1_model.joblib to the 768-dimensional TITAN slide features.
3. The output probability is the predicted probability of FH-dRCC.

## TITAN and CONCH weights
MahmoodLab TITAN and CONCH v1.5 pretrained weights are not redistributed in this repository.
Users should request official access through HuggingFace and install/configure TRIDENT before feature extraction.

## Research use only
This model is intended for research use only and is not intended for direct clinical diagnosis.

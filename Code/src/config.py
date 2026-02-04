import os
import torch
import numpy as np

# OpenSlide Configuration
OPENSLIDE_PATH = r"Path\to\your\openslide-bin-4.0.0.6-windows-x64\bin"

# Dataset Configuration
DATA_DIR = ".//data"
ANNOTATION_DIR = ".//annotations"

# Training Configuration
BATCH_SIZE = 384
LEARNING_RATE = 0.0001
NUM_EPOCHS = 100
THRESHOLD = 0.8
PATCH_SIZE = 512

# Model Configuration
MODEL_SAVE_PATH_NORMAL = "normal_tumor_model.pth"
MODEL_SAVE_PATH_RCC = "WCH-FHPM.pth"

def setup_env():
    torch.manual_seed(42)
    np.random.seed(42)

def get_openslide():
    """Import openslide with dll directory setup for Windows"""
    if hasattr(os, 'add_dll_directory') and os.path.exists(OPENSLIDE_PATH):
        with os.add_dll_directory(OPENSLIDE_PATH):
            import openslide
            return openslide
    else:
        try:
            import openslide
            return openslide
        except ImportError:
            print("OpenSlide not found or DLL path incorrect.")
            raise

import os
import cv2
import numpy as np
import json
from torch.utils.data import Dataset
from torchvision import transforms
from sklearn.model_selection import train_test_split
from . import config

# Initialize OpenSlide
openslide = config.get_openslide()

class RCCDataset(Dataset):
    def __init__(self, patches, labels, transform=None):
        self.patches = patches
        self.labels = labels
        self.transform = transform

    def __len__(self):
        return len(self.patches)

    def __getitem__(self, idx):
        patch = self.patches[idx]
        label = self.labels[idx]
        if self.transform:
            patch = self.transform(patch)
        return patch, label

data_transform = transforms.Compose([
    transforms.ToPILImage(),
    transforms.Resize((224, 224)),
    transforms.ToTensor(),
    transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225])
])

def extract_annotations_from_json(annotation_path):
    with open(annotation_path, 'r') as f:
        annotations_json = json.load(f)
    annotations = []
    # Check if 'features' key exists, handle different geojson structures if needed
    if 'features' in annotations_json:
        for annotation in annotations_json['features']:
            if annotation['geometry']['type'] == 'polygon':
                points = np.array(annotation['geometry']['coordinates'], dtype=np.int32)
                annotations.append({'type': 'polygon', 'points': points})
    return annotations

def extract_patches_from_svs(svs_path, annotation_path, patch_size=config.PATCH_SIZE, threshold=config.THRESHOLD):
    slide = openslide.OpenSlide(svs_path)
    width, height = slide.dimensions
    annotations = extract_annotations_from_json(annotation_path)
    mask = np.ones((height, width), dtype=np.uint8)
    for annotation in annotations:
        if annotation['type'] == 'polygon':
            cv2.fillPoly(mask, [annotation['points']], 0)
    
    patches = []
    labels = []
    rcc_patches = []
    rcc_labels = []
    
    for y in range(0, height, patch_size):
        for x in range(0, width, patch_size):
            # Check slide bounds
            if x + patch_size > width or y + patch_size > height:
                continue
                
            patch = slide.read_region((x, y), 0, (patch_size, patch_size))
            patch = np.array(patch)[:, :, :3]
            
            if patch.shape[0] != patch_size or patch.shape[1] != patch_size:
                continue
            
            gray = cv2.cvtColor(patch, cv2.COLOR_RGB2GRAY)
            # Filter white/background patches
            if np.mean(gray) > 240:
                continue
            
            mask_patch = mask[y:y+patch_size, x:x+patch_size]
            normal_ratio = np.sum(mask_patch==0) / (patch_size * patch_size)
            
            if normal_ratio >= threshold:
                patches.append(patch)
                labels.append(0)
            else:
                patches.append(patch)
                labels.append(1)
                if 'rcc' in svs_path.lower(): # Basic check, improved logic might be needed
                    rcc_patches.append(patch)
                    rcc_labels.append(1)
                else:
                    rcc_patches.append(patch)
                    rcc_labels.append(0)
    return patches, labels, rcc_patches, rcc_labels

def prepare_dataset(data_dir, annotation_dir):
    svs_files = sorted([f for f in os.listdir(data_dir) if f.endswith(".svs")])
    annotation_files = sorted([f for f in os.listdir(annotation_dir) if f.endswith(".geojson")])
    
    # Pair them up
    data_pairs = list(zip(svs_files, annotation_files))
    
    train_pairs, val_pairs = train_test_split(data_pairs, test_size=0.2, random_state=42)
    
    def process_pairs(pairs):
        all_patches = []
        all_labels = []
        all_rcc_patches = []
        all_rcc_labels = []
        
        for svs_file, annotation_file in pairs:
            svs_path = os.path.join(data_dir, svs_file)
            annotation_path = os.path.join(annotation_dir, annotation_file)
            patches, labels, rcc_patches, rcc_labels = extract_patches_from_svs(svs_path, annotation_path)
            all_patches.extend(patches)
            all_labels.extend(labels)
            all_rcc_patches.extend(rcc_patches)
            all_rcc_labels.extend(rcc_labels)
        return all_patches, all_labels, all_rcc_patches, all_rcc_labels

    print("Processing Training Pairs...")
    normal_train_patches, normal_train_labels, rcc_train_patches, rcc_train_labels = process_pairs(train_pairs)
    print("Processing Validation Pairs...")
    normal_val_patches, normal_val_labels, rcc_val_patches, rcc_val_labels = process_pairs(val_pairs)

    normal_train_dataset = RCCDataset(normal_train_patches, normal_train_labels, transform=data_transform)
    normal_val_dataset = RCCDataset(normal_val_patches, normal_val_labels, transform=data_transform)
    rcc_train_dataset = RCCDataset(rcc_train_patches, rcc_train_labels, transform=data_transform)
    rcc_val_dataset = RCCDataset(rcc_val_patches, rcc_val_labels, transform=data_transform)
    
    return normal_train_dataset, normal_val_dataset, rcc_train_dataset, rcc_val_dataset

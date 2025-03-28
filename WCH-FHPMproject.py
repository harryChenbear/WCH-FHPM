import os
import cv2
import numpy as np
from torch.utils.data import Dataset, DataLoader
import torch
from torch import nn
from torch.optim import AdamW
from torch.optim.lr_scheduler import CosineAnnealingLR
from torchvision import transforms
from tqdm import tqdm
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from transformers import ViTFeatureExtractor, ViTForImageClassification
import json
from torch.serialization import add_safe_globals

OPENSLIDE_PATH = r"Path\to\your\openslide-bin-4.0.0.6-windows-x64\bin"

if hasattr(os, 'add_dll_directory'):
    with os.add_dll_directory(OPENSLIDE_PATH):
        import openslide
else:
    import openslide

torch.manual_seed(42)
np.random.seed(42)

data_dir = ".//data"
annotation_dir = ".//annotations"

batch_size = 384
learning_rate = 0.0001
num_epochs = 100
threshold = 0.8

patch_size = 512

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
    transforms.ToTensor(),
    transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225])
])

def extract_annotations_from_json(annotation_path):
    with open(annotation_path, 'r') as f:
        annotations_json = json.load(f)
    annotations = []
    for annotation in annotations_json['features']:
        if annotation['geometry']['type'] == 'polygon':
            points = np.array(annotation['geometry']['coordinates'], dtype=np.int32)
            annotations.append({'type': 'polygon','points': points})
    return annotations

def extract_patches_from_svs(svs_path, annotation_path, patch_size=patch_size, threshold=threshold):
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
            patch = slide.read_region((x, y), 0, (patch_size, patch_size))
            patch = np.array(patch)[:, :, :3]
            if patch.shape[0] != patch_size or patch.shape[1] != patch_size:
                continue
            gray = cv2.cvtColor(patch, cv2.COLOR_RGB2GRAY)
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
                if 'rcc' in svs_path:
                    rcc_patches.append(patch)
                    rcc_labels.append(1)
                else:
                    rcc_patches.append(patch)
                    rcc_labels.append(0)
    return patches, labels, rcc_patches, rcc_labels

def prepare_dataset(data_dir, annotation_dir):
    svs_files = [f for f in os.listdir(data_dir) if f.endswith(".svs")]
    annotation_files = [f for f in os.listdir(annotation_dir) if f.endswith(".geojson")]
    all_patches = []
    all_labels = []
    all_rcc_patches = []
    all_rcc_labels = []
    for svs_file, annotation_file in zip(svs_files, annotation_files):
        svs_path = os.path.join(data_dir, svs_file)
        annotation_path = os.path.join(annotation_dir, annotation_file)
        patches, labels, rcc_patches, rcc_labels = extract_patches_from_svs(svs_path, annotation_path)
        all_patches.extend(patches)
        all_labels.extend(labels)
        all_rcc_patches.extend(rcc_patches)
        all_rcc_labels.extend(rcc_labels)
    normal_train_patches, normal_val_patches, normal_train_labels, normal_val_labels = train_test_split(all_patches, all_labels, test_size=0.2, random_state=42)
    rcc_train_patches, rcc_val_patches, rcc_train_labels, rcc_val_labels = train_test_split(all_rcc_patches, all_rcc_labels, test_size=0.2, random_state=42)
    normal_train_dataset = RCCDataset(normal_train_patches, normal_train_labels, transform=data_transform)
    normal_val_dataset = RCCDataset(normal_val_patches, normal_val_labels, transform=data_transform)
    rcc_train_dataset = RCCDataset(rcc_train_patches, rcc_train_labels, transform=data_transform)
    rcc_val_dataset = RCCDataset(rcc_val_patches, rcc_val_labels, transform=data_transform)
    return normal_train_dataset, normal_val_dataset, rcc_train_dataset, rcc_val_dataset

def train_normal_tumor_model(train_dataset, val_dataset, batch_size=batch_size, num_epochs=num_epochs, learning_rate=learning_rate, model_save_path="normal_tumor_model.pth"):
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, num_workers=4, pin_memory=False)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False, num_workers=4, pin_memory=False)
    model = ViTForImageClassification.from_pretrained("google/vit-base-patch16-224", num_labels=2, ignore_mismatched_sizes=True)
    optimizer = AdamW(model.parameters(), lr=learning_rate)
    scheduler = CosineAnnealingLR(optimizer, T_max=num_epochs, eta_min=0)
    criterion = nn.CrossEntropyLoss()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    for epoch in range(num_epochs):
        model.train()
        train_loss = 0.0
        train_corrects = 0
        for images, labels in tqdm(train_loader, desc=f"Epoch {epoch + 1}/{num_epochs}"):
            images = images.to(device)
            labels = labels.to(device)
            optimizer.zero_grad()
            outputs = model(images)
            loss = criterion(outputs.logits, labels)
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * images.size(0)
            train_corrects += torch.sum(torch.argmax(outputs.logits, dim=1) == labels).item()
        scheduler.step()
        model.eval()
        val_loss = 0.0
        val_corrects = 0
        with torch.no_grad():
            for images, labels in tqdm(val_loader, desc="Validation"):
                images = images.to(device)
                labels = labels.to(device)
                outputs = model(images)
                loss = criterion(outputs.logits, labels)
                val_loss += loss.item() * images.size(0)
                val_corrects += torch.sum(torch.argmax(outputs.logits, dim=1) == labels).item()
        train_loss = train_loss / len(train_dataset)
        train_acc = train_corrects / len(train_dataset)
        val_loss = val_loss / len(val_dataset)
        val_acc = val_corrects / len(val_dataset)
        print(f"Epoch {epoch + 1}/{num_epochs}")
        print(f"Train Loss: {train_loss:.4f} Train Acc: {train_acc:.4f}")
        print(f"Val Loss: {val_loss:.4f} Val Acc: {val_acc:.4f}")
    torch.save(model, model_save_path)
    return model

def train_rcc_model(train_dataset, val_dataset, batch_size=batch_size, num_epochs=num_epochs, learning_rate=learning_rate, model_save_path="rcc_model.pth"):
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, num_workers=4, pin_memory=False)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False, num_workers=4, pin_memory=False)
    model = ViTForImageClassification.from_pretrained("google/vit-base-patch16-224", num_labels=2, ignore_mismatched_sizes=True)
    optimizer = AdamW(model.parameters(), lr=learning_rate)
    scheduler = CosineAnnealingLR(optimizer, T_max=num_epochs, eta_min=0)
    criterion = nn.CrossEntropyLoss()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    for epoch in range(num_epochs):
        model.train()
        train_loss = 0.0
        train_corrects = 0
        for images, labels in tqdm(train_loader, desc=f"Epoch {epoch + 1}/{num_epochs}"):
            images = images.to(device)
            labels = labels.to(device)
            optimizer.zero_grad()
            outputs = model(images)
            loss = criterion(outputs.logits, labels)
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * images.size(0)
            train_corrects += torch.sum(torch.argmax(outputs.logits, dim=1) == labels).item()
        scheduler.step()
        model.eval()
        val_loss = 0.0
        val_corrects = 0
        with torch.no_grad():
            for images, labels in tqdm(val_loader, desc="Validation"):
                images = images.to(device)
                labels = labels.to(device)
                outputs = model(images)
                loss = criterion(outputs.logits, labels)
                val_loss += loss.item() * images.size(0)
                val_corrects += torch.sum(torch.argmax(outputs.logits, dim=1) == labels).item()
        train_loss = train_loss / len(train_dataset)
        train_acc = train_corrects / len(train_dataset)
        val_loss = val_loss / len(val_dataset)
        val_acc = val_corrects / len(val_dataset)
        print(f"Epoch {epoch + 1}/{num_epochs}")
        print(f"Train Loss: {train_loss:.4f} Train Acc: {train_acc:.4f}")
        print(f"Val Loss: {val_loss:.4f} Val Acc: {val_acc:.4f}")
    torch.save(model, model_save_path)
    return model

def evaluate_model(model, test_dataset, batch_size=batch_size):
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False, num_workers=4)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    model.eval()
    test_corrects = 0
    with torch.no_grad():
        for images, labels in tqdm(test_loader, desc="Evaluation"):
            images = images.to(device)
            labels = labels.to(device)
            outputs = model(images)
            test_corrects += torch.sum(torch.argmax(outputs.logits, dim=1) == labels).item()
    test_acc = test_corrects / len(test_dataset)
    print(f"Test Acc: {test_acc:.4f}")

def slide_level_prediction(normal_tumor_model, rcc_model, svs_path, annotation_path, patch_size=patch_size):
    patches, _, _, _ = extract_patches_from_svs(svs_path, annotation_path)
    patch_dataset = RCCDataset(patches, [0] * len(patches), transform=data_transform)
    patch_loader = DataLoader(patch_dataset, batch_size=batch_size, shuffle=False, num_workers=4)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    normal_tumor_model.to(device)
    rcc_model.to(device)
    normal_tumor_model.eval()
    rcc_model.eval()
    tumor_probs = []
    rcc_probs = []
    with torch.no_grad():
        for images, _ in tqdm(patch_loader, desc="Normal/Tumor Prediction"):
            images = images.to(device)
            outputs = normal_tumor_model(images)
            probs = torch.softmax(outputs.logits, dim=1).cpu().numpy()
            tumor_probs.extend(probs[:, 1])
        tumor_indices = [i for i, prob in enumerate(tumor_probs) if prob >= threshold]
        tumor_patches = [patches[i] for i in tumor_indices]
        if not tumor_patches:
            return 0.0
        tumor_dataset = RCCDataset(tumor_patches, [0] * len(tumor_patches), transform=data_transform)
        tumor_loader = DataLoader(tumor_dataset, batch_size=batch_size, shuffle=False, num_workers=4)
        for images, _ in tqdm(tumor_loader, desc="RCC Prediction"):
            images = images.to(device)
            outputs = rcc_model(images)
            probs = torch.softmax(outputs.logits, dim=1).cpu().numpy()
            rcc_probs.extend(probs[:, 1])
    print("RCC Patch Probabilities:")
    for prob in rcc_probs:
        print(f"{prob:.4f}")
    slide_rcc_prob = np.mean(rcc_probs) if rcc_probs else 0.0
    return slide_rcc_prob

def patient_level_prediction(normal_tumor_model, rcc_model, svs_paths, annotation_paths, patch_size=patch_size):
    slide_probs = []
    for svs_path, annotation_path in zip(svs_paths, annotation_paths):
        slide_prob = slide_level_prediction(normal_tumor_model, rcc_model, svs_path, annotation_path)
        slide_probs.append(slide_prob)
    patient_prob = np.mean(slide_probs) if len(slide_probs) > 0 else 0.0
    return patient_prob

def main():
    normal_tumor_sava_model = torch.load(".//normal_tumor_model.pth", weights_only=False)
    rcc_sava_model = torch.load(".//WCH-FHPM.pth", weights_only=False)
    svs_path = ".//data//test.svs"
    annotation_path = ".//data//test.geojson"
    slide_prob = slide_level_prediction(normal_tumor_sava_model, rcc_sava_model, svs_path, annotation_path)
    print(f"Slide RCC Probability: {slide_prob:.4f}")

if __name__ == "__main__":
    main()
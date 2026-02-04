import torch
import numpy as np
from torch.utils.data import DataLoader
from tqdm import tqdm
from . import config
from . import dataset

def evaluate_model(model, test_dataset, batch_size=config.BATCH_SIZE):
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

def slide_level_prediction(normal_tumor_model, rcc_model, svs_path, annotation_path, patch_size=config.PATCH_SIZE):
    patches, _, _, _ = dataset.extract_patches_from_svs(svs_path, annotation_path, patch_size=patch_size, threshold=config.THRESHOLD)
    
    # Check if patches were extracted
    if not patches:
        print(f"No valid patches found for {svs_path}")
        return 0.0

    patch_dataset = dataset.RCCDataset(patches, [0] * len(patches), transform=dataset.data_transform)
    patch_loader = DataLoader(patch_dataset, batch_size=config.BATCH_SIZE, shuffle=False, num_workers=4)
    
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
            
        tumor_indices = [i for i, prob in enumerate(tumor_probs) if prob >= config.THRESHOLD]
        tumor_patches = [patches[i] for i in tumor_indices]
        
        if not tumor_patches:
            return 0.0
            
        tumor_dataset = dataset.RCCDataset(tumor_patches, [0] * len(tumor_patches), transform=dataset.data_transform)
        tumor_loader = DataLoader(tumor_dataset, batch_size=config.BATCH_SIZE, shuffle=False, num_workers=4)
        
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

def patient_level_prediction(normal_tumor_model, rcc_model, svs_paths, annotation_paths, patch_size=config.PATCH_SIZE):
    slide_probs = []
    for svs_path, annotation_path in zip(svs_paths, annotation_paths):
        slide_prob = slide_level_prediction(normal_tumor_model, rcc_model, svs_path, annotation_path, patch_size)
        slide_probs.append(slide_prob)
    patient_prob = np.mean(slide_probs) if len(slide_probs) > 0 else 0.0
    return patient_prob

import torch
from torch import nn
from torch.optim import AdamW
from torch.optim.lr_scheduler import CosineAnnealingLR
from torch.utils.data import DataLoader
from transformers import ViTForImageClassification
from tqdm import tqdm
from . import config

def train_normal_tumor_model(train_dataset, val_dataset, 
                           batch_size=config.BATCH_SIZE, 
                           num_epochs=config.NUM_EPOCHS, 
                           learning_rate=config.LEARNING_RATE, 
                           model_save_path=config.MODEL_SAVE_PATH_NORMAL):
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, num_workers=4, pin_memory=False)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False, num_workers=4, pin_memory=False)
    
    model = ViTForImageClassification.from_pretrained("google/vit-base-patch16-224", num_labels=2, ignore_mismatched_sizes=True)
    optimizer = AdamW(model.parameters(), lr=learning_rate, weight_decay=0)
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

def train_rcc_model(train_dataset, val_dataset, 
                    batch_size=config.BATCH_SIZE, 
                    num_epochs=config.NUM_EPOCHS, 
                    learning_rate=config.LEARNING_RATE, 
                    model_save_path=config.MODEL_SAVE_PATH_RCC):
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, num_workers=4, pin_memory=False)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False, num_workers=4, pin_memory=False)
    
    model = ViTForImageClassification.from_pretrained("google/vit-base-patch16-224", num_labels=2, ignore_mismatched_sizes=True)
    optimizer = AdamW(model.parameters(), lr=learning_rate)
    scheduler = CosineAnnealingLR(optimizer, T_max=num_epochs, eta_min=0)
    criterion = nn.CrossEntropyLoss()
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    
    # ... Training loop similar to above ...
    # Note: Refactoring opportunity to make a generic train function
    # But for now, keeping close to original to avoid breaking custom logic if any.
    # The original had identical loops.
    
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

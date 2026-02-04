import torch
from src import config, inference, dataset, train

def main():
    # Setup environment
    config.setup_env()
    
    # Path configuration - Adjust these relative to where main.py is run
    # If running from Code/ directory:
    model_path_normal = ".//normal_tumor_model.pth"
    model_path_rcc = ".//WCH-FHPM.pth"
    svs_path = ".//data//test.svs"
    annotation_path = ".//data//test.geojson"
    
    # Check if models exist before loading, or train them if needed.
    # For this example, we proceed with the inference logic from the original file.
    
    try:
        normal_tumor_sava_model = torch.load(model_path_normal, weights_only=False)
        rcc_sava_model = torch.load(model_path_rcc, weights_only=False)
        
        slide_prob = inference.slide_level_prediction(
            normal_tumor_sava_model, 
            rcc_sava_model, 
            svs_path, 
            annotation_path
        )
        print(f"Slide RCC Probability: {slide_prob:.4f}")
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Please ensure model files and test data are available.")
        # Example on how to run training if needed:
        # train_normal, val_normal, train_rcc, val_rcc = dataset.prepare_dataset(config.DATA_DIR, config.ANNOTATION_DIR)
        # train.train_normal_tumor_model(train_normal, val_normal)
        # train.train_rcc_model(train_rcc, val_rcc)

if __name__ == "__main__":
    main()

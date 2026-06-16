from pathlib import Path
import torch


def torch_load(path, map_location="cpu"):
    try:
        return torch.load(path, map_location=map_location, weights_only=True)
    except TypeError:
        return torch.load(path, map_location=map_location)


def load_features(path, expected_dim=1024):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Feature file not found: {path}")
    x = torch_load(path, map_location="cpu")
    if isinstance(x, dict):
        if "features" in x:
            x = x["features"]
        else:
            raise ValueError(f"Unknown feature dict keys in {path}: {list(x.keys())}")
    if not torch.is_tensor(x):
        raise TypeError(f"Expected torch.Tensor in {path}, got {type(x)}")
    x = x.float()
    if x.ndim != 2:
        raise ValueError(f"Expected 2D tensor in {path}, got shape {tuple(x.shape)}")
    if expected_dim is not None and x.shape[1] != expected_dim:
        raise ValueError(f"Expected feature dim {expected_dim}, got {x.shape[1]} in {path}")
    return x

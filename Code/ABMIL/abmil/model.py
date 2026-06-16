import torch
import torch.nn as nn


class GatedABMIL(nn.Module):
    def __init__(self, input_dim=1024, hidden_dim=512, attn_dim=256, dropout=0.25):
        super().__init__()
        self.projector = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),
        )
        self.attn_V = nn.Sequential(
            nn.Linear(hidden_dim, attn_dim),
            nn.Tanh(),
        )
        self.attn_U = nn.Sequential(
            nn.Linear(hidden_dim, attn_dim),
            nn.Sigmoid(),
        )
        self.attn_w = nn.Linear(attn_dim, 1)
        self.classifier = nn.Linear(hidden_dim, 1)

    def forward(self, x):
        if x.ndim != 2:
            raise ValueError(f"Expected [N_patches, input_dim], got {tuple(x.shape)}")
        h = self.projector(x)
        a = self.attn_w(self.attn_V(h) * self.attn_U(h))
        a = torch.softmax(a, dim=0)
        bag_embedding = torch.sum(a * h, dim=0, keepdim=True)
        logit = self.classifier(bag_embedding).squeeze()
        return logit, a.squeeze(1)

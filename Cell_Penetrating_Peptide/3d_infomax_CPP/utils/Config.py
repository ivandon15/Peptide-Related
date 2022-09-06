import torch
import json
import torch.nn as nn
from torch import optim


class Config():
    def __init__(self):
        self.embed_dim = 256
        self.hidden_dim = 256
        self.num_layers = 2
        self.batch_first = True
        self.num_classes = 2
        self.batch_size = 32
        self.num_epochs = 100
        self.learning_rate = 1e-4
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        # for transformer
        self.d_ff = 2048
        self.n_heads = 8
        self.dropout = 0.1

        with open("data/aa_embedding.json", "r") as f:
            self.aa_num_embed = json.load(f)

        self.criterion = nn.CrossEntropyLoss()

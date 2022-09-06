import torch.nn as nn
import torch
from utils.Config import Config
from utils.util_function import embedding

config = Config()


class GRUCellModel((nn.Module)):
    def __init__(self, embed_dim, hidden_dim, num_layers, num_classes, batch_first=True):
        super(GRUCellModel, self).__init__()
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers
        self.gru = nn.GRUCell(embed_dim, hidden_dim, num_layers)
        self.fc = nn.Linear(hidden_dim, num_classes)

    def forward(self, x):
        ids = x.tolist()
        hiddens = []
        embed = embedding(x)
        # go through all the batch
        for b in range(x.size(0)):
            hidden = embed[b, 0, :].view(1, -1)
            for i, j in enumerate(ids[b]):
                if i == 0:
                    continue
                # sentence i, token j
                if j != 0:
                    hidden = self.gru(embed[b, i, :].view(1, -1), hidden)
                else:
                    break
            hiddens.append(hidden)

        logits = self.fc(torch.cat(hiddens, 0))
        return logits

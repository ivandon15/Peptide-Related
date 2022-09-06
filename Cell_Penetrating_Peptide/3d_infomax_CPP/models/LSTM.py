import torch.nn as nn
from torch.nn.utils.rnn import pack_padded_sequence
from torch.nn.utils.rnn import pad_packed_sequence
import torch
from utils.Config import Config
from utils.util_function import embedding

config = Config()


class LSTM(nn.Module):
    def __init__(self, embed_dim, hidden_dim, num_layers, num_classes, batch_first=True):
        super(LSTM, self).__init__()
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers
        self.lstm = nn.LSTM(embed_dim, hidden_dim, num_layers, batch_first=batch_first)
        self.fc = nn.Linear(hidden_dim, num_classes)
        self.batch_first = batch_first

    def forward(self, x, x_lens):
        x = embedding(x, config.aa_num_embed, config.embed_dim, config.device)
        x_packed = pack_padded_sequence(x, x_lens, batch_first=self.batch_first, enforce_sorted=False)

        output_packed, hidden = self.lstm(x_packed)  # out: tensor of shape (batch_size, seq_length, hidden_size)
        output_padded, _ = pad_packed_sequence(output_packed, batch_first=True)

        indices = torch.LongTensor(x_lens).unsqueeze(-1).repeat(1, config.hidden_dim).unsqueeze(1)

        out = torch.gather(output_padded.to(torch.device("cpu")), 1, indices - 1, out=None, sparse_grad=False).squeeze()
        out = self.fc(out.to(config.device))
        return out

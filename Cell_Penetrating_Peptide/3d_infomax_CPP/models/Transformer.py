from typing import Optional, Any, Tuple, Union
import torch
import torch.nn as nn
import math
from torch import optim
import copy

from utils.Config import Config
import torch.nn.functional as F

from utils.util_function import embedding

config = Config()


def clones(module, N):
    """Produce N identical layers."""
    return nn.ModuleList([copy.deepcopy(module) for _ in range(N)])


class PositionalEncoding(nn.Module):
    def __init__(self, d_model: int, dropout: float = 0.1, max_len: int = 5000):
        super(PositionalEncoding, self).__init__()
        self.dropout = nn.Dropout(p=dropout)

        # Compute the positional encodings once in log space.
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len).unsqueeze(1).float()
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * -(math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)
        self.register_buffer("pe", pe)

    def forward(self, x: torch.FloatTensor) -> torch.FloatTensor:
        """
        Args:
            x: `embeddings`, shape (batch, max_len, d_model)

        Returns:
            `encoder input`, shape (batch, max_len, d_model)
        """

        x = x + self.pe[:, : x.size(1)].to(config.device)
        return self.dropout(x)


class WarmupOptimizer:
    """Optim wrapper that implements rate.
    Warmup是在ResNet论文中提到的一种学习率预热的方法，它在训练开始的时候先选择使用一个较小的学习率，
    训练了一些 epoches,再修改为预先设置的学习率来进行训练
    """

    def __init__(self, base_optimizer: optim.Optimizer, d_model: int, scale_factor: float, warmup_steps: int):
        self.base_optimizer = base_optimizer
        self.warmup_steps = warmup_steps
        self.scale_factor = scale_factor
        self.d_model = d_model
        self._step = 0
        self._rate = 0

    def step(self):
        """Update parameters and rate"""
        self._step += 1
        # 调用底下的rate函数
        self._rate = self.rate()
        for p in self.base_optimizer.param_groups:
            p["lr"] = self._rate
        self.base_optimizer.step()

    def zero_grad(self):
        self.base_optimizer.zero_grad()

    def rate(self, step: Optional[int] = None) -> float:
        """Implement `lrate` above"""
        if step is None:
            step = self._step
        return self.scale_factor * self.d_model ** (-0.5) * min(step ** (-0.5), step * self.warmup_steps ** (-1.5))


class LayerNorm(nn.Module):
    def __init__(self, features: int, eps: float = 1e-6):
        # features = d_model
        super(LayerNorm, self).__init__()
        self.a = nn.Parameter(torch.ones(features).to(config.device))
        self.b = nn.Parameter(torch.zeros(features).to(config.device))
        self.eps = eps

    def forward(self, x: torch.FloatTensor) -> torch.FloatTensor:
        # 对embedding做mean，keepdim就是保持维度，让1维的不要squeeze
        mean = x.mean(-1, keepdim=True)
        std = x.std(-1, keepdim=True)
        return self.a * (x - mean) / (std + self.eps) + self.b


class ScaledDotProductAttention(nn.Module):
    def __init__(self):
        super(ScaledDotProductAttention, self).__init__()

    def forward(self, query: torch.FloatTensor, key: torch.FloatTensor, value: torch.FloatTensor,
                mask: Optional[torch.ByteTensor] = None, dropout: Optional[nn.Dropout] = None) -> Tuple[
        torch.Tensor, Any]:
        """
        Args:
            `query`: shape (batch_size, n_heads, max_len, d_q)
            `key`: shape (batch_size, n_heads, max_len, d_k)
            `value`: shape (batch_size, n_heads, max_len, d_v)
            `mask`: shape (batch_size, 1, 1, max_len)
            `dropout`: nn.Dropout
        Returns:
            `weighted value`: shape (batch_size, n_heads, max_len, d_v)
            `weight matrix`: shape (batch_size, n_heads, max_len, max_len)
        """
        d_k = query.size(-1)  # d_k = d_model / n_heads
        scores = torch.matmul(query, key.transpose(-2, -1)) / math.sqrt(d_k)  # B*H*L*L
        if mask is not None:
            # 对于mask==0的地方，通过masked_fil，将0变成-1e9
            # 因为attention它会对除了自己的每一个去计算一个score: L*L，如果是pad的地方
            # 计算完成了之后再去做mask
            scores = scores.masked_fill(mask.eq(0), -1e9)
        # p_attn就是乘上v的那个权重（attn值），做softmax应该导数第一个维度和第二个都可？
        p_attn = F.softmax(scores, dim=-1)  # B*H*L*L
        if dropout is not None:
            p_attn = dropout(p_attn)

        # 最后乘上对应的value值
        return torch.matmul(p_attn, value), p_attn


class MultiHeadAttention(nn.Module):
    def __init__(self, n_heads: int, d_model: int, dropout: float = 0.1):
        super(MultiHeadAttention, self).__init__()
        assert d_model % n_heads == 0
        # We assume d_v always equals d_k
        self.d_k = d_model // n_heads
        self.h = n_heads
        self.linears = clones(nn.Linear(d_model, d_model), 4)
        self.sdpa = ScaledDotProductAttention()
        self.attn = None
        self.dropout = nn.Dropout(p=dropout)

    def forward(self, query: torch.FloatTensor, key: torch.FloatTensor, value: torch.FloatTensor,
                mask: Optional[torch.ByteTensor] = None) -> torch.FloatTensor:
        """
        Args:
            `query`: shape (batch_size, max_len, d_model)
            `key`: shape (batch_size, max_len, d_model)
            `value`: shape (batch_size, max_len, d_model)
            `mask`: shape (batch_size, max_len)

        Returns:
            shape (batch_size, max_len, d_model)
        """
        if mask is not None:
            # Same mask applied to all h heads. B*1*1*L
            mask = mask.unsqueeze(1).unsqueeze(1)
        batch_size = query.size(0)
        # 1) Do all the linear projections in batch from d_model => h x d_k
        # 对进来的 q k v 都转换成 b -1 h d_k
        query, key, value = [l(x).view(batch_size, -1, self.h, self.d_k).transpose(1, 2) for l, x in
                             zip(self.linears, (query, key, value))]

        # 2) Apply attention on all the projected vectors in batch.
        # x: B x H x L x D_v
        x, self.attn = self.sdpa(query, key, value, mask=mask, dropout=self.dropout)

        # 3) "Concat" using a view and apply a final linear.
        # 原文好像是相加，但这里直接顺序拼接了
        # contiguous https://zhuanlan.zhihu.com/p/64551412，主要是为了在transpose之后还能使用view
        x = x.transpose(1, 2).contiguous().view(batch_size, -1, self.h * self.d_k)
        return self.linears[-1](x)


class FeedForward(nn.Module):
    def __init__(self, d_model: int, d_ff: int, dropout: float = 0.1):
        super(FeedForward, self).__init__()
        self.w_1 = nn.Linear(d_model, d_ff)
        self.w_2 = nn.Linear(d_ff, d_model)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.FloatTensor) -> torch.FloatTensor:
        """
        Args:
            `x`: shape (batch_size, max_len, d_model)
        Returns:
            same shape as input x
        """
        return self.w_2(self.dropout(F.relu(self.w_1(x))))


class EncoderLayer(nn.Module):
    """Encoder is made up of self-attn and feed forward"""

    def __init__(self, size: int, self_attn: MultiHeadAttention, feed_forward: FeedForward, dropout: float):
        super(EncoderLayer, self).__init__()
        self.self_attn = self_attn
        self.feed_forward = feed_forward
        # Multihead之后有一次 residual， ff后有一次
        self.sublayer = clones(SublayerConnection(size, dropout), 2)
        self.size = size

    def forward(self, x: torch.FloatTensor, mask: torch.ByteTensor) -> torch.FloatTensor:
        x = self.sublayer[0](x, lambda x: self.self_attn(x, x, x, mask))
        return self.sublayer[1](x, self.feed_forward)


class SublayerConnection(nn.Module):
    """
    A residual connection followed by a layer norm.
    Note for code simplicity the norm is first as opposed to last.
    """

    def __init__(self, size: int, dropout: float):
        super(SublayerConnection, self).__init__()
        self.norm = LayerNorm(size)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.FloatTensor, sublayer: Union[MultiHeadAttention, FeedForward]) -> torch.FloatTensor:
        """Apply residual connection to any sublayer with the same size."""
        return x + self.dropout(sublayer(self.norm(x)))


class Encoder(nn.Module):
    """Core encoder is a stack of N layers"""

    def __init__(self, layer: EncoderLayer, N: int):
        super(Encoder, self).__init__()
        self.layers = clones(layer, N)
        self.norm = LayerNorm(layer.size)

    def forward(self, x: torch.FloatTensor, mask: torch.ByteTensor) -> torch.FloatTensor:
        """Pass the input (and mask) through each layer in turn."""
        for layer in self.layers:
            x = layer(x, mask)
        return self.norm(x)


class TransformerEncoder(nn.Module):
    """The encoder of transformer

    Args:
        `n_layers`: number of stacked encoder layers
        `d_model`: model dimension
        `d_ff`: hidden dimension of feed forward layer
        `n_heads`: number of heads of self-attention
        `dropout`: dropout rate, default 0.1
    """

    def __init__(self, d_model: int, d_ff: int, n_heads: int = 1, n_layers: int = 1,
                 dropout: float = 0.1):
        super(TransformerEncoder, self).__init__()
        self.multi_headed_attention = MultiHeadAttention(n_heads, d_model, dropout)
        self.feed_forward = FeedForward(d_model, d_ff, dropout)
        self.encoder_layer = EncoderLayer(d_model, self.multi_headed_attention, self.feed_forward, dropout)
        self.encoder = Encoder(self.encoder_layer, n_layers)
        self.pe = PositionalEncoding(d_model)
        self.reset_parameters()

    def reset_parameters(self):
        for p in self.parameters():
            if p.dim() > 1:
                nn.init.xavier_uniform_(p)

    def forward(self, x: torch.FloatTensor, mask: torch.ByteTensor) -> torch.FloatTensor:
        # x是encode然后padded过的，然后需要先embedding一下，再计算positional encoding
        return self.encoder(self.pe(embedding(x, config.aa_num_embed, config.embed_dim, config.device)), mask).to(config.device)


class TransformerModel(nn.Module):
    def __init__(self, d_model: int, d_ff: int, num_classes, n_heads: int = 1, n_layers: int = 1,
                 dropout: float = 0.1):
        super(TransformerModel, self).__init__()
        self.transformer_encoder = TransformerEncoder(d_model, d_ff, n_heads, n_layers, dropout)
        self.relu = nn.ReLU()
        #         self.gelu = nn.GELU()
        self.fc = nn.Linear(d_model, num_classes)

    def forward(self, x: torch.FloatTensor, mask: torch.ByteTensor) -> torch.FloatTensor:
        # 这里的x就是encoded和padded
        logits = self.transformer_encoder(x, mask)
        logits = torch.mean(logits, dim=1)

        logits = self.fc(self.relu(logits))
        #         logits = self.fc(self.gelu(logits))

        return logits.to(config.device)
import torch.nn as nn
from torch.nn.utils.rnn import pad_sequence
import torch
from sklearn.metrics import precision_recall_fscore_support
import numpy as np
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix

from utils.Config import Config

config = Config()
def pad_collate(batch):
    (xx, yy) = zip(*batch)
    x_lens = [len(x) for x in xx]
    xx_pad = pad_sequence(xx, batch_first=True, padding_value=0)

    return xx_pad, x_lens, torch.tensor(yy).type(torch.LongTensor).to(config.device)


def pad_collate_transformer(batch):
    # 把batch中的xx和yy都zip起来，因为输出都是一样的，可以不用pad
    (xx, yy) = zip(*batch)
    xx_pad = pad_sequence(xx, batch_first=True, padding_value=0)
    mask = getmask(xx_pad, config.device)
    return xx_pad, mask, torch.tensor(yy).type(torch.LongTensor).to(config.device)


def embedding(x_pad, aa_num_embed, embed_dim, device):
    input_tensor = torch.randn(x_pad.size(0), x_pad.size(1), embed_dim)
    for b in range(x_pad.size(0)):
        for aa in range(x_pad.size(1)):
            input_tensor[b][aa] = torch.FloatTensor(aa_num_embed[str(x_pad[b][aa].item())])
    return input_tensor.to(device)


def getmask(x_pad, device):
    mask = torch.ones(x_pad.shape)
    mask[(x_pad == 0).squeeze()] = 0
    return mask.to(device)


def get_accuracy_from_logits(logits, labels):
    probs = torch.argmax(logits.squeeze(), dim=1)
    # print("predict: ",probs)
    # print("labels: ",labels)
    acc = (probs.squeeze() == labels).float().mean()
    return acc


def evaluate(model, dev_loader, ep):
    model.eval()

    mean_acc, mean_loss = 0, 0
    count = 0
    predicted = []
    true_label = []
    with torch.no_grad():
        for i, (x_pad, mask, yy) in enumerate(dev_loader):
            # caculate the accuracy and loss for every batch
            logits = model(x_pad.squeeze(), mask.squeeze())
            mean_loss += nn.functional.cross_entropy(logits.squeeze(), yy, reduction='mean')
            mean_acc += get_accuracy_from_logits(logits, yy)
            count += 1

            # record every p r f1
            true_label.extend(yy.tolist())
            probs = torch.argmax(logits.squeeze(), dim=1)
            tmp = (probs.squeeze() == yy).float().tolist()
            predicted.extend(tmp)

    p, r, f, _ = precision_recall_fscore_support(true_label, predicted, pos_label=1, average="macro")
    return mean_acc / count, mean_loss / count, p, r, f


def read_data(file):
    X = []
    with open(file, "r") as f:
        for line in f.readlines():
            X.append(line.strip())

    X = np.array(X)
    return X


def predict(model, test_loader):
    model.eval()

    logits = []
    with torch.no_grad():
        # go through each batch
        for i, (x_pad, mask, yy) in enumerate(test_loader):
            # compute logits
            logit = model(x_pad.squeeze(), mask.squeeze())
            logits.append(logit)

    # concatenate batch
    logits = torch.cat(logits, dim=0)
    probs = torch.argmax(logits.squeeze(), dim=1)

    return probs.cpu().numpy()


def result_info(label, predict):
    tn, fp, fn, tp = confusion_matrix(label, predict).ravel()
    specificity = tn / (tn + fp)
    p, r, f, _ = precision_recall_fscore_support(label, predict, pos_label=1, average="macro")
    return specificity, p, r, f, accuracy_score(label, predict)

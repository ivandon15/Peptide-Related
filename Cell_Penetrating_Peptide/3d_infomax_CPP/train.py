import os
import time

import torch
import torch.nn as nn
from torch import optim
import numpy as np
from torch.utils.data import DataLoader

from dataset.CPPDataset import CPPDataset
from models.Transformer import TransformerModel
from utils.Config import Config
from utils.util_function import get_accuracy_from_logits, evaluate, read_data, pad_collate_transformer


def train(model, optimizer, scheduler, train_loader, dev_loader, max_eps, saving_path):
    print("Starting Training...")
    best_acc = 0
    st = time.time()
    # every epoch
    for ep in range(max_eps):
        running_loss, running_acc = 0, 0

        model.train()
        for i, (x_pad, mask, yy) in enumerate(train_loader):
            # Clear gradients
            optimizer.zero_grad()
            logits = model(x_pad, mask)
            loss = nn.functional.cross_entropy(logits, yy)
            loss.backward()
            optimizer.step()

            running_acc += get_accuracy_from_logits(logits, yy)
            running_loss += loss.item()
        scheduler.step()
        print("iteration: {}, epoch: {}, loss: {}, accuracy: {}; Time taken (s): {}".format(i, ep, running_loss / len(
            train_loader), running_acc / len(train_loader), (time.time() - st)))
        st = time.time()

        dev_acc, dev_loss, p, r, f = evaluate(model, dev_loader, ep)
        print(
            "Epoch {} complete! Development Accuracy: {}; Development Loss: {}, Precision: {}, Recall: {}, F1: {}".format(
                ep, dev_acc, dev_loss, p, r, f))

        if dev_acc > best_acc:
            print("Best development accuracy improved from {} to {}, saving model...".format(best_acc, dev_acc))
            best_acc = dev_acc
            torch.save(model.state_dict(), os.path.join(saving_path, 'model_{}.pt'.format(ep)))


if __name__ == '__main__':
    print("Model starting...")
    y_train_path = "data/y_train.txt"
    y_dev_path = "data/y_dev.txt"
    y_test_path = "data/y_test.txt"
    X_train_path = "data/X_train.txt"
    X_dev_path = "data/X_dev.txt"
    X_test_path = "data/X_test.txt"

    config = Config()
    model_transformer = TransformerModel(config.embed_dim, config.d_ff, config.num_classes, config.n_heads,
                                         config.num_layers, config.dropout).to(config.device)
    optimizer = optim.Adam(model_transformer.parameters(), lr=config.learning_rate)
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=25, gamma=0.5)

    y_train = np.loadtxt(y_train_path)
    y_dev = np.loadtxt(y_dev_path)
    y_test = np.loadtxt(y_test_path)

    X_train, X_dev, X_test = read_data(X_train_path), read_data(X_dev_path), read_data(X_test_path)
    train_dataset = CPPDataset(X_train, y_train, config.device)
    dev_dataset = CPPDataset(X_dev, y_dev, config.device)
    test_dataset = CPPDataset(X_test, y_test, config.device)

    train_loader = DataLoader(dataset=train_dataset, batch_size=config.batch_size, shuffle=True,
                              collate_fn=pad_collate_transformer)
    dev_loader = DataLoader(dataset=dev_dataset, batch_size=config.batch_size, shuffle=True,
                            collate_fn=pad_collate_transformer)
    test_loader = DataLoader(dataset=test_dataset, batch_size=config.batch_size, collate_fn=pad_collate_transformer)

    train(model_transformer, optimizer, scheduler, train_loader, dev_loader, config.num_epochs, "parameters/")

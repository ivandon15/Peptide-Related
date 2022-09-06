import re
import torch.utils.data as data
import torch

aa_codes = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
    'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'K': 'LYS',
    'I': 'ILE', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
    'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
    'T': 'THR', 'V': 'VAL', 'Y': 'TYR', 'W': 'TRP'}

aa_nums = {
    'A': '1', 'C': '2', 'D': '3', 'E': '4', 'F': '5', 'G': '6',
    'H': '7', 'K': '8', 'I': '9', 'L': '10', 'M': '11', 'N': '12',
    'P': '13', 'Q': '14', 'R': '15', 'S': '16', 'T': '17', 'V': '18',
    'Y': '19', 'W': '20', '<PAD>': '0'}


def letter_to_num(string, dict_):
    patt = re.compile('[' + ''.join(dict_.keys()) + ']')
    # lambda m:m.group(0) used for go through the string one by one
    num_string = patt.sub(lambda m: dict_[m.group(0)] + ' ', string)
    num = [int(i) for i in num_string.split()]
    return num


class CPPDataset(data.Dataset):
    def __init__(self, seqs, labels, device):
        self.seqs = seqs
        self.labels = labels
        self.device = device

    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, idx):
        seq_encoding_tensor = torch.tensor(letter_to_num(self.seqs[idx], aa_nums))
        label_tensor = self.labels[idx]
        return seq_encoding_tensor.to(self.device), label_tensor

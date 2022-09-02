import argparse
import os
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from pymol import cmd
import csv
import glob
import math
import pandas as pd


def getSeq(file):
    for record in SeqIO.parse(file, "pdb-atom"):
        seq = str(record.seq)
        return seq


def calculateRMSD(origin_file, predict_file):
    """
    Calculate the RMSD between each two pdb files.
    :return: CA RMSD value
    """
    origin_name = pathlib.Path(origin_file).stem
    predict_name = pathlib.Path(predict_file).stem
    cmd.load(origin_file)
    cmd.load(predict_file)

    cmd.select("originlayer", "model " + origin_name)
    cmd.select("predictlayer", "model " + predict_name)

    output = cmd.align("predictlayer////CA", "originlayer////CA")
    cmd.delete("all")
    return output[0]


def write2csv(file, output):
    """Write the output to csv files"""
    with open(file, "w", newline='') as f:
        wrt = csv.writer(f)
        wrt.writerow(["id1", "id2", "len1", "len2", "lowest1", "lowest2", "seq1_rmsd", "seq2_rmsd", "true_s1_s2_rmsd",
                      "predict_s1_s2_rmsd"])
        wrt.writerows(output)


def lowestEnergy(files):
    """Find the lowest energy conformation among files"""
    lowest = math.inf
    for file in files:
        with open(file, 'r') as f:
            score = float(f.readlines()[-3].split()[1])
            if score < lowest:
                lowest = score
                lowest_file = file
    return lowest_file


def compare(origin_path, pair_list_file, predict_path):
    """
    Cross Comparing, including the original structure and prediction structure with same sequence,
    original structure with mutation sequence pairs and prediction structure (lowest energy) with
    mutation sequence pairs.
    :param origin_path: the original pdb files
    :param pair_list_file: mutation sequence pair files, in txt file, and one pair a row, with ", " split.
    :param predict_path: the predicted pdb files
    :return: a list of list results
    """
    pair_list = []

    with open(pair_list_file, "r") as f:
        try:
            for line in f:
                pair = line.strip().split(",")
                check = pair[1]
                pair_list.append(pair)
        except Exception:
            print("Please make sure your pair list file's format!")
            sys.exit()

    output = []
    # the predict_path may contains several different conformation in one pdbid, need
    # to find the lowest
    for pair in pair_list:
        id1, id2 = pair[0].strip(), pair[1].strip()
        origin_file_1 = os.path.join(origin_path, id1 + ".pdb")
        origin_file_2 = os.path.join(origin_path, id2 + ".pdb")

        first_files = glob.glob(os.path.join(predict_path, id1 + "*.pdb"))
        second_files = glob.glob(os.path.join(predict_path, id2 + "*.pdb"))
        oneline = [id1, id2, len(getSeq(origin_file_1)), len(getSeq(origin_file_2))]

        lowest_file_1 = lowestEnergy(first_files)
        lowest_file_2 = lowestEnergy(second_files)

        # record the file name of the lowest energy conformation file
        oneline.append(pathlib.Path(lowest_file_1).stem)
        oneline.append(pathlib.Path(lowest_file_2).stem)

        # cross comparing
        oneline.append(calculateRMSD(origin_file_1, lowest_file_1))
        oneline.append(calculateRMSD(origin_file_2, lowest_file_2))
        oneline.append(calculateRMSD(origin_file_1, origin_file_2))
        oneline.append(calculateRMSD(lowest_file_1, lowest_file_2))
        output.append(oneline)

    return output


def plotResult(csv_file):
    """
    Visualizing the CSV files
    """
    plt.style.use('seaborn-whitegrid')
    df = pd.read_csv(csv_file)

    id1, id2 = df['id1'].tolist(), df['id2'].tolist()
    seq1_rmsd, seq2_rmsd = df['seq1_rmsd'].tolist(), df['seq2_rmsd'].tolist()
    true_s1_s2_rmsd, predict_s1_s2_rmsd = df['true_s1_s2_rmsd'].tolist(), df['predict_s1_s2_rmsd'].tolist()

    labels = []
    for id in zip(id1, id2):
        labels.append(id[0] + "\n" + id[1])

    plt.figure(figsize=(18, 6), dpi=80)

    # four sets of data
    plt.plot()

    # x axis stands for ids
    x = np.arange(len(labels))
    # bars' width
    width = 0.2

    plt.bar(x - 1.5 * width, seq1_rmsd, width - 0.05, label='row_1_seq_rmsd', color="#5B9BD5")
    plt.bar(x - 0.5 * width, seq2_rmsd, width - 0.05, label='row_2_seq_rmsd', color="#ED7D31")
    plt.bar(x + 0.5 * width, true_s1_s2_rmsd, width - 0.05, label='two_seq_origin_rmsd', color="#A5A5A5")
    plt.bar(x + 1.5 * width, predict_s1_s2_rmsd, width - 0.05, label='two_seq_predict_rmsd', color="#FFC000")
    plt.ylabel('RMSD')
    plt.title("Four RMSD in Sequence Pairs (length <= 20) using Rosetta")
    plt.xticks(x, labels=labels)
    plt.legend()

    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Analyzing the predicted structure")
    parser.add_argument('PDBPATH',
                        metavar='path',
                        type=str,
                        help='the path of the folder of original pdb files')
    parser.add_argument('PAIRLIST',
                        metavar='file',
                        type=str,
                        help='the txt file name of pairlist')
    parser.add_argument('PREPATH',
                        metavar='path',
                        type=str,
                        help='the output path of the folder of predicted files')
    parser.add_argument('CSVFILE',
                        metavar='file',
                        type=str,
                        help='the csv file name')
    args = parser.parse_args()

    if not os.path.isdir(args.PDBPATH) or not os.path.isdir(args.PREPATH):
        print('The path specified does not exist')
        sys.exit()
    elif not os.path.exists(args.PAIRLIST):
        print('The pair list file does not exist')
        sys.exit()
    else:
        print("Comparing the predicted structure and the original ones...")
        output = compare(args.PDBPATH,args.PAIRLIST,args.PREPATH)
        print("DONE in comparing!")

        print("Saving the results in to csv file...")
        write2csv(args.CSVFILE,output)
        print("DONE in saving!")

        print("Visualizing the results...")
        plotResult(args.CSVFILE)
        print("DONE in visualizing!")
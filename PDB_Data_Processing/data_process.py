"""
This file contains several useful functions related to peptide.
Including:
1. Get sequence according to PDB file (with one chain)
2. Get the distance between 2 specific residues by specific atom type
3. Get every pair of residue's C-alpha/C-beta distance in one sequence
4. Get each phi and psi of one sequence
5. Calculate the similarity using Blosum62/90 between 2 sequences
6. Convert the second structure files provided by DSSP to dictionary
7. Extract consecutive second structure
8. Calculate two pdbs' RMSD (C-alpha alignment)
"""
import json
import math
import os
import pathlib
from collections import defaultdict
import warnings
from pymol import cmd
from Bio import SeqIO
from Bio.PDB import PDBParser
from tqdm import tqdm

warnings.filterwarnings("ignore")

# used for calculating similarity
aa_3l = {'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C': 4, 'Q': 5, 'E': 6, 'G': 7, 'H': 8, 'I': 9, 'L': 10,
         'K': 11, 'M': 12, 'F': 13, 'P': 14, 'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19, 'B': 20,
         'Z': 21, 'X': 22, '*': 23}
#            A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
blosum90 = [[5, -2, -2, -3, -1, -1, -1, 0, -2, -2, -2, -1, -2, -3, -1, 1, 0, -4, -3, -1, -2, -1, -1, -6],
            [-2, 6, -1, -3, -5, 1, -1, -3, 0, -4, -3, 2, -2, -4, -3, -1, -2, -4, -3, -3, -2, 0, -2, -6],
            [-2, -1, 7, 1, -4, 0, -1, -1, 0, -4, -4, 0, -3, -4, -3, 0, 0, -5, -3, -4, 4, -1, -2, -6],
            [-3, -3, 1, 7, -5, -1, 1, -2, -2, -5, -5, -1, -4, -5, -3, -1, -2, -6, -4, -5, 4, 0, -2, -6],
            [-1, -5, -4, -5, 9, -4, -6, -4, -5, -2, -2, -4, -2, -3, -4, -2, -2, -4, -4, -2, -4, -5, -3, -6],
            [-1, 1, 0, -1, -4, 7, 2, -3, 1, -4, -3, 1, 0, -4, -2, -1, -1, -3, -3, -3, -1, 4, -1, -6],
            [-1, -1, -1, 1, -6, 2, 6, -3, -1, -4, -4, 0, -3, -5, -2, -1, -1, -5, -4, -3, 0, 4, -2, -6],
            [0, -3, -1, -2, -4, -3, -3, 6, -3, -5, -5, -2, -4, -5, -3, -1, -3, -4, -5, -5, -2, -3, -2, -6],
            [-2, 0, 0, -2, -5, 1, -1, -3, 8, -4, -4, -1, -3, -2, -3, -2, -2, -3, 1, -4, -1, 0, -2, -6],
            [-2, -4, -4, -5, -2, -4, -4, -5, -4, 5, 1, -4, 1, -1, -4, -3, -1, -4, -2, 3, -5, -4, -2, -6],
            [-2, -3, -4, -5, -2, -3, -4, -5, -4, 1, 5, -3, 2, 0, -4, -3, -2, -3, -2, 0, -5, -4, -2, -6],
            [-1, 2, 0, -1, -4, 1, 0, -2, -1, -4, -3, 6, -2, -4, -2, -1, -1, -5, -3, -3, -1, 1, -1, -6],
            [-2, -2, -3, -4, -2, 0, -3, -4, -3, 1, 2, -2, 7, -1, -3, -2, -1, -2, -2, 0, -4, -2, -1, -6],
            [-3, -4, -4, -5, -3, -4, -5, -5, -2, -1, 0, -4, -1, 7, -4, -3, -3, 0, 3, -2, -4, -4, -2, -6],
            [-1, -3, -3, -3, -4, -2, -2, -3, -3, -4, -4, -2, -3, -4, 8, -2, -2, -5, -4, -3, -3, -2, -2, -6],
            [1, -1, 0, -1, -2, -1, -1, -1, -2, -3, -3, -1, -2, -3, -2, 5, 1, -4, -3, -2, 0, -1, -1, -6],
            [0, -2, 0, -2, -2, -1, -1, -3, -2, -1, -2, -1, -1, -3, -2, 1, 6, -4, -2, -1, -1, -1, -1, -6],
            [-4, -4, -5, -6, -4, -3, -5, -4, -3, -4, -3, -5, -2, 0, -5, -4, -4, 11, 2, -3, -6, -4, -3, -6],
            [-3, -3, -3, -4, -4, -3, -4, -5, 1, -2, -2, -3, -2, 3, -4, -3, -2, 2, 8, -3, -4, -3, -2, -6],
            [-1, -3, -4, -5, -2, -3, -3, -5, -4, 3, 0, -3, 0, -2, -3, -2, -1, -3, -3, 5, -4, -3, -2, -6],
            [-2, -2, 4, 4, -4, -1, 0, -2, -1, -5, -5, -1, -4, -4, -3, 0, -1, -6, -4, -4, 4, 0, -2, -6],
            [-1, 0, -1, 0, -5, 4, 4, -3, 0, -4, -4, 1, -2, -4, -2, -1, -1, -4, -3, -3, 0, 4, -1, -6],
            [-1, -2, -2, -2, -3, -1, -2, -2, -2, -2, -2, -1, -1, -2, -2, -1, -1, -3, -2, -2, -2, -1, -2, -6],
            [-6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, 1]]

blosum62 = [[4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -4],
            [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -4],
            [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -4],
            [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4],
            [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4],
            [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1, -4],
            [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4],
            [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -4],
            [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4],
            [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4],
            [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4],
            [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -4],
            [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4],
            [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4],
            [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -4],
            [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -4],
            [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -4],
            [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2, -4],
            [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4],
            [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4],
            [-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1, -4],
            [-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4],
            [0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4],
            [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1]]


def getSeq(file):
    """
    Get sequence according to PDB file
    :param file: pdb file
    :return: a string
    """
    for record in SeqIO.parse(file, "pdb-atom"):
        seq = str(record.seq)
        return seq


def get2ReisdueDistance(path, pdb_chain, pos1, pos2, atom_name):
    """
    Get the distance between specific residue, the distance is between the
    same atoms in two residue. For example, atom_name is CA, the you can get
    the distance between two CA atom for the residues
    :param path: pdbfile's path (not file, but folder path)
    :param pdb_chain: pdbfile's name
    :param pos1: int, residue 1's position
    :param pos2: int
    :param atom_name: CA, CB, SG....
    :return:
    """
    cmd.load(os.path.join(path, pdb_chain + ".pdb"))
    cmd.remove('solvent')
    # select two residues in pymol
    cmd.select('resi1', 'resi ' + str(pos1))
    cmd.select('resi2', 'resi ' + str(pos2))

    atoms_1 = []
    atoms_2 = []
    for atom_c in cmd.get_model("resi1").atom:
        if atom_c.name == atom_name:
            atoms_1.append(atom_c)
    for atom_c in cmd.get_model("resi2").atom:
        if atom_c.name == atom_name:
            atoms_2.append(atom_c)
    cmd.delete('all')

    # prevent the situation that one residue contains two same atoms
    if len(atoms_2) != 1 and len(atoms_1) != 1:
        return None
    distance = getAtomDistance(atoms_1[0], atoms_2[0])
    return distance


def getPairDistance(pdbfile, atom_name):
    """
    Get the distance between every pair of residue's C-alpha/C-beta distance
    :param pdbfile: pdb file
    :param atom_name: CA or CB
    :return: a list of lists, if there is no CB, use CA instead
    """

    p = PDBParser()
    structure = p.get_structure("X", pdbfile)

    # get the sequence length
    seq_len = len(getSeq(pdbfile))

    # initial the matrix
    matrix = [[10000 for _ in range(seq_len)] for _ in range(seq_len)]
    for model in structure:
        # only one chain in each pdb
        for chain in model:
            i, j = 0, 0
            # TODO: could reduce the calculation
            for residue_i in chain:
                for residue_j in chain:
                    # every aa must have one CA but not CB
                    try:
                        atom_i = residue_i[atom_name]
                        atom_j = residue_j[atom_name]
                        matrix[i][j] = atom_i - atom_j
                    except KeyError:

                        # use CA to instead CB
                        atom_i = residue_i["CA"]
                        atom_j = residue_j["CA"]
                        matrix[i][j] = atom_i - atom_j
                        j += 1
                        continue
                    j += 1
                i += 1
                j = 0

        # only use first model
        break

    return matrix


def getTorsion(pdbfile):
    """
    Used to get the dihedral psi and phi of the peptide
    :param pdbfile: pdb file
    :return: a string of dihedral "phi1, psi1; phi2, psi2; ..."
    """
    p = PDBParser()
    structure = p.get_structure("X", pdbfile)

    phi_psi = ""
    for model in structure:
        # get the internal_coordinates to calculate the phi psi
        model.atom_to_internal_coordinates()
        for r in model.get_residues():
            if r.internal_coord:
                # Ca-N -> phi, Ca-C ->psi
                # [phi, psi], the first aa does not has phi,
                # last one does not has psi
                phi_psi = phi_psi + "%s, %s; " % (str(r.internal_coord.get_angle('phi')),
                                                  str(r.internal_coord.get_angle('psi')))
        break

    return phi_psi


def getAtomDistance(atom1, atom2):
    """
    Get two atoms' distance
    :param atom1: atom object in pymol
    :param atom2: atom object in pymol
    :return: float, unit: Å
    """
    pos1 = atom1.coord[0], atom1.coord[1], atom1.coord[2]
    pos2 = atom2.coord[0], atom2.coord[1], atom2.coord[2]
    distance = math.sqrt(sum(list(tuple(map(lambda i, j: pow(i - j, 2), pos1, pos2)))))
    return distance


def getBlosumSim(aa1, aa2, blosum=62):
    """
    Calculate the blosum similarity between two amino acid
    :param aa1: string, one character
    :param aa2:
    :param blosum: int, 62/90
    :return: float
    """
    if aa1 not in aa_3l or aa2 not in aa_3l:
        return 0

    if aa1 == aa2:
        return 1

    idx1 = aa_3l[aa1]
    idx2 = aa_3l[aa2]
    b = blosum62[idx1][idx2]
    if blosum != 62:
        b = blosum90[idx1][idx2]

    # 3 is the highest score for non-identical substitutions, so substract 4 to get into range [-10, -1]
    b = abs(b - 4)

    # map to (0, 1], 1 means exactly same
    b = 1. - (b / 10.0)
    return b


def editDistance(seq1, seq2, blosum=62):
    """
    Modified edit distance to caluclate the cost of changing one seq to another
    :param seq1: string
    :param seq2:
    :param blosum: int, 62/90
    :return: float, costs
    """
    len1 = len(seq1)
    len2 = len(seq2)

    DP = [[0 for _ in range(len1 + 1)]
          for _ in range(2)]

    for i in range(0, len1 + 1):
        DP[0][i] = i
    for i in range(1, len2 + 1):
        for j in range(0, len1 + 1):
            if j == 0:
                # if seq2 is empty, we should initialize the DP
                DP[i % 2][j] = i
            elif seq1[j - 1] == seq2[i - 1]:
                # if the same, then pick the previous information
                DP[i % 2][j] = DP[(i - 1) % 2][j - 1]
            else:
                mini = min(DP[(i - 1) % 2][j],
                           min(DP[i % 2][j - 1],
                               DP[(i - 1) % 2][j - 1]))
                # substitute is the priority
                if mini == DP[(i - 1) % 2][j - 1]:
                    change_cost = 1 - getBlosumSim(seq1[j - 1], seq2[i - 1], blosum=blosum)
                    DP[i % 2][j] = change_cost + mini
                else:
                    DP[i % 2][j] = 1 + mini
    return DP[len2 % 2][len1]


def similar(src_seq, tgt_seq, blosum=62):
    """
    The function to calculate the similarity of two sequences
    :param src_seq: string
    :param tgt_seq:
    :param blosum: int
    :return:
    """
    length = max(len(src_seq), len(tgt_seq))
    sim_value = 1 - editDistance(src_seq, tgt_seq, blosum) / float(length)

    return sim_value


def ssfile2dict(ssfile_txt):
    """
    Convert the second structure files provided by DSSP to dictionary
    :param ssfile_txt: ss file
    :return: a dictionary
            {
                "121p_A": {
                        "1": "-",
                        "2": "E",
                        "3": "E",
                        "4": "E"
                        },
                "5wcv_A": {
                        "1": "-",
                        "2": "E",
                        "3": "E",
                        "4": "E"
                        }
            }
    """
    dict = {}
    below_header = False
    with open(ssfile_txt, "r") as f:
        for line in f:
            if line.strip().startswith("#"):
                below_header = True
                continue
            if below_header and line.split():
                num = int(line.strip().split()[0])
                dict[num] = line[16] if line[16] != " " else "-"
    return dict


def dict2json(output_json, output):
    with open(output, 'w') as fp:
        json.dump(output_json, fp, indent=4)
    print("Saved to JSON file.")


def json2dict(json_file):
    with open(json_file, "r") as f:
        dictionary = json.load(f)
    return dictionary


def getEparts(ssjson):
    """
    Find segments with at least 5 consecutive E(second structure E)'s in the sequence
    :param ssjson: a json file that contains several pdbs' second structure (you can use
                   the functions above to generate from dssp files)
    :return: a dictionary
            {"pdb_chain": [(2,8),(45,50)]}
            which means in pdb_chain, from 2 to 8 residues are E second structure.
    """
    # 把每个pdb中的连续5个及以上的E端拿出来
    ss_map = json2dict(ssjson)
    epart_map = defaultdict(list)

    for pdb_chain, pdb_ss_map in ss_map.items():
        previous = ""
        left = 0
        for position, ss in pdb_ss_map.items():
            position = int(position)
            if ss == "E" and ss != previous:
                left = position
            if ss != "E" and previous == "E":
                # 有至少五个连续的E再放进来
                if position - left > 4:
                    epart_map[pdb_chain].append((left, position - 1))
            previous = ss
    return epart_map


def getPos1Pos2(pairs):
    """
    :param pairs:
    :return:
    """
    pos1_pos2s = []
    # 两个E边界中间的氨基酸个数大于等于1，并且小于等于25
    for i in range(1, len(pairs)):
        if 1 <= pairs[i][0] - pairs[i - 1][1] - 1 <= 25:
            pos1_pos2s.append([pairs[i - 1][1], pairs[i][0]])
    return pos1_pos2s


def screenEpart(pdb_path, epart_json, output_file):
    """
    Used to find backbone which two residues are in different E parts,
    and the number of amino acids between is less than 25 and their CB
    atoms' distance is less than 10.
    :param pdb_path:
    :param epart_json:
    :param output_file:
    :return:
    """
    epart_map = json2dict(epart_json)

    for pdb_chain, pairs in tqdm(epart_map.items()):
        pos1_pos2s = getPos1Pos2(pairs)
        for pos1_pos2 in pos1_pos2s:
            candidate_pdb = []
            cb_distance = get2ReisdueDistance(pdb_path, pdb_chain, pos1_pos2[0], pos1_pos2[1], "CB")
            if cb_distance and cb_distance <= 10:
                candidate_pdb.append(pdb_chain)
                candidate_pdb.append("{}-{}".format(pos1_pos2[0], pos1_pos2[1]))
                candidate_pdb.append(str(cb_distance))
                candidate_pdb.append(getSeq(os.path.join(pdb_path, pdb_chain + ".pdb")))

                with open(output_file, "a") as f:
                    f.write(", ".join(candidate_pdb) + "\n")


def calculateRMSD(origin_file, predict_file):
    """
    Calculate two pdbs' RMSD (C-alpha alignment)
    :param origin_file: pdb file
    :param predict_file: pdb file
    :return: rmsd value
    """
    origin_name = pathlib.Path(origin_file).stem
    predict_name = pathlib.Path(predict_file).stem
    print(origin_file, predict_file)
    cmd.load(origin_file)
    cmd.load(predict_file)

    cmd.select("originlayer", "model " + origin_name)
    cmd.select("predictlayer", "model " + predict_name)

    output = cmd.align("predictlayer////CA", "originlayer////CA")
    cmd.delete("all")
    return output[0]

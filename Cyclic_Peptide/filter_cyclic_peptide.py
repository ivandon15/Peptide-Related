"""
This file is used to find whether a sequence is a cyclic peptide through PDB files.
"""
import re
from Bio import SeqIO
from pymol import cmd
from PDB_Data_Processing.data_process import getAtomDistance


def getCyclic(pdbfile, ss="", check_type="disulfide"):
    """
    Find all the cyclic position of one pdb file
    :param pdbfile: pdb file
    :param ss: string, string of second structures,
                       you can use ssfile2dict function from
                       PDB_Data_Processing/data_process.py to convert by using "".join(list(diction.values()))
                       when check_type is "bond", you need to provide the ss string
    :param check_type: string
           1. "disulfide": only find the disulfide bonds, check if the distance
                           between two CYSs' SG atoms are less
                           than 3A (ref. https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3668-8)
           2. "bond": find all the cyclic bonds (including disulfide ones), check if
                      the distance between any atoms in two residues are less than
                      3A. And exclude the hydrogen bonds using second structures.
    """
    pairs = []

    cmd.load(pdbfile)
    cmd.remove('solvent')

    if check_type == "disulfide":
        # get all the cys
        cmd.select('allCys', 'resn cys')

        # get all the S atom in CYS
        cys_s_atoms = []
        for atom_c in cmd.get_model("allCys").atom:
            if atom_c.name == "SG":
                cys_s_atoms.append(atom_c)

        for i in range(len(cys_s_atoms)):
            for j in range(i + 1, len(cys_s_atoms)):
                # check if each pair cys's SG atoms' distance is less than 3 A
                distance = getAtomDistance(cys_s_atoms[i], cys_s_atoms[j])
                if distance <= 3:
                    # add residue number
                    pairs.append(cys_s_atoms[i].resi + "-" + cys_s_atoms[j].resi)

    elif check_type == "bond":
        if ss == "":
            print("Please enter corresponding second structure")
            return

        seq_len = 0
        for record in SeqIO.parse(pdbfile, "pdb-atom"):
            seq_len = len(record.seq)
        cmd.select('allAA', ' resi 1-' + str(seq_len))

        # get all the cys position
        allposition = set([str(i) for i in range(1, seq_len + 1)])

        for aa in allposition:

            # create multiple layers with each cys
            cmd.select('AA' + aa, 'resi ' + aa)

            # get the atoms that directly bonded with current aa
            cmd.select('AA' + aa + 'nearby', 'neighbor AA' + aa)

            # except for the previous and next
            # consider head tail connection
            if len(cmd.get_model('AA' + aa + 'nearby').atom) >= 3 or (
                    (int(aa) == 1 or int(aa) == int(seq_len)) and len(
                cmd.get_model('AA' + aa + 'nearby').atom) == 2):
                for atom in cmd.get_model('AA' + aa + 'nearby').atom:
                    pair = []
                    try:
                        # not neighbors
                        if int(atom.resi) + 1 != int(aa) and int(atom.resi) - 1 != int(aa):
                            # not ACE?
                            if int(aa) <= int(seq_len) and int(atom.resi) <= int(seq_len):
                                pair.append(aa)
                                pair.append(atom.resi)
                                pair = sorted(pair)
                                # not already found and not hydrogen bonds
                                if "-".join(pair) not in pairs and not hBond(pair, ss):
                                    pairs.append("-".join(sorted(pair)))
                    except ValueError:
                        continue
    cmd.delete('all')
    cmd.reinitialize()


def splitHGI(ss):
    """
    Split a string exclude H, G, I second structures.
    For example, "--GGHHH-EEIIGGG-HHH" ->
    -> ["-","-","GGHHH","-","E","E","IIGGG","-","HHH"]
    :param ss: the string of second structures
    :return: a dictionary help to check if two residues
            are in the same helical parts
    """
    temp = re.split("([^HGI])", ss)
    res = []
    for i in temp:
        if i != "":
            res.append(i)
    current = 0
    HGI_map = {}
    group = 0
    for i in res:
        # go through the list ["-","-","GGHHH","-","E","E","IIGGG","-","HHH"]
        if i.count("H") + i.count("G") + i.count("I") != 0:
            # means current item is helical part
            for j in range(len(i)):
                # current is the index of list.
                # take current=2 ("GGHHH"), go through the GGHHH
                # set index 2,3,4,5,6 to the same number, which means
                # there are the same helical part
                HGI_map[current + j] = group
            group += 1
        current += 1
    return HGI_map


def hBond(pair, ss):
    """
    Check if pair of residues are in the same helical part
    :param pair: a list of position [int,int]
    :param ss: a string second structures of corresponding sequence
    :return: True if in the same helical part (they may be hydrogen bond,
            helical parts contain H, G, I)
    """
    aa1 = pair[0]
    aa2 = pair[1]

    # if not all two residues are in helical parts, then it definitely not a h bond
    if not (ss[aa1] in ["H", "G", "I"] and ss[aa2] in ["H", "G", "I"]):
        return False
    else:
        # both of them are in helical parts, now we need to confirm whether
        # they are in different helical parts
        HGI_map = splitHGI(ss)
        if HGI_map[aa1] != HGI_map[aa2]:
            return False
        else:
            return True

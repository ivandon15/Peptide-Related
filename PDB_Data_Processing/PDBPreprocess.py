"""
This file contains:
1. PDB files download: Extract pdbid from fasta file
2. PDB files modified (extract specific chain)
3. PDB files filter: remove the empty modified files
4. Renumber the atom and residue number

Attention!!
Since the program is using PyMOL, please
download pymol from Schrödinger first!
"""
import os
import pathlib
import urllib.request
from urllib.error import URLError
from Bio import SeqIO
from tqdm import tqdm
from Bio.PDB import PDBParser, PDBIO, Select, is_aa
import warnings

warnings.filterwarnings('ignore')


class ChainResidueSelect(Select):
    """
    Override the Select class in Bio, it has four methods. More details
    plz check 11.1.7 Writing PDB files in
    http://biopython.org/DIST/docs/tutorial/Tutorial.html

    Change the Constructor, add "remain" param which means the chain id
    that we need to remain.
    """

    def __init__(self, remain):
        self.remain = remain

    def accept_chain(self, chain):
        if chain.get_id() == self.remain:
            return True
        else:
            return False

    # uncomment this if you want the peptides only contains
    # standard amino acides
    # if you just want exclude that contains small molecule, delete `standard` 
    # def accept_residue(self, residue,standard=True):
    #     # filter standard aa
    #     return is_aa(residue)


class Pdb:
    """
    A PDB Preprocessing class, including:
    1. download PDB files according to fasta file from RCSB website
    2. extract the specific chain (NOTE: modify the original downloaded pdb files!)
    3. remove the PDB files with only one chain that contains AB residue (like ACYS/BCYS in the pdb files)
    4. renumber the PDB files with only one chain
    """

    def __init__(self, fasta_file="", download_path=""):
        """
        :param fasta_file: original input fasta file
        :param download_path: the folder path you want to save the downloaded pdb files
        """
        self.fasta_file = fasta_file
        self.download_path = download_path
        self.pdblist = []

        files = os.listdir(self.download_path)
        for file in files:
            self.pdblist.append(pathlib.Path(file).stem)

    def extractID(self):
        """
        Extract pdb id and their split info from fasta file and save txt file into export pdb file.
        If you want to output the list to file, you may uncomment the following block.
        >6t9m.BBB|lcl|split_3
        GPAMK
        """
        print("Start extracting PDB ids...")
        for seq_record in SeqIO.parse(self.fasta_file, "fasta"):
            # TODO: you may change your symbol according to your fasta file
            ids = seq_record.id.split("|")
            self.pdblist.append(ids[0])

        # uncomment it if you want save the list to file
        # with open("pdblist.txt", 'w') as f:
        #     for i in pdblist:
        #         f.write("%s\n" % i)
        #     print('Done')

        print("Extracted %d pdb ids with chains." % self.get_list_len())
        return self.pdblist

    def download(self):
        """
        Download the pdb file using pdblist from rcsb website
        :param download_path: output file path
        """
        print("Start downloading...")
        failed = []
        for item in tqdm(self.pdblist):
            # TODO: you may change your symbol. here I was using 1a1p_A
            pdbid, chain = item.split("_")

            try:
                # If you are Linux user, you may use this
                # pdb_url = "https://files.rcsb.org/view/%s.pdb" % pdb_id
                # pdb_file = os.path.join(pdb_path, "%s.pdb" % pdb_id)
                # if not os.path.exists(pdb_file) or os.path.getsize(pdb_file) == 0:
                #     os.system("wget -O %s %s" % (pdb_file, pdb_url))

                urllib.request.urlretrieve('http://files.rcsb.org/download/%s.pdb' % pdbid,
                                           self.download_path + "%s_%s.pdb" % (pdbid, chain))
            except urllib.error.HTTPError as e:
                # throw the HTTP error first
                # if something wrong with downloading, then remove the id in the list
                failed.append(item)
                print("Something wrong with the server.")
                print('Error code: ', e.code)
                print('Error PDB ID: ', pdbid)
                continue
            except urllib.error.URLError as e:
                failed.append(item)
                print("Cannot reach the server.")
                print('Reason: ', e.reason)
                print('Error PDB ID: ', pdbid)
                continue
        for item in failed:
            self.pdblist.remove(item)

        print("Successfully download %d pdb files with chains." % self.get_list_len())

        # now the pdblist only contain valid pdb file
        return self.pdblist

    def isABResidue(self, file):
        """
        Since some peptides may contain something like ACYS/BCYS,
        this function is used to check if the peptide contain A/B residues
        :param file: PDB file
        :return: True/False
        """
        is_AB = False
        with open(file, 'r') as f:
            for line in f:
                if line.startswith("ATOM") and (line[16] == 'A' or line[16] == 'B'):
                    is_AB = True
                    break
        return is_AB

    def modifiedPDBfile(self):
        """
        Modified PDB file in the download_path to fit the pdblist (remain one of the chain and only remain ATOM)
        :param download_path: the folder contains pdb file
        """

        print("Start modifying the pdb files...")
        for item in tqdm(self.pdblist):
            pdbid = item

            # get PDB parser
            parser = PDBParser(PERMISSIVE=1)
            # get the original pdb structure if the file exists
            if os.path.exists(self.download_path + pdbid + '.pdb'):
                if os.path.getsize(self.download_path + pdbid + '.pdb') == 0:
                    os.remove(self.download_path + pdbid + '.pdb')
                else:
                    try:
                        structure = parser.get_structure(pdbid, self.download_path + pdbid + '.pdb')
                        # write new PDB file
                        io = PDBIO()
                        io.set_structure(structure)
                        # TODO: you may change your symbol
                        io.save(self.download_path + pdbid + '.pdb',
                                ChainResidueSelect(pdbid.split("_")[1]))
                    except Exception:
                        os.remove(self.download_path + pdbid + '.pdb')
        print("Before modified, there are %d PDB files." % self.get_list_len())
        print("Checking empty PDB files...")

        # check invalid(empty) pdb files
        removed = []
        self.pdblist = []
        files = os.listdir(self.download_path)
        for file in files:
            self.pdblist.append(pathlib.Path(file).stem)

        for item in self.pdblist:
            pdbid = item
            parser = PDBParser(PERMISSIVE=1)
            structure = parser.get_structure(pdbid, self.download_path + pdbid + '.pdb')
            if len(structure) == 0:
                removed.append(item)
        for item in removed:
            self.pdblist.remove(item)
            os.remove(self.download_path + item.split(", ")[0] + '.pdb')

        print("PDB files have been modified. And remove the empty PDB files, remain %d PDB files" % self.get_list_len())

    def renumber(self):
        """
        Renumber the atom and residues on the download path files
        """

        print("Start renumbering the residue...")
        files = os.listdir(self.download_path)

        for file in tqdm(files):
            out = list()
            with open(self.download_path + file, "r") as f:
                atom_no = 1
                residue_no = 0
                # record the previous residue number
                previous_residue_no = ""
                for line in f:
                    if line.startswith(('ATOM', 'HETATM', 'TER')):
                        # ATOM or HETATM line
                        if len(line.split()) > 5:
                            # get the current residue number, may be int maybe 1C
                            origin_residue_no = line[22:27]
                            if previous_residue_no != origin_residue_no:
                                # count the residue number
                                residue_no += 1
                                previous_residue_no = origin_residue_no

                        # rewrite atom number
                        atom_num = str(atom_no)
                        while len(atom_num) < 5:
                            atom_num = ' ' + atom_num
                        line = '%s%s%s' % (line[:6], atom_num, line[11:])
                        atom_no += 1

                        # rewrite residue number
                        residue_num = str(residue_no)
                        while len(residue_num) < 4:
                            residue_num = ' ' + residue_num
                        new_row = '%s%s' % (line[:22], residue_num)
                        while len(new_row) < 29:
                            new_row += ' '
                        xcoord = line[30:38].strip()
                        while len(xcoord) < 9:
                            xcoord = ' ' + xcoord
                        line = '%s%s%s' % (new_row, xcoord, line[38:])

                        # if multiple chain or model, need to reset the count
                        if line.startswith('TER'):
                            atom_no = 1
                            residue_no = 0
                            previous_residue_no = ""
                    out.append(line)
            # rewrite the pdb file
            with open(self.download_path + file, 'w') as f:
                for line in out:
                    f.write(line)
        print("Renumber process is done.")

    def get_list_len(self):
        return len(self.pdblist)

    def filterABResidue(self):
        """
        Filter out the peptides that contains AB Residue, if you want to keep it,
        you may not use this function
        """
        print("Start filtering out the AB atom...")
        files = os.listdir(self.download_path)
        for file in tqdm(files):
            if self.isABResidue(self.download_path + file):
                os.remove(self.download_path + file)
        self.pdblist = []
        files = os.listdir(self.download_path)
        for file in files:
            self.pdblist.append(pathlib.Path(file).stem)
        print("After AB atom filtering, we got %d files left." % self.get_list_len())

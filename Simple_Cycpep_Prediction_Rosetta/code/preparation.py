import argparse
import os
import pathlib
import sys

from Bio import SeqIO
from tqdm import tqdm

aa_codes = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
    'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'K': 'LYS',
    'I': 'ILE', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
    'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
    'T': 'THR', 'V': 'VAL', 'Y': 'TYR', 'W': 'TRP'}


def getTXTSeq(file):
    """Generate sequence from pdb files"""
    for record in SeqIO.parse(file, "pdb-atom"):
        seq = []
        for aa in list(str(record.seq)):
            seq.append(aa_codes[aa])
        return " ".join(seq)


def pdbs2txts(pdb_path, txt_path):
    """Convert PDB files to txt files"""
    pdb_files = os.listdir(pdb_path)
    for file in tqdm(pdb_files):
        with open(os.path.join(txt_path, pathlib.Path(file).stem + ".txt"), "w") as f:
            f.write(getTXTSeq(os.path.join(pdb_path, file)))


def fasta2txts(fasta_file, txt_path):
    """Convert fasta file to txt files"""
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        pdbid = seq_record.id
        with open(os.path.join(txt_path, pdbid + ".txt"), "w") as f:
            f.write(getTXTSeq(str(seq_record.seq)))


def predictShell(output_shell_file, output_silent_path, input_seq_path, input_pdb_path=""):
    """
    Generate the prediction shell.
    :param output_shell_file:  the output shell file
    :param output_silent_path: the path of output silent files
    :param input_seq_path: the path of input txt sequence files
    :param input_pdb_path: the path of input pdb files(optional), if needed, plz uncomment the last cmd.
    """
    seq_files = os.listdir(input_seq_path)
    for file in seq_files:
        with open(output_shell_file, "a") as f:
            command = "mpirun -np 20 simple_cycpep_predict.mpi.linuxgccrelease " \
                      "-cyclic_peptide:sequence_file " \
                      + os.path.join(input_seq_path, file) + \
                      " -nstruct 100 " \
                      "-cyclic_peptide:genkic_closure_attempts 100 " \
                      "-cyclic_peptide:rama_cutoff 2.0 " \
                      "-cyclic_peptide:MPI_output_fraction 0.01" \
                      "-cyclic_peptide:genkic_min_solution_count 3 " \
                      "-cyclic_peptide:use_rama_filter true " \
                      "-symmetric_gly_tables true " \
                      "-cyclic_peptide:cyclization_type " \
                      "\"terminal_disulfide\" " \
                      "-cyclic_peptide:MPI_auto_2level_distribution " \
                      "-cyclic_peptide:MPI_batchsize_by_level 125 " \
                      "-cyclic_peptide:compute_pnear_to_this_fract 1.0 " \
                      "-cyclic_peptide:MPI_pnear_lambda 1.0 " \
                      "-cyclic_peptide:MPI_pnear_kbt 0.62 " \
                      "-cyclic_peptide:default_rama_sampling_table " \
                      "\"flat_symm_dl_aa_ramatable\" " \
                      "-min_genkic_hbonds 0 " \
                      "-min_final_hbonds 0 " \
                      "-mute all -unmute protocols.cyclic_peptide_predict.SimpleCycpepPredictApplication_MPI_summary " \
                      "-out:file:silent " + os.path.join(output_silent_path, pathlib.Path(file).stem + ".silent")
            # "-in:file:native "+os.path.join(input_pdb_path,pathlib.Path(file).stem + ".pdb")
            f.write(command + ";\n")


def extractPDBShell(output_shell_file, input_silent_path):
    """
    Extract silent files to pdb files, and rename the pdb files' name.
    There is no output path, so plz execute the shell file in the folder that you want to save the output files.
    :param output_shell_file: the output shell file
    :param input_silent_path: the path of input silent files
    """
    input_silent_files = os.listdir(input_silent_path)
    for file in input_silent_files:
        with open(output_shell_file, "a") as f:
            command = "extract_pdbs.mpi.linuxgccrelease " \
                      "-in::file::silent " \
                      + os.path.join(input_silent_path, file) + \
                      ".silent " \
                      "-in:auto_setup_metals\n" \
                      "for file in result_*.pdb; do mv -v \"$file\" \"${file/result/" \
                      + pathlib.Path(file).stem + "}\"; done;"
            f.write(command + "\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Using Rosetta to do the predictions.")
    # if you are using pdb files to convert to sequence files, you may uncomment the following comments
    # and comment the FASTAFILE.

    # parser.add_argument('PDBPATH',
    #                     metavar='path',
    #                     type=str,
    #                     help='the path of the folder of pdb files')
    parser.add_argument('FASTAFILE',
                        metavar='file',
                        type=str,
                        help='the file name of fasta files')

    parser.add_argument('TXTPATH',
                        metavar='path',
                        type=str,
                        help='the output path of the folder of txt files')
    parser.add_argument('SILENTPATH',
                        metavar='path',
                        type=str,
                        help='the output path of the folder of silent files')
    parser.add_argument('PRESHELL',
                        metavar='file',
                        type=str,
                        help='the file name of prediction shell')

    parser.add_argument('XTRCSHELL',
                        metavar='file',
                        type=str,
                        help='the file name of extracting pdb shell')
    args = parser.parse_args()

    if not os.path.isdir(args.PDBPATH) or not os.path.isdir(args.TXTPATH) or not os.path.isdir(args.SILENTPATH):
        print('The path specified does not exist')
        sys.exit()
    else:
        # if you are using pdb files to convert to sequence files, you may uncomment the following comments
        # and comment the uncommented ones.

        # print("Converting PDB files to txt sequence files...")
        print("Converting FASTA file to txt sequence files...")
        fasta2txts(args.FASTAFILE, args.TXTPATH)
        print("DONE in converting FASTA files to txt sequence files!")
        # pdbs2txts(args.PDBPATH, args.TXTPATH)
        # print("DONE in converting PDB files to txt sequence files!")

        print("Generating Prediction shell...")
        predictShell(args.PRESHELL, args.SILENTPATH, args.TXTPATH, args.PDBPATH)
        print("DONE in generating prediction shell!")

        print("Generating Extracting PDB shell...")
        extractPDBShell(args.XTRCSHELL)
        print("DONE in extracting PDB shell!")

# Peptide Related
Contains some peptide data utils tool and some models or tasks related to peptides.

## Cyclic Peptide Related
1. filter_cyclic_peptide.py: You can use this py to find the cyclic postion of one sequence by using its PDB file

## PDB Data Processing
1. PDBPreprocess.py: You can use this py to deal with..
    1. PDB files download: Extract pdbid from fasta file
    2. PDB files modified: extract specific chain
    3. PDB files filter: remove the empty modified files
    4. Renumber the atom and residue number
2. data_process.py:

    1. Get sequence according to PDB file   (with one chain)
    2. Get the distance between 2 specific  residues by specific atom type
    3. Get every pair of residue's C-alpha/ C-beta distance in one sequence
    4. Get each phi and psi of one sequence
    5. Calculate the similarity using Blosum62/ 90 between 2 sequences
    6. Convert the second structure files   provided by DSSP to dictionary
    7. Extract consecutive second structure
    8. Calculate two pdbs' RMSD (C-alpha alignment)

## Simple Cycpep Prediction Rosetta
Usage about simple_cycpep_predict from Rosetta(in chinese), more details please see the readme in this folder.

## Cell Penetrating Peptide
Using 3D Infomax Pre-trained model to predict CPP.
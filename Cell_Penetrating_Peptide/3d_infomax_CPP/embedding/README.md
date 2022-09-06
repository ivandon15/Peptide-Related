# Use 3D Infomax to get Embedding
We use 3D Infomax to get the embedding of each molecule.
## Download Code
Please download 3D Infomax code from https://github.com/HannesStark/3DInfomax
## 2D inference
Currently, 3D Infomax only got a pre-trained model with 2D input. And the input will be converted to mol object through RDKit, so you may use nearly every format of molecular as inputs.

The original code is using **smiles** as input.

(Every format mentioned below has a sample in the *sample* folder.)

Inputs:
1. **smiles**: put all the smiles that you want to generate the embeddings into `dataset/inference_smiles.txt` line by line. And run 
    ```
    python inference.py --config=configs_clean/fingerprint_inference.yml
    ```
    Since the smile format is the original format for the inference of 3D Infomax, there is no other code needed to be modified.

2. **PDB**:
    1. Put all the PDB files' absolute path in the `dataset/inference_pdb.txt` line by line.
    2. Modify the second line in `configs_clean/fingerprint_inference.yml`, change the smile_txt_path to ``dataset/inference_pdb.txt``.
     You may not change the parameter name.
    3. Go to `datasets/inference_dataset.py`, find class `InferenceDataset`, and you should modify the first for loop in order to read every pdb file.
      And you can use  `mol = Chem.MolFromPDBFile(pdb_path)` to get the mol object.

3. **SDF**:
    1. Modify the second line in `configs_clean/fingerprint_inference.yml`, change the smile_txt_path to all sdf file path.
     You may not change the parameter name.
    2. Use `mols = Chem.SDMolSupplier(smiles_txt_path)` instead of reading the content in the path.
    3. Modify the for loop to go through the `mols` above.

4. **Sequence**:
    1. Similar to smiles, and need to use `mol = Chem.MolFromSequence(seq)` to get mol object. Here `seq` is a string object.
    
## Other
`sample/aa_embedding.json` is using 2D net inference to get the embedding of the 20 standard amino acids.
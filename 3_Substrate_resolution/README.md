# SKiD: A Structure-Oriented Kinetics Database of Enzyme-Substrate Interactions
# Authors: Sowmya Ramaswamy Krishnan, Nishtha Pandey, Rajgopal Srinivasan, Arijit Roy*

# Step 3: Substrate resolution
1. map_smiles_v1.py - Program to map SMILES to IUPAC names of substrates collected from the BRENDA database using OPSIN and PubChemPy libraries
2. generate_mol_files.py - Program to generate MMFF94-optimized 3D coordinates (in SDF format with explicit hydrogens) for the substrates from their SMILES using RDKit and OpenBabel libraries
3. 3_resolve_inchi_keys_v1.py - Program to take the SMILES mapping file from the first step and further map the InChIKeys for them using RDKit for redundancy resolution

# Requirements - Preferably a conda environment with all these packages installed
* python>=v3.10.13
* pandas>=2.2.3
* numpy>=1.26.4
* py2opsin
* pubchempy
* openbabel==3.1.1
* rdkit-pypi>=2023.9.6

# Data
Sample input and output files and any data used from external databases (PDB, UniProtKB, EMBL Cofactor database etc.,) are provided under `sample_input` and `sample_output` folders, respectively. The BRENDA (v2023) JSON file and the PDB files of the enzyme-substrate complexes used during the calculations are not included and must be downloaded prior to running the scripts.

# Code usage - Sample commands
```
* python -u 1_map_smiles_v1.py ./sample_input/Unique_ligands_v1.csv ./sample_output/Ligands_all_final_v1.csv <unmapped ligands for manual resolution>
* python -u 2_generate_mol_files.py ./sample_output/Ligands_all_final_v1.csv <path to store ligand SDF files>
* python -u 3_resolve_inchi_keys_v1.py ./sample_output/Ligands_all_final_v1.csv ./sample_output/Ligands_with_InchIKey_v1.csv
```

# Miscellaneous
All inputs to the codes are provided as command-line arguments. Any output folder required to write the results (Ex: Organisms, Proteins etc) must be created by the user before running the scripts. Only paths to stand-alone programs should be modified in the code to ensure reproducibility. Modification of any other parameters is at your own risk.

# License: Creative Commons Attribution Non Commercial No Derivatives 4.0 International












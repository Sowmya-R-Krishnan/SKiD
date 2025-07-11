# SKiD: A Structure-Oriented Kinetics Database of Enzyme-Substrate Interactions
# Authors: Sowmya Ramaswamy Krishnan, Nishtha Pandey, Rajgopal Srinivasan, Arijit Roy*

# Step 4: Omission of unresolved entries
1. omit_apo_proteins_v1.py - Program to remove entries from dataset involving apo proteins whose binding sites could not be automatically/manually resolved
2. omit_small_molecules_v1.py - Program to remove entries from dataset involving small molecular substrates whose SMILES could not be mapped successfully or belonging to one of the following categories: peptides, proteins, nucleic acids, polymers, carbohydrates, peptide-sugar conjugates

# Requirements - Preferably a conda environment with all these packages installed
* python>=v3.10.13
* pandas>=2.2.3
* numpy>=1.26.4

# Data
Sample input and output files and any data used from external databases (PDB, UniProtKB, EMBL Cofactor database etc.,) are provided under `sample_input` and `sample_output` folders, respectively. The BRENDA (v2023) JSON file and the PDB files of the enzyme-substrate complexes used during the calculations are not included and must be downloaded prior to running the scripts.

# Code usage - Sample commands
```
* python -u 1_omit_apo_proteins_v1.py ../1_BRENDA_parsing_and_structuring/sample_output/BRENDA_kcat_data_all_v1.csv ../2_Binding_site_resolution/sample_output/BRENDA_kcat_unique_enzymes_v3.csv ../2_Binding_site_resolution/sample_input/BRENDA_enzyme_binding_sites_all_v1.csv ./sample_output/BRENDA_kcat_after_omissions_final_v1.csv
* python -u 1_omit_apo_proteins_v1.py ../1_BRENDA_parsing_and_structuring/sample_output/BRENDA_Km_data_all_v1.csv ../2_Binding_site_resolution/sample_output/BRENDA_Km_unique_enzymes_v3.csv ../2_Binding_site_resolution/sample_input/BRENDA_enzyme_binding_sites_all_v1.csv ./sample_output/BRENDA_Km_after_omissions_final_v1.csv
* python -u 2_omit_small_molecules_v1.py ./sample_output/BRENDA_kcat_after_omissions_final_v1.csv ./sample_input/BRENDA_omitted_ligands_v1.csv <file to save removed entries> <file to save final kcat dataset>
* python -u 2_omit_small_molecules_v1.py ./sample_output/BRENDA_Km_after_omissions_final_v1.csv ./sample_input/BRENDA_omitted_ligands_v1.csv <file to save removed entries> <file to save final Km dataset>
```

# Miscellaneous
All inputs to the codes are provided as command-line arguments. Any output folder required to write the results (Ex: Organisms, Proteins etc) must be created by the user before running the scripts. Only paths to stand-alone programs should be modified in the code to ensure reproducibility. Modification of any other parameters is at your own risk.

# License: Creative Commons Attribution Non Commercial No Derivatives 4.0 International












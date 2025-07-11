# SKiD: A Structure-Oriented Kinetics Dataset of Enzyme-Substrate Interactions
# Authors: Sowmya Ramaswamy Krishnan, Nishtha Pandey, Rajgopal Srinivasan, Arijit Roy*

# Step 2: Binding site resolution
1. get_unique_enzymes_v1.py - Get a list of all unique enzymes present in the dataset collated from BRENDA in Step 1
2. extract_pdb_chain_mapping_v1.py - Program to extract the number of chains and chain IDs covered by each PDB ID mapped to the dataset
3. categorize_enzyme_binding_sites_v1.py - Program to split the PDB IDs mapped to the dataset into 4 categories (Substrate-only, Cofactor-only, Substrate+Cofactor, Apo structures)
4. correct_ec_chain_mapping_v1.py - Program to find the chain(s) matching the EC number of the enzyme in PDB - to handle cases where multimers have only one/few chains from the enzyme
5. map_holo_to_apo_sites_v1.py - Program to use MSAs from Clustal Omega to map holo to apo sites and cofactor-only sites in the dataset

# Requirements - Preferably a conda environment with all these packages installed
* python>=v3.10.13
* pandas>=2.2.3
* numpy>=1.26.4
* biopython>=1.83

# Other additional stand-alone programs required
* Clustal Omega (pre-compiled binary) - http://www.clustal.org/omega/

# Data
Sample input and output files and any data used from external databases (PDB, UniProtKB, EMBL Cofactor database etc.,) are provided under `sample_input` and `sample_output` folders, respectively. The BRENDA (v2023) JSON file and the PDB files of the enzyme-substrate complexes used during the calculations are not included and must be downloaded prior to running the scripts.

# Code usage - Sample commands
```
* python -u 1_get_unique_enzymes_v1.py ../1_BRENDA_parsing_and_structuring/sample_output/BRENDA_kcat_data_all_v1.csv ./sample_output/BRENDA_kcat_unique_enzymes_v1.csv
* python -u 1_get_unique_enzymes_v1.py ../1_BRENDA_parsing_and_structuring/sample_output/BRENDA_Km_data_all_v1.csv ./sample_output/BRENDA_Km_unique_enzymes_v1.csv
* python -u 2_extract_pdb_chain_mapping_v1.py <path to all PDB IDs in kcat dataset> ./sample_output/kcat_PDB_chain_mapping_v1.csv
* python -u 2_extract_pdb_chain_mapping_v1.py <path to all PDB IDs in Km dataset> ./sample_output/Km_PDB_chain_mapping_v1.csv
* python -u 3_categorize_enzyme_binding_sites_v1.py ./sample_output/BRENDA_kcat_unique_enzymes_v1.csv ./sample_input/PDBe_cofactors.json <path to all PDB IDs in kcat dataset> ./sample_output/BRENDA_kcat_unique_enzymes_v2.csv
* python -u 3_categorize_enzyme_binding_sites_v1.py ./sample_output/BRENDA_Km_unique_enzymes_v1.csv ./sample_input/PDBe_cofactors.json <path to all PDB IDs in Km dataset> ./sample_output/BRENDA_Km_unique_enzymes_v2.csv
* python -u 4_correct_ec_chain_mapping_v1.py ./sample_input/EC_with_entity_ID_v1.csv ./sample_output/kcat_PDB_chain_mapping_v1.csv ../1_BRENDA_parsing_and_structuring/sample_output/BRENDA_kcat_data_all_v1.csv ./sample_output/BRENDA_kcat_unique_enzymes_v3.csv
* python -u 4_correct_ec_chain_mapping_v1.py ./sample_input/EC_with_entity_ID_v1.csv ./sample_output/Km_PDB_chain_mapping_v1.csv ../1_BRENDA_parsing_and_structuring/sample_output/BRENDA_Km_data_all_v1.csv ./sample_output/BRENDA_Km_unique_enzymes_v3.csv
* python -u 5_map_holo_to_apo_sites_v1.py ./sample_input/BRENDA_kcat_apo_only_v1.csv <path to all PDB IDs in kcat dataset> ./sample_input/kcat_MSA/ ./sample_input/BRENDA_enzyme_binding_sites_all_v1.csv <output file>
* python -u 5_map_holo_to_apo_sites_v1.py ./sample_input/BRENDA_Km_apo_only_v1.csv <path to all PDB IDs in Km dataset> ./sample_input/Km_MSA/ ./sample_input/BRENDA_enzyme_binding_sites_all_v1.csv <output file>
```

# Miscellaneous
All inputs to the codes are provided as command-line arguments. Any output folder required to write the results (Ex: Organisms, Proteins etc) must be created by the user before running the scripts. Only paths to stand-alone programs should be modified in the code to ensure reproducibility. Modification of any other parameters is at your own risk.

# License: Creative Commons Attribution Non Commercial No Derivatives 4.0 International












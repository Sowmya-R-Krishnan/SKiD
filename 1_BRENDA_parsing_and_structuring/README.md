# SKiD: A Structure-Oriented Kinetics Database of Enzyme-Substrate Interactions
# Authors: Sowmya Ramaswamy Krishnan, Nishtha Pandey, Rajgopal Srinivasan, Arijit Roy*

# Step 1: BRENDA parsing and structuring
1. parse_brenda_json_v1.py - Script to parse the BRENDA JSON file and separate the organisms, proteins, references, kcat values and Km values for each EC number. BRENDA JSON file must be downloaded from this link: https://www.brenda-enzymes.org/download.php
2. collate_brenda_data_v1.py - Script to combine all details from BRENDA into a single CSV file for each EC number, along with UniProtKB ID -> PDB mapping.
3. resolve_mutations_BRENDA_v1.py - Script to use regular expressions on BRENDA comments section and extract mutation(s) present in the records
4. correct_multi_organism_data_v1.py - Script to split entries with multiple organisms mapped to same record and map the related organism information
5. fill_pdb_for_multi_organism_no_uniprot_v1.py - Script to resolve PDB IDs where UniProtKB ID is not mapped to the entry in BRENDA

# Requirements - Preferably a conda environment with all these packages installed
* python>=v3.10.13
* pandas>=2.2.3
* numpy>=1.26.4
* json5>=0.9.25

# Data
Sample input and output files and any data used from external databases (PDB, UniProtKB, EMBL Cofactor database etc.,) are provided under `sample_input` and `sample_output` folders, respectively. The BRENDA (v2023) JSON file and the PDB files of the enzyme-substrate complexes used during the calculations are not included and must be downloaded prior to running the scripts.

# Code usage - Sample commands
```
* python -u 1_parse_brenda_json_v1.py <path to JSON file> ./sample_output/Organisms/ ./sample_output/Proteins/ ./sample_output/References/ ./sample_output/kcat_values/ ./sample_output/Km_values/
* python -u 2_collate_brenda_data_v1.py ./sample_input/uniprotkb_database_pdb_2024_01_10.tsv ./sample_output/Proteins/ ./sample_output/Organisms/ ./sample_output/References/ ./sample_output/Km_values/ ./sample_output/kcat_values/ ./sample_output/After_collation/
* python -u 3_resolve_mutations_BRENDA_v1.py ./sample_output/After_collation/ ./sample_output/After_collation/ <output file name>
* python -u 4_correct_multi_organism_data_v1.py <output file from previous step> ./sample_output/Organisms/ ./sample_output/Proteins/ ./sample_output/References/ ./sample_input/uniprotkb_database_pdb_2024_01_10.tsv ./sample_output/BRENDA_kcat_data_all_v1.csv ./sample_output/BRENDA_Km_data_all_v1.csv <output file name>
* python -u 5_fill_pdb_for_multi_organism_no_uniprot_v1.py <output file from previous step> ./sample_input/Organism_taxonomy_mapping_v1.csv ./sample_output/Multi_organism_data_no_uniprot_resolved_Km_v1.csv ./sample_output/Multi_organism_data_no_uniprot_resolved_kcat_v1.csv <output file for unresolved entries>
```

# Miscellaneous
All inputs to the codes are provided as command-line arguments. Any output folder required to write the results (Ex: Organisms, Proteins etc) must be created by the user before running the scripts. Only paths to stand-alone programs should be modified in the code to ensure reproducibility. Modification of any other parameters is at your own risk.

# License: Creative Commons Attribution Non Commercial No Derivatives 4.0 International












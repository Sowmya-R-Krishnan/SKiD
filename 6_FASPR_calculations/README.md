# SKiD: A Structure-Oriented Kinetics Database of Enzyme-Substrate Interactions
# Authors: Sowmya Ramaswamy Krishnan, Nishtha Pandey, Rajgopal Srinivasan, Arijit Roy*

# Step 6: FASPR calculations
1. separate_mut_data_complete_v1.py - Program to separate the BRENDA dataset into 3 subsets - WT, Binding site mutations and Non-binding site mutations
2. prepare_FASPR_input_v1.py - Program to extract mutations in a structured format for running FASPR calculations
3. separate_chains_and_seqs_v1.py - Program to create FASPR input sequence files for local structure optimization mode from the input files prepared in previous step and the wild-type PDB structures
4. run_FASPR_v1.py - Program to run FASPR calculations and obtain optimized models of the mutant enzymes
5. extract_FASPR_errors_v1.py - Extract cases where FASPR could not model the mutation due to missing atoms in backbone, chain mapping errors, mutation mapping errors not captured earlier, non-canonical residue at the mutation site etc. For each unmodelled structure, manual verification of the error was performed to resolve the error automatically using FASPR
6. resolve_FASPR_errors_v1.py - Program to resolve manually labelled errors using a combination of OpenMM PDBFixer and FASPR. This program specifically addresses missing backbone atoms, non-canonical to canonical residue conversion and chain ID mismatch issues

# Requirements - Preferably a conda environment with all these packages installed
* python>=v3.10.13
* pandas>=2.2.3
* numpy>=1.26.4
* biopython>=1.83
* openmm>=8.1.1
* pdbfixer>=1.9

# Other additional stand-alone programs required
* FASPR (CLI compiled from source) - https://github.com/tommyhuangthu/FASPR
* PDB2PQR (CLI compiled from source) - https://github.com/Electrostatics/pdb2pqr
* PyMOL (Python bindings) - https://pymol.org/

# Data
Sample input and output files and any data used from external databases (PDB, UniProtKB, EMBL Cofactor database etc.,) are provided under `sample_input` and `sample_output` folders, respectively. The BRENDA (v2023) JSON file and the PDB files of the enzyme-substrate complexes used during the calculations are not included and must be downloaded prior to running the scripts.

# Code usage - Sample commands
```
* python -u 1_separate_mut_data_complete_v1.py ./sample_input/BRENDA_kcat_PDB_level_mut_mapping_v1.csv ../4_Omission_of_unresolved_entries/sample_output/BRENDA_kcat_after_omissions_final_v1.csv ../2_Binding_site_resolution/sample_input/BRENDA_enzyme_binding_sites_all_v1.csv ../5_Mutation_mapping_from_BRENDA_to_PDB/sample_input/BRENDA_kcat_apo_only_v1.csv ../5_Mutation_mapping_from_BRENDA_to_PDB/sample_input/Cofactor_only_structures_v1.csv ../5_Mutation_mapping_from_BRENDA_to_PDB/sample_output/BRENDA_kcat_unique_enzymes_wt_resolved_v1.csv ./sample_output/BRENDA_kcat_wildtype_data_v1.csv ./sample_output/BRENDA_kcat_binding_site_mutations_data_v1.csv ./sample_output/BRENDA_kcat_non_binding_site_mutations_data_v1.csv
* python -u 1_separate_mut_data_complete_v1.py ./sample_input/BRENDA_Km_PDB_level_mut_mapping_v1.csv ../4_Omission_of_unresolved_entries/sample_output/BRENDA_Km_after_omissions_final_v1.csv ../2_Binding_site_resolution/sample_input/BRENDA_enzyme_binding_sites_all_v1.csv ../5_Mutation_mapping_from_BRENDA_to_PDB/sample_input/BRENDA_Km_apo_only_v1.csv ../5_Mutation_mapping_from_BRENDA_to_PDB/sample_input/Cofactor_only_structures_v1.csv ../5_Mutation_mapping_from_BRENDA_to_PDB/sample_output/BRENDA_Km_unique_enzymes_wt_resolved_v1.csv ./sample_output/BRENDA_Km_wildtype_data_v1.csv ./sample_output/BRENDA_Km_binding_site_mutations_data_v1.csv ./sample_output/BRENDA_Km_non_binding_site_mutations_data_v1.csv
* python -u 2_prepare_FASPR_input_v1.py ./sample_output/BRENDA_kcat_binding_site_mutations_data_v1.csv ./sample_input/BRENDA_kcat_PDB_level_mut_mapping_v1.csv ../2_Binding_site_resolution/sample_input/BRENDA_enzyme_binding_sites_all_v1.csv ./sample_output/kcat_FASPR_input_v1.csv 
* python -u 2_prepare_FASPR_input_v1.py ./sample_output/BRENDA_kcat_non_binding_site_mutations_data_v1.csv ./sample_input/BRENDA_kcat_PDB_level_mut_mapping_v1.csv ../2_Binding_site_resolution/sample_input/BRENDA_enzyme_binding_sites_all_v1.csv ./sample_output/kcat_FASPR_input_v1.csv
* python -u 2_prepare_FASPR_input_v1.py ./sample_output/BRENDA_Km_binding_site_mutations_data_v1.csv ./sample_input/BRENDA_Km_PDB_level_mut_mapping_v1.csv ../2_Binding_site_resolution/sample_input/BRENDA_enzyme_binding_sites_all_v1.csv ./sample_output/Km_FASPR_input_v1.csv
* python -u 2_prepare_FASPR_input_v1.py ./sample_output/BRENDA_Km_non_binding_site_mutations_data_v1.csv ./sample_input/BRENDA_Km_PDB_level_mut_mapping_v1.csv ../2_Binding_site_resolution/sample_input/BRENDA_enzyme_binding_sites_all_v1.csv ./sample_output/Km_FASPR_input_v1.csv
* python -u 3_separate_chains_and_seqs_v1.py <wild-type PDB ID> <path to the wild-type PDB structure> <mutation string> ./sample_output/FASPR_input/
* python -u 4_run_FASPR_v1.py ./sample_input/kcat_FASPR_input_v1.csv ./sample_input/FASPR_input/ <path to save mutant structures from FASPR>
* python -u 4_run_FASPR_v1.py ./sample_input/Km_FASPR_input_v1.csv ./sample_input/FASPR_input/ <path to save mutant structures from FASPR>
* python -u 5_extract_FASPR_errors_v1.py ./sample_input/kcat_FASPR_input_v1.csv <path to mutant structures from FASPR> ./sample_input/FASPR_errors_v1.csv
* python -u 5_extract_FASPR_errors_v1.py ./sample_input/Km_FASPR_input_v1.csv <path to mutant structures from FASPR> ./sample_input/FASPR_errors_v2.csv
* python -u 6_resolve_FASPR_errors_v1.py ./sample_input/FASPR_errors_v1.csv <path to wild-type PDB files> ./sample_output/FASPR_input/ <path to save mutant structures from FASPR>
* python -u 6_resolve_FASPR_errors_v1.py ./sample_input/FASPR_errors_v2.csv <path to wild-type PDB files> ./sample_output/FASPR_input/ <path to save mutant structures from FASPR>
```

# Miscellaneous
All inputs to the codes are provided as command-line arguments. Any output folder required to write the results (Ex: Organisms, Proteins etc) must be created by the user before running the scripts. Only paths to stand-alone programs should be modified in the code to ensure reproducibility. Modification of any other parameters is at your own risk.

# License: Creative Commons Attribution Non Commercial No Derivatives 4.0 International












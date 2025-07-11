# SKiD: A Structure-Oriented Kinetics Database of Enzyme-Substrate Interactions
# Authors: Sowmya Ramaswamy Krishnan, Nishtha Pandey, Rajgopal Srinivasan, Arijit Roy*

# Step 5: Mutation mapping from BRENDA to PDB
1. download_SIFTs_mapping_v1.py - Program to automate download of EMBL SIFTS mapping XML files for each PDB ID in the kcat and Km datasets
2. separate_wt_mutants_v1.py - Program to separate wild-type enzyme-substrate pairs from the mutants for mutation mapping
3. map_brenda_mut_to_bs_v1.py - Program to use SIFTS numbering schemes to get an initial mapping of mutations from BRENDA (UniProtKB) to PDB
4. correct_mut_mapping_errors_v1.py - After manual verification and labelling of error types observed in the initial mutation mapping, automated correction of a subset of error types is done through this script

# Requirements - Preferably a conda environment with all these packages installed
* python>=v3.10.13
* pandas>=2.2.3
* numpy>=1.26.4
* biopython>=1.83
* libxml2>=2.12.6

# Data
Sample input and output files and any data used from external databases (PDB, UniProtKB, EMBL Cofactor database etc.,) are provided under `sample_input` and `sample_output` folders, respectively. The BRENDA (v2023) JSON file and the PDB files of the enzyme-substrate complexes used during the calculations are not included and must be downloaded prior to running the scripts.

# Code usage - Sample commands
```
* python -u 1_download_SIFTs_mapping_v1.py <path to all PDB IDs in kcat dataset> <path to all PDB IDs in Km dataset> <path to store the downloaded SIFTS XML files>
* python -u 2_separate_wt_mutants_v1.py ../2_Binding_site_resolution/sample_output/BRENDA_kcat_unique_enzymes_v3.csv ./sample_input/PDB_resolution_year_data_all_v1.csv ./sample_input/BRENDA_kcat_apo_only_v1.csv <path to PDB files in kcat dataset> ./sample_output/BRENDA_kcat_unique_enzymes_wt_resolved_v1.csv
* python -u 2_separate_wt_mutants_v1.py ../2_Binding_site_resolution/sample_output/BRENDA_Km_unique_enzymes_v3.csv ./sample_input/PDB_resolution_year_data_all_v1.csv ./sample_input/BRENDA_Km_apo_only_v1.csv <path to PDB files in Km dataset> ./sample_output/BRENDA_Km_unique_enzymes_wt_resolved_v1.csv
* python -u 3_map_brenda_mut_to_bs_v1.py ../4_Omission_of_unresolved_entries/BRENDA_kcat_after_omissions_final_v1.csv ../2_Binding_site_resolution/BRENDA_enzyme_binding_sites_all_v1.csv <path to SIFTS XML files> ./sample_output/BRENDA_kcat_bs_mutations_v1.csv 
* python -u 3_map_brenda_mut_to_bs_v1.py ../4_Omission_of_unresolved_entries/BRENDA_Km_after_omissions_final_v1.csv ../2_Binding_site_resolution/BRENDA_enzyme_binding_sites_all_v1.csv <path to SIFTS XML files> ./sample_output/BRENDA_Km_bs_mutations_v1.csv
* python -u 4_correct_mut_mapping_errors_v1.py ./sample_output/BRENDA_kcat_bs_mutations_v1.csv ../2_Binding_site_resolution/sample_output/BRENDA_kcat_unique_enzymes_v3.csv ../2_Binding_site_resolution/BRENDA_enzyme_binding_sites_all_v1.csv <path to SIFTS XML files> ./sample_output/BRENDA_kcat_bs_mutations_v1.csv ./sample_input/Cofactor_only_structures_v1.csv ./sample_input/BRENDA_kcat_apo_only_v1.csv <manual error type mapping to BRENDA mutations>
* python -u 4_correct_mut_mapping_errors_v1.py ./sample_output/BRENDA_Km_bs_mutations_v1.csv ../2_Binding_site_resolution/sample_output/BRENDA_Km_unique_enzymes_v3.csv ../2_Binding_site_resolution/BRENDA_enzyme_binding_sites_all_v1.csv <path to SIFTS XML files> ./sample_output/BRENDA_Km_bs_mutations_v1.csv ./sample_input/Cofactor_only_structures_v1.csv ./sample_input/BRENDA_Km_apo_only_v1.csv <manual error type mapping to BRENDA mutations>
```

# Miscellaneous
All inputs to the codes are provided as command-line arguments. Any output folder required to write the results (Ex: Organisms, Proteins etc) must be created by the user before running the scripts. Only paths to stand-alone programs should be modified in the code to ensure reproducibility. Modification of any other parameters is at your own risk.

# License: Creative Commons Attribution Non Commercial No Derivatives 4.0 International












# SKiD: A Structure-Oriented Kinetics Database of Enzyme-Substrate Interactions
# Authors: Sowmya Ramaswamy Krishnan, Nishtha Pandey, Rajgopal Srinivasan, Arijit Roy*

# Step 7: Pre-processing and docking scripts
1. prepare_datasets_final_v1.py - Program to prepare docking dataset based on availability of protein and ligand structure files
2. extract_cofactors_v1.py - Program to separate substrates and cofactors from PDB files 
3. docking_preprocessing_v1.py - Program to retain cofactors in PDB files (WT) and copy cofactors to modelled PDB structures (mutants). The cofactor and substrate mapping from previous step is used as input after thorough manual verification by experts
4. extract_binding_sites_v1.py - Program to extract binding sites for all 4 structure categories using either the bound ligands or through superimposition of homologous structures identified (using PyMOL). The superimposition mode (super vs align) in PyMOL was verified manually for each apo and cofactor-only structure by experts
5. separate_ref_ligands_v1.py - Once the binding sites are mapped to all 4 structure categories, the reference ligand for setting the docking grid (--autobox_ligand in GNINA) is separated as a SDF file through this program
6. protonate_propka_v1.py - Program to use PROPKA tool from PDB2PQR package for protonation of the enzyme structure as per the experimental pH value reported in BRENDA
7. run_docking_calc_v1.py - Program to automate all docking calculations with the pre-processed files using GNINA (multi-pose sampling with CNN-rescoring). The GNINA executable (pre-compiled binary) is assumed to be available in the same directory for this program to work

# Requirements - Preferably a conda environment with all these packages installed
* python>=v3.10.13
* pandas>=2.2.3
* numpy>=1.26.4
* rdkit-pypi>=2023.9.6
* biopython>=1.83
* openmm>=8.1.1
* pdbfixer>=1.9
* json5>=0.9.25

# Other additional stand-alone programs required
* GNINA (pre-compiled binary) - https://github.com/gnina/gnina/releases/tag/v1.3
* PDB2PQR (CLI compiled from source) - https://github.com/Electrostatics/pdb2pqr
* PyMOL (Python bindings) - https://pymol.org/

# Data
Sample input and output files and any data used from external databases (PDB, UniProtKB, EMBL Cofactor database etc.,) are provided under `sample_input` and `sample_output` folders, respectively. The BRENDA (v2023) JSON file and the PDB files of the enzyme-substrate complexes used during the calculations are not included and must be downloaded prior to running the scripts.

# Code usage - Sample commands
```
* python -u 1_prepare_datasets_final_v1.py ../5_Omission_of_unresolved_entries/sample_output/BRENDA_kcat_after_omissions_final_v1.csv ../6_FASPR_calculations/sample_input/BRENDA_kcat_PDB_level_mut_mapping_v1.csv <path to PDB files> <path to ligand SDF files> ./sample_output/BRENDA_kcat_docking_all_v1.csv <path to final kcat dataset>
* python -u 1_prepare_datasets_final_v1.py ../5_Omission_of_unresolved_entries/sample_output/BRENDA_Km_after_omissions_final_v1.csv ../6_FASPR_calculations/sample_input/BRENDA_Km_PDB_level_mut_mapping_v1.csv <path to PDB files> <path to ligand SDF files> ./sample_output/BRENDA_Km_docking_all_v1.csv <path to final Km dataset>
* python -u 2_extract_cofactors_v1.py ./sample_output/BRENDA_kcat_docking_all_v1.csv ../2_Binding_site_resolution/sample_output/BRENDA_kcat_unique_enzymes_v3.csv <path to all PDB IDs in kcat dataset> ../2_Binding_site_resolution/sample_input/PDBe_cofactors.json ./sample_input/Cofactor_substrate_mapping_all_v1.csv
* python -u 2_extract_cofactors_v1.py ./sample_output/BRENDA_Km_docking_all_v1.csv ../2_Binding_site_resolution/sample_output/BRENDA_Km_unique_enzymes_v3.csv <path to all PDB IDs in Km dataset> ../2_Binding_site_resolution/sample_input/PDBe_cofactors.json ./sample_input/Cofactor_substrate_mapping_all_v1.csv
* python -u 3_docking_preprocessing_v1.py ./sample_output/BRENDA_kcat_docking_all_v1.csv ./sample_input/Cofactor_substrate_mapping_all_v1.csv <path to all WT PDB structures> <path to all mutant PDB structures from FASPR> <path to write the processed PDB structures> ./sample_output/BRENDA_kcat_docking_all_v2.csv
* python -u 3_docking_preprocessing_v1.py ./sample_output/BRENDA_Km_docking_all_v1.csv ./sample_input/Cofactor_substrate_mapping_all_v1.csv <path to all WT PDB structures> <path to all mutant PDB structures from FASPR> <path to write the processed PDB structures> ./sample_output/BRENDA_Km_docking_all_v2.csv
* python -u 4_extract_binding_sites_v1.py ./sample_input/Cofactor_substrate_mapping_all_v1.csv ./sample_input/Cofactor_only_list_final_v1.csv ./sample_input/Apo_list_final_v1.csv <path to all processed PDB files> <path to all apo PDB structures> <path to all cofactor-only PDB structures> ./sample_output/Prepped_v1/
* python -u 5_separate_ref_ligands_v1.py ./sample_output/BRENDA_kcat_docking_all_v2.csv ./sample_output/Prepped_v1/ ./sample_output/Prepped_v2/ ./sample_input/Cofactor_substrate_mapping_all_v1.csv ./sample_input/Apo_list_final_v1.csv ./sample_input/Cofactor_only_list_final_v1.csv
* python -u 5_separate_ref_ligands_v1.py ./sample_output/BRENDA_Km_docking_all_v2.csv ./sample_output/Prepped_v1/ ./sample_output/Prepped_v2/ ./sample_input/Cofactor_substrate_mapping_all_v1.csv ./sample_input/Apo_list_final_v1.csv ./sample_input/Cofactor_only_list_final_v1.csv
* python -u 6_protonate_propka_v1.py <path to final kcat dataset> ./sample_output/BRENDA_kcat_docking_all_v2.csv ./sample_output/Prepped_v2/ ./sample_output/Protonated/ <path to write log files from PROPKA>
* python -u 6_protonate_propka_v1.py <path to final Km dataset> ./sample_output/BRENDA_Km_docking_all_v2.csv ./sample_output/Prepped_v2/ ./sample_output/Protonated/ <path to write log files from PROPKA>
* python -u 7_run_docking_calc_v1.py <path to final kcat dataset> ./sample_output/Protonated/ <path to all ligand SDF files generated in Step 3 (Substrate resolution)> <path to write docking output files>
* python -u 7_run_docking_calc_v1.py <path to final Km dataset> ./sample_output/Protonated/ <path to all ligand SDF files generated in Step 3 (Substrate resolution)> <path to write docking output files>
```

# Miscellaneous
All inputs to the codes are provided as command-line arguments. Any output folder required to write the results (Ex: Organisms, Proteins etc) must be created by the user before running the scripts. Only paths to stand-alone programs should be modified in the code to ensure reproducibility. Modification of any other parameters is at your own risk.

# License: Creative Commons Attribution Non Commercial No Derivatives 4.0 International












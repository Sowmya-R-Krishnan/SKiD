# SKiD: A Structure-Oriented Kinetics Database of Enzyme-Substrate Interactions
# Authors: Sowmya Ramaswamy Krishnan, Nishtha Pandey, Rajgopal Srinivasan, Arijit Roy*

# Step 8: Dataset finalization
1. outlier_analysis_v1.py - Program to extract duplicate entries for geometric mean calculations, perform log-transformation of the kcat and Km values, and perform outlier analysis (3*stdev) for dataset finalization. The duplicate entries are carefully analyzed manually and entries with value > 0.5(maximum value) are discarded from geometric mean calculations.

# Requirements - Preferably a conda environment with all these packages installed
* python>=v3.10.13
* pandas>=2.2.3
* numpy>=1.26.4
* pickle

# Data
Sample input and output files and any data used from external databases (PDB, UniProtKB, EMBL Cofactor database etc.,) are provided under `sample_input` and `sample_output` folders, respectively. The BRENDA (v2023) JSON file and the PDB files of the enzyme-substrate complexes used during the calculations are not included and must be downloaded prior to running the scripts.

# Code usage - Sample commands
```
* python -u 1_outlier_analysis_v1.py <path to final kcat dataset> <path to successfully docked entries> <path to failed entries> <path to log-transformed kcat dataset> ./sample_output/kcat_all_data_logscale_final_v1.csv
* python -u 1_outlier_analysis_v1.py <path to final Km dataset> <path to successfully docked entries> <path to failed entries> <path to log-transformed Km dataset> ./sample_output/Km_all_data_logscale_final_v1.csv
```

# Miscellaneous
All inputs to the codes are provided as command-line arguments. Any output folder required to write the results (Ex: Organisms, Proteins etc) must be created by the user before running the scripts. Only paths to stand-alone programs should be modified in the code to ensure reproducibility. Modification of any other parameters is at your own risk.

# License: Creative Commons Attribution Non Commercial No Derivatives 4.0 International












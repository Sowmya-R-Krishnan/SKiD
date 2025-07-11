# SKiD: A Structure-oriented Kinetics Dataset of enzyme-substrate interactions
Authors: Sowmya Ramaswamy Krishnan, Nishtha Pandey, Rajgopal Srinivasan, Arijit Roy*

Automation scripts used for construction of the SKiD dataset arranged into a pipeline. Breaks are present in between the automated processes wherein manual data curation and error resolution was performed by domain experts. The final datasets obtained through this pipeline are available from Zenodo: https://doi.org/10.5281/zenodo.15355031

# Requirements - Preferably a conda environment with all these packages installed
* python>=v3.10.13
* pandas>=2.2.3
* numpy>=1.26.4
* pickle
* matplotlib
* seaborn
* rdkit-pypi>=2023.9.6
* openbabel==3.1.1
* json5>=0.9.25
* biopython>=1.83
* openmm>=8.1.1
* pdbfixer>=1.9
* py2opsin
* pubchempy
* libxml2>=2.12.6

# Other additional stand-alone programs required
* Clustal Omega (pre-compiled binary) - http://www.clustal.org/omega/
* NCBI BLAST+ (CLI from apt) - sudo apt install ncbi-blast+
* GNINA (pre-compiled binary) - https://github.com/gnina/gnina/releases/tag/v1.3
* FASPR (CLI compiled from source) - https://github.com/tommyhuangthu/FASPR
* PDB2PQR (CLI compiled from source) - https://github.com/Electrostatics/pdb2pqr
* PyMOL (Python bindings) - https://pymol.org/

# Usage disclaimer
This is a minimal version of the SKiD dataset curation source code necessary to reproduce the automated steps of the pipeline. Wherever possible, the manual resolution results are also provided. Any changes made to the source code (except paths to stand-alone programs) are done at your own risk. The authors will not be liable to any discrepancies observed in the results due to changes made to the source code.

# Data
Sample input and output files and any data used from external datasets (PDB, UniProtKB, EMBL Cofactor dataset etc.,) are provided under each sub-folder. The BRENDA (v2023) JSON file and the PDB files of the enzyme-substrate complexes used during the calculations are not included and must be downloaded prior to running the scripts.

# Code usage
Detailed instructions on how to use the codes are provided under each sub-folder of the repository. For any queries related to code usage, contact the corresponding author for more information.

# Order of navigation
1. BRENDA parsing and structuring
2. Binding site resolution
3. Substrate resolution
4. Omission of unresolved entries
5. Mutation mapping from BRENDA to PDB
6. FASPR calculations
7. Preprocessing and docking scripts
8. Dataset finalization

# Copyright Notice
SKiD code repository is a TCS proprietary resource and should be used for academic purposes only. The contents of this repository should not be used for any commercial purpose without the consent of ALL the authors involved. By downloading and utilizing the scripts, the user consents that any and all Intellectual Property derived from the SKiD code repository is fully owned by TCS in the associated jurisdictions. SKiD code repository usage without citation will be considered illegal.

# Contact Us
For further queries related to code usage, please write to us: roy.arijit3@tcs.com

# Citation
Please cite this article if you use the codes in this repository for your research: 

# License: Creative Commons Attribution Non Commercial No Derivatives 4.0 International

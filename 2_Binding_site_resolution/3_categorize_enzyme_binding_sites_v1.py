#Program to categorize binding sites in PDB structures

import os
import sys
import pandas as pd
import numpy as np
import Bio
import json
import warnings
from collections import defaultdict
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, PDBList

infile = sys.argv[1]
cofactor_map = sys.argv[2]
pdb_path = sys.argv[3]
outfile = sys.argv[4]

df = pd.read_csv(infile, sep="\t", header=0)
#print(df.columns)  #['EC_number', 'UniProt_ID', 'PDB_IDs', 'Organism_name']

with open(cofactor_map, 'r') as f:
	cofactor_codes = json.load(f)

cofactor_list = list(cofactor_codes.keys())

#Program to check if the ligand (HETATM group) is a cofactor or substrate
def check_cofactor(ec_num, het_id):
	cof_found = 0
	for cof in cofactor_codes.keys():
		substrates = cof[0]['cofactors']
		ec_list = cof[0]['EC']
		if((het_id in substrates) and (ec_num in ec_list)):
			cof_found = cof_found + 1
			break
		else:
			continue
	if(cof_found>0):
		return True  #Input ligand is a cofactor
	else:
		return False  #Input ligand is a substrate/metal ion/inhibitor/activator

for i, row in df.iterrows():
	substrate_pockets = []
	cofactor_pockets = []
	joint_pockets = []
	apo_pockets = []

	pdb_list = str(row["PDB_IDs"]).split(",")  #XXX: Takes into account cases where the PDB ID is a number (XeXX)
	for j, pdb_id in enumerate(pdb_list):
		with warnings.catch_warnings():
			warnings.simplefilter('ignore')
			pdb = PDBParser().get_structure(pdb_id, pdb_path+pdb_id+".pdb")
			io = PDBIO()
			io.set_structure(pdb)

			for model in pdb:
				for chain in model:
					substrate_check = 0
					cofactor_check = 0
					for resi in chain:
						if(resi.id[0].startswith("H_")):
							het_pdb_id = resi.id[0].split("_")[1]
							if(check_cofactor(row["EC_number"], het_pdb_id)):
								cofactor_check = cofactor_check + 1
							else:
								substrate_check = substrate_check + 1

					if(substrate_check>0 and cofactor_check==0):
						substrate_pockets.append(pdb_id+"_"+str(chain.id))
					elif(substrate_check==0 and cofactor_check>0):
						cofactor_pockets.append(pdb_id+"_"+str(chain.id))
					elif(substrate_check>0 and cofactor_check>0):
						joint_pockets.append(pdb_id+"_"+str(chain.id))
					else:
						apo_pockets.append(pdb_id+"_"+str(chain.id))

				break #XXX: For NMR structures, use only the first model for the analysis

		print(substrate_pockets, cofactor_pockets, joint_pockets, apo_pockets)

		break
	break
						
































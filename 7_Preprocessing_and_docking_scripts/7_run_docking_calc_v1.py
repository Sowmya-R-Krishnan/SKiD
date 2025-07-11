#Program to automate docking of Platinum database ligands with their respective receptors

import os
import sys
import csv
import pandas as pd
import numpy as np

infile = sys.argv[1]
protpath = sys.argv[2]
ligpath = sys.argv[3]
outpath = sys.argv[4]

dock_df = pd.read_csv(infile, sep="\t", header=0)
#print(dock_df.columns)  #['EC_number', 'Substrate', 'UniProt_ID', 'Protein_file', 'Organism_name', 'Km_value', 'Mutant', 'Mutation', 'pH', 'Temperature', 'References', 'Site_type', 'Substrate_SMILES', 'Mol_file']

for i, row in dock_df.iterrows():
	prot_name = row["Protein_file"]
	pdb_id = str(prot_name[0:4])
	
	try:
		ph_val = float(row["pH"])
	except:
		ph_val = "nil"

	prot_file = prot_name.replace(".pdb", "")+"_"+str(ph_val)+".pdb"
	#If protonation has failed, need to use the standard pH file
	if(os.path.exists(protpath+prot_file) and os.stat(protpath+prot_file).st_size == 0):
		prot_file = prot_name.replace(".pdb", "")+"_nil.pdb"
	
	lig_file = row["Mol_file"]
	
	ref_lig = ""
	for x in os.listdir(protpath):
		if(x.endswith("_orig.sdf") and x.startswith(pdb_id)):
			ref_lig = x
			break
	
	if(os.path.exists(protpath+prot_file) and os.path.exists(ligpath+lig_file)):
		os.system("./gnina -r "+protpath+prot_file+" -l "+ligpath+lig_file+" --autobox_ligand="+protpath+ref_lig+" -o "+outpath+prot_file.replace(".pdb","")+"_"+lig_file.replace(".sdf","")+"_out.sdf --log="+outpath+prot_file.replace(".pdb","")+"_"+lig_file.replace(".sdf","")+"_log.txt --seed=123456 --num_modes=10")
		print(prot_file, lig_file)




























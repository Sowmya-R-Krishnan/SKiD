#Parse BRENDA datasets and create a single file for every EC number with all data available

import os
import sys
import csv
import pandas as pd
import numpy as np
from collections import defaultdict

uniprot_pdb_path = sys.argv[1]
protpath = sys.argv[2]
orgpath = sys.argv[3]
refpath = sys.argv[4]
km_path = sys.argv[5]
kcat_path = sys.argv[6]
outpath = sys.argv[7]

uniprot_pdb_map = defaultdict(list)
uniprot_df = pd.read_csv(uniprot_pdb_path, sep="\t", header=0)

for i, row in uniprot_df.iterrows():
	uniprot_pdb_map[str(row["Entry"])] = row["PDB"].split(";")[:-1]

for fname in os.listdir(protpath):
	ec_num = fname.replace(".csv", "")
	print(ec_num)

	prot_df = pd.read_csv(protpath+fname, sep="\t", header=0, on_bad_lines='skip')
	ref_df = pd.read_csv(refpath+fname, sep="\t", header=0, on_bad_lines='skip')
	km_df = pd.read_csv(km_path+fname, sep="\t", header=0, on_bad_lines='skip')

	out1 = open(outpath+ec_num+"_Km.csv", "w")
	print('Substrate\tProtein_details\tPDB_IDs\tOrganism_name\tKm_value\tConditions\tReferences', file=out1)

	for i, row in km_df.iterrows():
		if(row['Protein_ID']!="-----"):
			prot_data = prot_df[prot_df['Protein_ID']==int(row['Protein_ID'])]
			if(row['Km_value']!="-----"):
				for j, r2 in prot_data.iterrows():
					uniprot = r2['Accession']
					pdb_data = list(set(uniprot_pdb_map[uniprot]))

					#TODO: Fix UniProt -> PDB mapping deficiency in BRENDA
					if(len(pdb_data)>0):
						pdb_ids = ",".join(list(pdb_data))
						print(row['Substrate']+"\t"+row['Protein_details']+"\t"+pdb_ids+"\t"+row['Organism_name']+"\t"+str(row['Km_value'])+"\t"+row['Conditions']+"\t"+str(row['References']), file=out1)
					else:  #Even UniProt is useless! TODO: Save these later to separate file/folder
						continue

	out1.close()  #XXX: Very important to avoid pandas empty file error

	if(os.path.exists(kcat_path+fname)):
		kcat_df = pd.read_csv(kcat_path+fname, sep="\t", header=0, on_bad_lines='skip')
		out2 = open(outpath+ec_num+"_kcat.csv", "w")
		print('Substrate\tProtein_details\tPDB_IDs\tOrganism_name\tkcat_value\tConditions\tReferences', file=out2)

		for i, row in kcat_df.iterrows():
			if(row["Protein_ID"]!="-----"):
				prot_data = prot_df[prot_df['Protein_ID']==int(row['Protein_ID'])]
				if(row['kcat_value']!="-----"):
					for j, r2 in prot_data.iterrows():
						uniprot = r2['Accession']
						pdb_data = list(set(uniprot_pdb_map[uniprot]))

						#TODO: Fix UniProt -> PDB mapping deficiency in BRENDA
						if(len(pdb_data)>0):
							pdb_ids = ",".join(list(pdb_data))
							print(row['Substrate']+"\t"+row['Protein_details']+"\t"+pdb_ids+"\t"+row['Organism_name']+"\t"+str(row['kcat_value'])+"\t"+row['Conditions']+"\t"+str(row['References']), file=out2)
						else:  #Even UniProt is useless! TODO: Save these later to separate file/folder
							continue

		out2.close()  #XXX: Very important to avoid pandas empty file error

		kcat_check = pd.read_csv(outpath+ec_num+"_kcat.csv", sep="\t", header=0)
		if(len(kcat_check.index)==0):  #Empty file
			os.remove(outpath+ec_num+"_kcat.csv")

	km_check = pd.read_csv(outpath+ec_num+"_Km.csv", sep="\t", header=0)
	if(len(km_check.index)==0):  #Empty file
		os.remove(outpath+ec_num+"_Km.csv")



































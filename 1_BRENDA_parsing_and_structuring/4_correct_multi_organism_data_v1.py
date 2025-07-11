#Program to separate the data for multi-organism cases to map the PDB IDs correctly

import os
import sys
import csv
import re
import pandas as pd
import numpy as np
from collections import defaultdict

data = sys.argv[1]
organism_path = sys.argv[2]
protein_path = sys.argv[3]
ref_path = sys.argv[4]
uniprot_pdb_path = sys.argv[5]
kcat_outfile = sys.argv[6]
km_outfile = sys.argv[7]
no_uniprot_outfile = sys.argv[8]

input_df = pd.read_csv(data, sep="\t", header=0)  #EC_number, Substrate, Protein_details, PDB_IDs, Organism_name, Kinetic_value, Conditions, References, Kinetic_type

uniprot_pdb_map = defaultdict(list)
uniprot_df = pd.read_csv(uniprot_pdb_path, sep="\t", header=0)
for i, row in uniprot_df.iterrows():
	uniprot_pdb_map[str(row["Entry"])] = row["PDB"].split(";")[:-1]

out1 = open(kcat_outfile, "w")
out2 = open(km_outfile, "w")
out3 = open(no_uniprot_outfile, "w")

print("EC_number\tSubstrate\tProtein_details\tPDB_IDs\tOrganism_name\tkcat_value\tConditions\tReferences", file=out1)
print("EC_number\tSubstrate\tProtein_details\tPDB_IDs\tOrganism_name\tKm_value\tConditions\tReferences", file=out2)
print("EC_number\tSubstrate\tOrganism_name\tKinetic_value\tConditions\tReferences\tKinetic_type", file=out3)

for i, row in input_df.iterrows():
	ec_num = row["EC_number"]
	orgn_df = pd.read_csv(organism_path+ec_num+".csv", sep="\t", header=0, on_bad_lines='skip')  
	prot_df = pd.read_csv(protein_path+ec_num+".csv", sep="\t", header=0, on_bad_lines='skip')
	ref_df = pd.read_csv(ref_path+ec_num+".csv", sep="\t", header=0, on_bad_lines='skip')

	substrate = row["Substrate"]
	kinetic_value = row["Kinetic_value"]
	refs = str(row["References"]).split(",")
	kinetic_type = row["Kinetic_type"]
	organisms = row["Conditions"].split("; ")

	ref_dict = {}
	for ref_id in refs:
		try:
			ref_sub = ref_df[ref_df["ID"]==int(ref_id)]
			ref_dict[int(ref_id)] = ref_sub["Title"].item()+","+str(ref_sub["Year"].item())
		except:
			ref_dict[int(ref_id)] = ""
			continue

	for i, orgn in enumerate(organisms):
		print(ec_num, orgn)
		try:
			orgn_matches = re.search("#{1}(\d+)#{1}", orgn)  
			orgn_id = orgn_matches.group(1)
			orgn_id_list = [orgn_id]
		except:
			try: 
				orgn_matches = re.search("#{1}([0-9\,]{1,})#{1}", orgn)
				orgn_id_list = orgn_matches.group(1).split(",")
				orgn_id_list = [int(x) for x in orgn_id_list]
			except:
				print(ec_num+"\t"+substrate+"\t"+orgn_name+"\t"+str(kinetic_value)+"\t"+orgn+"\t"+str(ref_dict)+"\t"+kinetic_type, file=out3)
				continue

		for orgn_id in orgn_id_list:
			orgn_name_df = orgn_df[orgn_df["Organism_ID"]==int(orgn_id)]
			orgn_name = orgn_name_df["Organism_name"].item()

			uniprot_id = "-----"
			try:
				protein_df = prot_df[prot_df["Protein_ID"]==int(orgn_id)]
				uniprot_id = protein_df["Accession"].item()
				pdb_ids = list(set(uniprot_pdb_map[uniprot_id]))
				pdb_str = ",".join(list(pdb_ids))
				if(kinetic_type=="kcat"):
					print(ec_num+"\t"+substrate+"\t"+str(uniprot_id)+"\t"+str(pdb_str)+"\t"+orgn_name+"\t"+str(kinetic_value)+"\t"+orgn+"\t"+str(ref_dict), file=out1)
				else:
					print(ec_num+"\t"+substrate+"\t"+str(uniprot_id)+"\t"+str(pdb_str)+"\t"+orgn_name+"\t"+str(kinetic_value)+"\t"+orgn+"\t"+str(ref_dict), file=out2)
			except:
				print(ec_num+"\t"+substrate+"\t"+orgn_name+"\t"+str(kinetic_value)+"\t"+orgn+"\t"+str(ref_dict)+"\t"+kinetic_type, file=out3)
				continue






























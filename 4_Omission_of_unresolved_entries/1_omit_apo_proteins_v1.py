#Program to compute the final datasets (Km and kcat) after omitting the apo proteins without PDB mapping

import sys
import csv
import pandas as pd
import numpy as np
import ast
from ast import literal_eval

data = sys.argv[1]
enzyme_list = sys.argv[2]
binding_sites = sys.argv[3]
outfile1 = sys.argv[4]
outfile2 = sys.argv[5]

main_df = pd.read_csv(data, sep="\t", header=0)
enzyme_df = pd.read_csv(enzyme_list, sep="\t", header=0)
bsites = pd.read_csv(binding_sites, sep="\t", header=None)

out1 = open(outfile1, "w")
out2 = open(outfile2, "w")

print("EC_number\tSubstrate\tProtein_details\tPDB_IDs\tOrganism_name\tKm_value\tMutant\tMutation\tpH\tTemperature\tReferences", file=out1)
print("EC_number\tSubstrate\tProtein_details\tPDB_IDs\tOrganism_name\tKm_value\tMutant\tMutation\tpH\tTemperature\tReferences\tSite_type", file=out2)

for i, row in main_df.iterrows():
	site_type = ""

	try:
		uniprot_dict = literal_eval(row["Protein_details"])
		uniprot_list = uniprot_dict["accessions"]
	except:
		uniprot_list = [row["Protein_details"]]

	ec_num = row["EC_number"]
	orgn_name = row["Organism_name"]

	pdb_df = enzyme_df[enzyme_df["EC_number"]==ec_num]

	pdb_found = 0
	for j, r2 in pdb_df.iterrows():
		if(r2["UniProt_ID"] in uniprot_list):
			pdb_found = pdb_found + 1
			sub_only = r2["Substrate_only"].split(",")
			cof_only = r2["Cofactor_only"].split(",")
			sub_and_cof = r2["Substrate+Cofactor"].split(",")
			apo = r2["Apo_structures"].split(",")

			check_valid_site = 0
			#Substrate+Cofactor > Substrate > Cofactor > Apo
			if(sub_and_cof[0]!="-"):
				sites_available = list(bsites[0])
				for site in sub_and_cof:
					if(site in sites_available):
						site_type = "Substrate+Cofactor"
						check_valid_site = 1
						break

			elif(sub_only[0]!="-"):
				sites_available = list(bsites[0])
				for site in sub_only:
					if(site in sites_available):
						site_type = "Substrate_only"
						check_valid_site = 1
						break

			elif(cof_only[0]!="-"):
				sites_available = list(bsites[0])
				for site in cof_only:
					if(site in sites_available):
						site_type = "Cofactor_only"
						check_valid_site = 1
						break

			elif(apo[0]!="-"):
				sites_available = list(bsites[0])
				for site in apo:
					if(site in sites_available):
						site_type = "Apo site"
						check_valid_site = 1
						break

			else:  #No structures mapped
				contents = row.to_numpy().flatten().tolist()
				contents = [str(x) for x in contents]
				print("\t".join(contents), file=out1)

			if(check_valid_site==1):  #Structures and binding sites mapped successfully
				contents = row.to_numpy().flatten().tolist()
				contents = [str(x) for x in contents]
				print("\t".join(contents)+"\t"+site_type, file=out2)
			else:  #Structures mapped, but binding sites not mapped successfully
				contents = row.to_numpy().flatten().tolist()
				contents = [str(x) for x in contents]
				print("\t".join(contents), file=out1)

			break  #Reduces time to check UniProt match when first entry itself matches to the list. Ensure no duplicates in file
		else:
			continue

	if(pdb_found==0):  #To handle cases where multiple EC numbers have the same UniProt ID mapped in BRENDA
		pdb_list = row["PDB_IDs"].split(",")
		site_found = 0
		for site in pdb_list:
			sites_available = list(bsites[0])
			for s in sites_available:
				if(s.startswith(site)):
					site_found = site_found + 1
					break

		if(site_found>0):
			print("\t".join(contents)+"\tAmbiguous", file=out2)
		else:
			print(ec_num+"\t"+str(uniprot_list)+"\t"+orgn_name)
			print("\t".join(contents), file=out1)




























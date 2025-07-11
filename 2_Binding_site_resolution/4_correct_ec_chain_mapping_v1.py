#Program to extract PDB-chain pairs which are actually mapped to the EC number of interest

import sys
import csv
import pandas as pd
import numpy as np

ec_map = sys.argv[1]
chain_map = sys.argv[2]
infile = sys.argv[3]
outfile = sys.argv[4]

ec_df = pd.read_csv(ec_map, sep=",", header=0)
chain_df = pd.read_csv(chain_map, sep="\t", header=0)
data_df = pd.read_csv(infile, sep="\t", header=0)

out = open(outfile, "w")
print("EC_number\tUniProt_ID\tOrganism_name\tSubstrate_only\tCofactor_only\tSubstrate+Cofactor\tApo_structures", file=out)

for i, row in data_df.iterrows():
	new_sub_sites = []
	new_cof_sites = []
	new_sub_cof_sites = []
	new_apo_sites = []

	sub_sites = row["Substrate_only"].split(",")
	cof_sites = row["Cofactor_only"].split(",")
	sub_cof_sites = row["Substrate+Cofactor"].split(",")
	apo_sites = row["Apo_structures"].split(",")

	if(sub_sites[0]=="-"):
		new_sub_sites.append("-")
	else:
		for site in sub_sites:
			pdb_id = site.split("_")[0]
			chain_id = site.split("_")[1]

			chains_all = chain_df[chain_df["PDB_ID"]==str.upper(pdb_id)]["Chains"].item()
			chains_list = chains_all.split(",")
			#print(chains_list)

			ec_sub_df = ec_df[ec_df["PDB ID"]==pdb_id]

			ec_found = 0
			entity_id = 0
			entity_list = []

			entity_init = 1
			for j, entry in ec_sub_df.iterrows():
				for e in range(int(entry["Total Number of polymer Entity Instances (Chains) per Entity"])):
					entity_list.append(entity_init)

				entity_init = entity_init + 1

				if(entry["EC Number"]=="-"):
					continue
				else:
					ec_found = ec_found + 1
					entity_id = int(entry["Entity ID"])
			
			if(ec_found==0):
				new_sub_sites.append(site)
			else:
				for eid, cid in zip(entity_list, chains_list):
					if(entity_id==eid and cid==chain_id):
						new_sub_sites.append(site)
						break

	if(cof_sites[0]=="-"):
		new_cof_sites.append("-")
	else:
		for site in cof_sites:
			pdb_id = site.split("_")[0]
			chain_id = site.split("_")[1]

			chains_all = chain_df[chain_df["PDB_ID"]==str.upper(pdb_id)]["Chains"].item()
			chains_list = chains_all.split(",")
			
			ec_sub_df = ec_df[ec_df["PDB ID"]==pdb_id]

			ec_found = 0
			entity_id = 0
			entity_list = []

			entity_init = 1
			for j, entry in ec_sub_df.iterrows():
				for e in range(int(entry["Total Number of polymer Entity Instances (Chains) per Entity"])):
					entity_list.append(entity_init)

				entity_init = entity_init + 1

				if(entry["EC Number"]=="-"):
					continue
				else:
					ec_found = ec_found + 1
					entity_id = int(entry["Entity ID"])
			
			if(ec_found==0):
				new_cof_sites.append(site)
			else:
				for eid, cid in zip(entity_list, chains_list):
					if(entity_id==eid and cid==chain_id):
						new_cof_sites.append(site)
						break

	if(sub_cof_sites[0]=="-"):
		new_sub_cof_sites.append("-")
	else:
		for site in sub_cof_sites:
			pdb_id = site.split("_")[0]
			chain_id = site.split("_")[1]

			chains_all = chain_df[chain_df["PDB_ID"]==str.upper(pdb_id)]["Chains"].item()
			chains_list = chains_all.split(",")
			
			ec_sub_df = ec_df[ec_df["PDB ID"]==pdb_id]

			ec_found = 0
			entity_id = 0
			entity_list = []

			entity_init = 1
			for j, entry in ec_sub_df.iterrows():
				for e in range(int(entry["Total Number of polymer Entity Instances (Chains) per Entity"])):
					entity_list.append(entity_init)

				entity_init = entity_init + 1

				if(entry["EC Number"]=="-"):
					continue
				else:
					ec_found = ec_found + 1
					entity_id = int(entry["Entity ID"])
			
			if(ec_found==0):
				new_sub_cof_sites.append(site)
			else:
				for eid, cid in zip(entity_list, chains_list):
					if(entity_id==eid and cid==chain_id):
						new_sub_cof_sites.append(site)
						break	
					

	if(apo_sites[0]=="-"):
		new_apo_sites.append("-")
	else:
		for site in apo_sites:
			pdb_id = site.split("_")[0]
			chain_id = site.split("_")[1]

			chains_all = chain_df[chain_df["PDB_ID"]==str.upper(pdb_id)]["Chains"].item()
			chains_list = chains_all.split(",")
			
			ec_sub_df = ec_df[ec_df["PDB ID"]==pdb_id]

			ec_found = 0
			entity_id = 0
			entity_list = []

			entity_init = 1
			for j, entry in ec_sub_df.iterrows():
				for e in range(int(entry["Total Number of polymer Entity Instances (Chains) per Entity"])):
					entity_list.append(entity_init)

				entity_init = entity_init + 1

				if(entry["EC Number"]=="-"):
					continue
				else:
					ec_found = ec_found + 1
					entity_id = int(entry["Entity ID"])
			
			if(ec_found==0):
				new_apo_sites.append(site)
			else:
				for eid, cid in zip(entity_list, chains_list):
					if(entity_id==eid and cid==chain_id):
						new_apo_sites.append(site)
						break

	#print(new_sub_sites, new_cof_sites, new_sub_cof_sites, new_apo_sites)
	print(str(row["EC_number"])+"\t"+str(row["UniProt_ID"])+"\t"+str(row["Organism_name"])+"\t"+",".join(new_sub_sites)+"\t"+",".join(new_cof_sites)+"\t"+",".join(new_sub_cof_sites)+"\t"+",".join(new_apo_sites), file=out)

	#break





























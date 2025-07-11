#Redundancy removal, log-conversion and outlier analysis for the Km and kcat datasets

import os
import sys
import csv
import pandas as pd
import numpy as np
import statistics
from collections import defaultdict
from statistics import geometric_mean
import ast
from ast import literal_eval
import pickle
import math

infile = sys.argv[1]
pass_list = sys.argv[2]
fail_list = sys.argv[3]
outfile1 = sys.argv[4]
outfile2 = sys.argv[5]
outfile3 = sys.argv[6]

#----------------------------------------------------------------------------------------------------------------------------------------
#Stage 1: Redundancy removal - Geometric mean approach
df = pd.read_csv(infile, sep="\t", header=0)
#print(df.columns)  #['EC_number', 'Substrate', 'UniProt_ID', 'Protein_file', 'Organism_name', 'Km_value', 'Mutant', 'Mutation', 'pH', 'Temperature', 'References', 'Site_type', 'Substrate_SMILES', 'Mol_file']

pass_df = pd.read_csv(pass_list, sep="\t", header=0)
fail_df = pd.read_csv(fail_list, sep="\t", header=0)


duplicate_entries = defaultdict(list)
for i, r1 in df.iterrows():
	for j, r2 in df.iterrows():
		if(i!=j):
			#XXX: Protein file name takes care of matching mutation, pH, temperature between the duplicated rows
			if((r1["Protein_file"]==r2["Protein_file"]) and (r1["Mol_file"]==r2["Mol_file"]) and r1["pH"]==r2["pH"] and r1["Temperature"]==r2["Temperature"] and (r1["Km_value"]!=r2["Km_value"])):
				if(r1["Km_value"] not in duplicate_entries[(r1["EC_number"], r1["Organism_name"], r1["Protein_file"], r1["Mol_file"], r1["pH"], r1["Temperature"])]):
					duplicate_entries[(r1["EC_number"], r1["Organism_name"], r1["Protein_file"], r1["Mol_file"], r1["pH"], r1["Temperature"])].append((r1["Km_value"], r1["References"]))
				if(r2["Km_value"] not in duplicate_entries[(r1["EC_number"], r1["Organism_name"], r1["Protein_file"], r1["Mol_file"], r1["pH"], r1["Temperature"])]):
					duplicate_entries[(r1["EC_number"], r1["Organism_name"], r1["Protein_file"], r1["Mol_file"], r1["pH"], r1["Temperature"])].append((r2["Km_value"], r2["References"]))
	print(i)

with open("./Outlier_analysis_results/Km_WT_duplicate_entries.pkl", "wb") as f:
	pickle.dump(duplicate_entries, file=f)

pass_entries = defaultdict(list)
fail_entries = defaultdict(list)

for i, row in pass_df.iterrows():
	pass_entries[(row["Protein_file"], row["Mol_file"], row["pH"], row["Temperature"])] = literal_eval(row["Km_value_with_ref"])
	
for i, row in fail_df.iterrows():
	#print(row)
	if(row["Value_to_be_retained_for_model"]!="-----"):
		fail_entries[(row["Protein_file"], row["Mol_file"], row["pH"], row["Temperature"])] = literal_eval(row["Value_to_be_retained_for_model"])

nr_df = pd.DataFrame(columns=list(df.columns))
done_dups = []

for i, row in df.iterrows():
	if((row["Protein_file"], row["Mol_file"], row["pH"], row["Temperature"]) in list(pass_entries.keys())):
		if((row["Protein_file"], row["Mol_file"], row["pH"], row["Temperature"], row["References"]) not in done_dups):
			Km_list = pass_entries[(row["Protein_file"], row["Mol_file"], row["pH"], row["Temperature"])]
			km_val = [x[0] for x in Km_list]
			new_Km = geometric_mean(km_val)
			
			row["Km_value"] = new_Km
			row_df = pd.DataFrame([row])
			nr_df = pd.concat([nr_df, row_df], ignore_index=True)
			done_dups.append((row["Protein_file"], row["Mol_file"], row["pH"], row["Temperature"], row["References"]))
			#print(Km_list, new_Km)
		else:
			continue
			
	elif((row["Protein_file"], row["Mol_file"], row["pH"], row["Temperature"]) in list(fail_entries.keys())):
		if((row["Protein_file"], row["Mol_file"], row["pH"], row["Temperature"], row["References"]) not in done_dups):
			Km_list = fail_entries[(row["Protein_file"], row["Mol_file"], row["pH"], row["Temperature"])]
			try:
				km_val = [x[0] for x in Km_list]
			except:
				km_val = [Km_list[0]]
			#print(Km_list, km_val)
				
			new_Km = geometric_mean(km_val)
			
			row["Km_value"] = new_Km
			row_df = pd.DataFrame([row])
			nr_df = pd.concat([nr_df, row_df], ignore_index=True)
			done_dups.append((row["Protein_file"], row["Mol_file"], row["pH"], row["Temperature"], row["References"]))
		else:
			continue
			
	elif((row["EC_number"], row["Organism_name"], row["Protein_file"], row["Mol_file"], row["pH"], row["Temperature"]) in list(duplicate_entries.keys())):
		#print(row["EC_number"], row["Organism_name"], row["Protein_file"], row["Mol_file"], row["pH"], row["Temperature"])
		continue  #Omit the entry	
		
	else:
		row_df = pd.DataFrame([row])
		nr_df = pd.concat([nr_df, row_df], ignore_index=True)
		
print(nr_df)
nr_df.to_csv(outfile1, sep="\t", header=True, index=False)

#----------------------------------------------------------------------------------------------------------------------------------------
#Stage 2: Log-transformation of Km and kcat values

df = pd.read_csv(infile, sep="\t", header=0)
#print(df.columns)  #['EC_number', 'Substrate', 'UniProt_ID', 'Protein_file', 'Organism_name', 'Km_value', 'Mutant', 'Mutation', 'pH', 'Temperature', 'References', 'Site_type', 'Substrate_SMILES', 'Mol_file']

col_list = list(df.columns)
col_list.append("pKm_value")
log_df = pd.DataFrame(columns=col_list)
for i, row in df.iterrows():
	row["pKm_value"] = math.log10(float(row["Km_value"]))
	row_df = pd.DataFrame([row])
	log_df = pd.concat([log_df, row_df], ignore_index=True)
	
print(log_df)
log_df.to_csv(outfile2, sep="\t", header=True, index=False)

#----------------------------------------------------------------------------------------------------------------------------------------
#Stage 3: Outlier removal based on log-scale Km and kcat values - TODO: Confirm ranges based on literature and discussion
df = pd.read_csv(infile, sep="\t", header=0)
#print(df.columns)  #['EC_number', 'Substrate', 'UniProt_ID', 'Protein_file', 'Organism_name', 'Km_value', 'Mutant', 'Mutation', 'pH', 'Temperature', 'References', 'Site_type', 'Substrate_SMILES', 'Mol_file', pkcat_value]

data = list(df["pKm_value"])
mean = np.mean(data)
stdev = np.std(data)
stdev_tripled = 3*stdev

#pass_dataset = 0
pass_df = pd.DataFrame(columns=list(df.columns))
fail_df = pd.DataFrame(columns=list(df.columns))
del_count = 0
for i, val in enumerate(data):
	if(val>=mean-stdev_tripled and val<=mean+stdev_tripled):
		#pass_dataset = pass_dataset + 1
		row_df = pd.DataFrame([df.iloc[i]])
		pass_df = pd.concat([pass_df, row_df], ignore_index=True)
	else:
		del_count = del_count + 1
		row_df = pd.DataFrame([df.iloc[i]])
		fail_df = pd.concat([fail_df, row_df], ignore_index=True)
		
print(del_count)
fail_df.to_csv(outfile3, sep="\t", header=True, index=False)



























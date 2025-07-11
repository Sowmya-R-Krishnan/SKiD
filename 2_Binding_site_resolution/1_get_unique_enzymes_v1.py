#Program to get the list of unique enzymes present in the curated BRENDA database

import sys
import csv
import pandas as pd
import numpy as np
import ast
from ast import literal_eval

infile = sys.argv[1]
outfile = sys.argv[2]

df = pd.read_csv(infile, sep="\t", header=0)

out = open(outfile, "w")
print("EC_number\tUniProt_ID\tPDB_IDs\tOrganism_name", file=out)

uniprot_list = []
orgn_names = []

for i, row in df.iterrows():
	orgn_name = row["Organism_name"]
	if(orgn_name not in orgn_names):
		orgn_names.append(orgn_name)

	try:
		prot_det = literal_eval(row["Protein_details"])
		uniprot_id = prot_det["accessions"][0]
	except:
		uniprot_id = row["Protein_details"]

	if(uniprot_id not in uniprot_list):
		uniprot_list.append(uniprot_id)
		print(row["EC_number"]+"\t"+uniprot_id+"\t"+row["PDB_IDs"]+"\t"+row["Organism_name"], file=out)

print("No. of unique proteins: "+str(len(uniprot_list)))
print("No. of unique organisms: "+str(len(orgn_names)))

































#Program to compute the final datasets (Km and kcat) after omitting the problematic small molecules (without SMILES, polymers, proteins, peptides, nucleic acids)

import sys
import csv
import pandas as pd
import numpy as np
import ast
from ast import literal_eval

data = sys.argv[1]
mol_list = sys.argv[2]
outfile1 = sys.argv[3]
outfile2 = sys.argv[4]

main_df = pd.read_csv(data, sep="\t", header=0)
#print(main_df.columns)  #['EC_number', 'Substrate', 'Protein_details', 'PDB_IDs', 'Organism_name', 'Km_value', 'Mutant', 'Mutation', 'pH', 'Temperature', 'References', 'Site_type']

mol_df = pd.read_csv(mol_list, sep="\t", header=0)
#print(mol_df.columns)  #['SMILES', 'Substrate', 'Mol_ID']

out1 = open(outfile1, "w")
out2 = open(outfile2, "w")

print("EC_number\tSubstrate\tProtein_details\tPDB_IDs\tOrganism_name\tKm_value\tMutant\tMutation\tpH\tTemperature\tReferences\tSite_type", file=out1)
print("EC_number\tSubstrate\tProtein_details\tPDB_IDs\tOrganism_name\tKm_value\tMutant\tMutation\tpH\tTemperature\tReferences\tSite_type\tSubstrate_SMILES\tMol_ID", file=out2)

for i, row in main_df.iterrows():
	mol_names = list(mol_df["Substrate"])

	if(row["Substrate"] in mol_names):
		substrate_df = mol_df[mol_df["Substrate"]==row["Substrate"]]
		smiles = substrate_df.iloc[0]["SMILES"]
		mol_id = substrate_df.iloc[0]["Mol_ID"]

		contents = row.to_numpy().flatten().tolist()
		contents = [str(x) for x in contents]
		print("\t".join(contents)+"\t"+smiles+"\t"+str(mol_id), file=out2)

	else:
		contents = row.to_numpy().flatten().tolist()
		contents = [str(x) for x in contents]
		print("\t".join(contents), file=out1)














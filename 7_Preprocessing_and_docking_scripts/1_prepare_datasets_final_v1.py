#Program to prepare final Km and kcat datasets for docking calculations

import os
import sys
import csv
import re
import pandas as pd
import numpy as np
from collections import defaultdict

input_data = sys.argv[1]
mut_map = sys.argv[2]
checkpath1 = sys.argv[3]
checkpath2 = sys.argv[4]
outfile1 = sys.argv[5]
outfile2 = sys.argv[6]

data_df = pd.read_csv(input_data, sep="\t", header=0)
#print(data_df.columns)  #['EC_number', 'Substrate', 'UniProt_ID', 'WT_PDB', 'Mutant_PDB', 'Mutant_PDB_only', 'Organism_name', 'kcat_value', 'Mutant', 'Mutation', 'pH', 'Temperature', 'References', 'Site_type', 'Substrate_SMILES', 'Mol_ID']

three_to_one = {"ALA":"A", "GLY":"G", "VAL":"V", "LEU":"L", "ILE":"I", "PHE":"F", "TRP":"W", "TYR":"Y", "HIS":"H", "ARG":"R", "LYS":"K", "ASN":"N", "GLN":"Q", "ASP":"D", "GLU":"E", "SER":"S", "CYS":"C", "THR":"T", "MET":"M", "PRO":"P", "OCS":"C", "CSO":"C", "MSE":"M", "MHS":"H", "CME":"C", "CSD":"C", "TSY":"C", "CSS":"C", "SCY":"C", "PTR":"Y", "DDZ":"A", "JJJ":"C", "MHO":"M", "TPO":"T", "SEP":"S", "KPI":"K", "SLZ":"K", "KPF":"K", "KPY":"K", "74P":"K", "VPV":"K", "KGC":"K", "KYQ":"K", "LYF":"K", "CIR":"R", "SNC":"C", "CSP":"C", "LLP":"K", "KCX":"K"}

pdb_mutmap = defaultdict(dict)
mut_df = pd.read_csv(mut_map, sep="\t", header=0)
for i, row in mut_df.iterrows():
	pdb = row["PDB_ID"]
	if(pdb!="-" and row["Filename"]!="-"):
		pdb_mutmap[pdb][row["Mutation"]] = row["Filename"]

protnames = os.listdir(checkpath1)
lignames = os.listdir(checkpath2)

out1 = open(outfile1, "w")
out2 = open(outfile2, "w")

print("Protein_file\tLigand_file\tData_type", file=out1)
print("EC_number\tSubstrate\tUniProt_ID\tProtein_file\tOrganism_name\tKm_value\tMutant\tMutation\tpH\tTemperature\tReferences\tSite_type\tSubstrate_SMILES\tMol_file", file=out2)

missing_pdb = []
for i, row in data_df.iterrows():
	ligcheck = 0
	protcheck = 0
	
	ligname = ""
	protname = ""
	prot_type = ""
	
	if(row["Mutant"]=="yes"):
		prot_type = "MUT"
		if("mol_"+str(row["Mol_ID"])+".sdf" in lignames):
			ligcheck = 1
			ligname = "mol_"+str(row["Mol_ID"])+".sdf"
		else:
			print(str(row["Mol_ID"])+" not found.")
			
		if(row["WT_PDB"]!="-"):
			pdb_id = row["WT_PDB"]
		else:
			pdb_id = row["Mutant_PDB"]
		
		try:
			protname = pdb_mutmap[pdb_id][row["Mutation"]]
			protcheck = 1
		except:
			if(pdb_id!="-" and pdb_id not in missing_pdb):
				missing_pdb.append(pdb_id)
				
			continue
	else:
		prot_type = "WT"
		prots = [x for x in protnames if x.startswith(row["WT_PDB"]) or x.startswith(row["Mutant_PDB"])]
		if("mol_"+str(row["Mol_ID"])+".sdf" in lignames):
			ligcheck = 1
			ligname = "mol_"+str(row["Mol_ID"])+".sdf"
		else:
			print(str(row["Mol_ID"])+" not found.")
			
		if(len(prots)>0):
			if(row["WT_PDB"]!="-"):
				protname = row["WT_PDB"]+".pdb"
			else:
				protname = row["Mutant_PDB"]+".pdb"
				
			protcheck = 1
		else:
			#print(row["WT_PDB"], row["Mutant_PDB"])
			if(row["WT_PDB"]!="-" and row["WT_PDB"] not in missing_pdb):
				missing_pdb.append(row["WT_PDB"])
			if(row["Mutant_PDB"]!="-" and row["Mutant_PDB"] not in missing_pdb):
				missing_pdb.append(row["Mutant_PDB"])
			
	if(ligcheck==1 and protcheck==1):
		print(str(protname)+"\t"+ligname+"\t"+prot_type, file=out1)
		print(row["EC_number"]+"\t"+row["Substrate"]+"\t"+row["UniProt_ID"]+"\t"+str(protname)+"\t"+row["Organism_name"]+"\t"+str(row["Km_value"])+"\t"+row["Mutant"]+"\t"+row["Mutation"]+"\t"+str(row["pH"])+"\t"+str(row["Temperature"])+"\t"+str(row["References"])+"\t"+row["Site_type"]+"\t"+row["Substrate_SMILES"]+"\t"+ligname, file=out2)
			
print(",".join(missing_pdb))





























#Program to extract the wild-type and mutant PDB structures ordered by their resolution for each enzyme

import sys
import csv
import pandas as pd
import numpy as np
import Bio
import warnings
import re

from Bio import PDB
from Bio.PDB import PDBIO, PDBParser

infile = sys.argv[1]
resolution_data = sys.argv[2]
apo_mapping = sys.argv[3]
inpath = sys.argv[4]
outfile = sys.argv[5]

df = pd.read_csv(infile, sep="\t", header=0)
res_df = pd.read_csv(resolution_data, sep=",", header=0)
apo_df = pd.read_csv(apo_mapping, sep="\t", header=0)

apo_holo_map = {}
for i, row in apo_df.iterrows():
	apo_holo_map[row["Apo_structures"]] = row["Holo_mapping"]

out = open(outfile, "w")
print("EC_number\tUniProt_ID\tOrganism_name\tWT_structures\tMutant_structures\tStructure_type", file=out)

for i, row in df.iterrows():
	struct_type = ""
	#print(row['EC_number'], row['UniProt_ID'])
	pdb_list = []
	pdbs = []
	if(row['Substrate+Cofactor']!="-"):
		pdbs = row['Substrate+Cofactor'].split(",")
		pdb_list = [x.split("_")[0] for x in pdbs]
		struct_type = "Substrate+Cofactor"
	elif(row['Substrate_only']!="-"):
		pdbs = row['Substrate_only'].split(",")
		pdb_list = [x.split("_")[0] for x in pdbs]
		struct_type = "Substrate only"
	elif(row['Cofactor_only']!="-"):
		pdbs = row['Cofactor_only'].split(",")
		pdb_list = [x.split("_")[0] for x in pdbs]
		struct_type = "Cofactor only"
	elif(row['Apo_structures']!="-"):
		try:
			pdbs = row['Apo_structures']
			holo_list = [apo_holo_map[pdbs]]
			pdb_list = [x.split("_")[0] for x in holo_list]	
			struct_type = "Apo only"	
		except:
			continue
	else:
		continue

	pdb_list = list(set(pdb_list))
	res_sub_df = res_df[res_df["PDB ID"].isin(pdb_list)]

	best_struct = []
	if(len(pdb_list)==1):
		best_struct.append(pdb_list[0])
	else:
		res_sub_df.sort_values(by="Resolution (Å)", axis=0, inplace=True)
		best_struct = list(res_sub_df["PDB ID"])

	#Extracting wild-type vs mutant information from PDB files
	wt_struct = []
	mut_struct = []
	for j, struct in enumerate(best_struct):
		mut = 0
		contents = []
		
		with open(inpath+struct+".pdb") as pfile:
			for line in pfile.readlines():
				line = line.strip()
				if(line.endswith(" MUTATION: YES") or line.endswith(" ENGINEERED: YES")):  #Can be any COMPND header
					mut = 1
				contents.append(line)
				
		if(mut==0):
			wt_struct.append(str(struct))
		else:
			mutations = []
			for line in contents:
				if(line.startswith("SEQADV")):
					res = line.split(" ")
					res = [x for x in res if x!=""]
					if((res[-1].endswith("MUTATION") or res[-1].endswith("VARIANT")) and res[7] in ["GLY","ALA","VAL","LEU","ILE","PHE","TRP","TYR","CYS","SER","THR","MET","ARG","HIS","LYS","ASP","ASN","GLU","GLN","PRO"]):
						mut_resi = (res[2]+"_"+res[4]+"_"+res[3], res[6], res[7]+"_"+res[8])
						mutations.append(mut_resi)
					elif((res[-1].endswith("MUTATION") or res[-1].endswith("VARIANT")) and res[7] not in ["GLY","ALA","VAL","LEU","ILE","PHE","TRP","TYR","CYS","SER","THR","MET","ARG","HIS","LYS","ASP","ASN","GLU","GLN","PRO"]):
						mut_resi = (res[2]+"_"+res[4]+"_"+res[3], res[6])		
					elif(res[-1].endswith("INSERTION")):
						mut_resi = ("INS_"+res[2]+"_"+res[4]+"_"+res[3], res[6])
						mutations.append(mut_resi)
					elif(res[-1].endswith("DELETION")):
						mut_resi = (res[-4], "DEL_"+res[-3]+"_"+res[-2])
						mutations.append(mut_resi)
					else:  #EXPRESSION TAG, ...
						continue

			mut_struct.append(str(struct))
			
	if(len(mut_struct)==0):
		mut_struct.append("-")
	if(len(wt_struct)==0):
		wt_struct.append("-")
		
	print(row['EC_number']+"\t"+row['UniProt_ID']+"\t"+row['Organism_name']+"\t"+",".join(wt_struct)+"\t"+",".join(mut_struct)+"\t"+struct_type, file=out)




























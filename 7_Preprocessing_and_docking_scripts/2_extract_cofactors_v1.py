#Program to extract ligand names (substrate, cofactor) etc from PDB files

import os
import sys
import csv
import pandas as pd
import json
import numpy as np
from Bio import PDB
from Bio.PDB import PDBIO, PDBParser, Select
from collections import defaultdict
import warnings

data = sys.argv[1]
ec_mapping = sys.argv[2]
inpath = sys.argv[3]
cofactors = sys.argv[4]
outfile = sys.argv[5]

df = pd.read_csv(data, sep="\t", header=0)
#print(df.columns)  #['Protein_file', 'Ligand_file', 'Data_type']

ec_df = pd.read_csv(ec_mapping, sep="\t", header=0)
#print(ec_df.columns)  #['EC_number', 'UniProt_ID', 'Organism_name', 'WT_structures', 'Mutant_structures', 'Structure_type']

with open(cofactors) as f:
	cofactor_data = json.load(f)

cofactor_dict = defaultdict(list)
for key, value in cofactor_data.items():
	ec_list = value[0]["EC"]
	cofs = value[0]["cofactors"]
	
	for ec in ec_list:
		cofactor_dict[ec].extend(cofs)
		
out = open(outfile, "w")
print("EC_number\tPDB_ID\tStructure_type\tSubstrates\tCofactors", file=out)
		
ions = ['NCO', 'NA', 'RHD', 'FE', 'SR', 'PO4', 'ZN', 'MG', 'AG', 'K', 'AU3', 'NI', 'SE4', 'ACT', 'IRI', 'IR', 'CL', 'LU', 'TB', 'CA', 'CS', 'TL', 'SO4', 'CU', 'CAC', 'AU', 'BR', 'SM', 'CD', 'MN', 'NH4', 'IR3', 'CO', 'HG', '3CO', 'OS', 'CON', 'BA', 'FE2', 'PB', 'SIN', 'MES', 'GOL', 'MPD', 'BME', '1PE', 'EPE', 'OXY', 'O', 'NO3', 'FES', 'SF4', 'MOS', 'TRS', 'YT3', 'LI', 'DIO', 'MRD', 'CO3', 'F', 'LI']

nc_resi = ["OCS", "CSO", "MSE", "MHS", "CME", "CSD", "TSY", "CSS", "SCY", "PTR", "DDZ", "JJJ", "MHO", "TPO", "SEP", "KPI", "SLZ", "KPF", "KPY", "74P", "VPV", "KGC", "KYQ", "LYF", "CIR", "SNC", "CSP", "LLP", "KCX", "MLY", "MLZ", "OHI", "HTR", "MTY", "OAS", "CSX", "OMT", "TRQ", "CXM", "SCH", "ALY", "TYX", "CAS", "IAS", "SNN", "NEP", "D4B", "2DA", "YCM", "CSR", "HS8", "A0A", "ALS", "TYS", "C3Y", "ASL", "ACE", "DTY", "OMY", "OMZ", "PCA"]

done = []

for i, row in df.iterrows():
	try:
		pdb_id = str.upper(row["Protein_file"][0:4])
		if(pdb_id not in done):
			struct_type = ""
			
			for j, r in ec_df.iterrows():
				if(pdb_id in r["WT_structures"].split(",") or pdb_id in r["Mutant_structures"].split(",")):
					ec_num = r["EC_number"]
					break
			
			#print(pdb, ec_num)
			cof_list = []
			try:
				cof_list = cofactor_dict[ec_num]
				if(len(cof_list)==0):
					cof_list.append("NA")
			except:
				cof_list.append("NA")
				
			#print(pdb_id, cof_list)
			
			cof_dict = defaultdict(list)
			sub_dict = defaultdict(list)
			with warnings.catch_warnings():
				warnings.simplefilter('ignore')
				pdb = PDBParser().get_structure(pdb_id, inpath+pdb_id+".pdb") #PDBParser module to get structure
				io = PDBIO()
				io.set_structure(pdb)
				
				for model in pdb:
					for chain in model:
						for resi in chain:
							if(resi.id[0].startswith("H_")):
								het_id = resi.resname
								if(cof_list[0]!="NA" and het_id not in ions and het_id not in nc_resi):
									if(het_id in cof_list):
										cof_dict[pdb_id].append((chain.id, het_id))
									else:
										sub_dict[pdb_id].append((chain.id, het_id))
								elif(het_id not in ions and het_id not in nc_resi):
									sub_dict[pdb_id].append((chain.id, het_id))
								else:
									continue
				
				if(len(cof_dict.keys())>0 and len(sub_dict.keys())>0):
					struct_type = "Substrate+Cofactor"
				elif(len(cof_dict.keys())==0 and len(sub_dict.keys())>0):
					struct_type = "Substrate only"
				elif(len(cof_dict.keys())>0 and len(sub_dict.keys())==0):
					struct_type = "Cofactor only"
				else:
					struct_type = "Apo structure"		
				#print(cof_dict, sub_dict)
								
				print(ec_num+"\t"+str(pdb_id)+"\t"+struct_type+"\t"+str(sub_dict[pdb_id])+"\t"+str(cof_dict[pdb_id]), file=out)
				done.append(pdb_id)
	except:
		continue
	
	#break































#Extracting all information for mutation mapping to PDB structures in BRENDA dataset

import sys
import csv
import pandas as pd
import numpy as np
from collections import defaultdict
import xml.etree.ElementTree as ET
import ast
from ast import literal_eval

data = sys.argv[1]
pdb_final_data = sys.argv[2]
bs_data = sys.argv[3]
sifts_path = sys.argv[4]
outfile = sys.argv[5]
cofactor_resolve_df = sys.argv[6]
apo_resolve_df = sys.argv[7]
error_corrections = sys.argv[8]

df = pd.read_csv(data, sep="\t", header=0)
#print(df.columns)  #['EC_number', 'Substrate', 'Protein_details', 'PDB_IDs', 'Organism_name', 'kcat_value', 'Mutant', 'Mutation', 'pH', 'Temperature', 'References', 'Site_type', 'Substrate_SMILES', 'Mol_ID']

pdb_mapping_df = pd.read_csv(pdb_final_data, sep="\t", header=0)
#print(pdb_mapping_df.columns)  #[EC_number, UniProt_ID	Organism_name	Substrate_only	Cofactor_only	Substrate+Cofactor	Apo_structures]

cofactor_df = pd.read_csv(cofactor_resolve_df, sep="\t", header=0)  #XXX: Consider only proximal cofactor sites and not distal sites
#print(cofactor_df.columns)  #[EC_number, UniProt_ID, Organism_name, Cofactor_only, Substrate_site, Reference_PDB, Comments]
cofactor_sub_df = cofactor_df[cofactor_df["Substrate_site"]=="Proximal"]

apo_df = pd.read_csv(apo_resolve_df, sep="\t", header=0)
#print(apo_df.columns)  #[EC_number, UniProt_ID, Organism_name, Apo_structures, Holo_mapping]

error_res_df = pd.read_csv(error_corrections, sep="\t", header=None)
#for i, row in error_res_df.iterrows():
	#print(row[0])
	#print(int(row[1]))

bs_df = pd.read_csv(bs_data, sep="\t", header=None)
three_to_one = {"GLY":"G", "ALA":"A", "VAL":"V", "LEU":"L", "ILE":"I", "PRO":"P", "MET":"M", "SER":"S", "TYR":"Y", "TRP":"W", "PHE":"F", "THR":"T", "HIS":"H", "LYS":"K", "ARG":"R", "ASP":"D", "ASN":"N", "GLU":"E", "GLN":"Q", "CYS":"C"}
one_to_three = dict((c,i) for i, c in three_to_one.items())

#--------------------------------------------------------------SIFTS MAPPING---------------------------------------------------------------
#Get SIFTS mapping between author-provided residue numbering and UniProt residue numbering
def get_sifts_mapping_v1(pdb_id, chain_target, uniprot_id, sifts_path):
	mod_res = {"SME":"MET", "TOX":"TRP"}
	
	chain_mapping = defaultdict(dict)
	pdb_id = str.lower(pdb_id)
	try:
		tree = ET.parse(sifts_path+pdb_id+'.xml')
		root = tree.getroot()
		parent_map = {c: p for p in tree.iter() for c in p}
		#print(parent_map)
		
		#Entity tags in the file correspond to unique chains in the PDB file
		for child in root:
			if(child.tag.endswith("entity")):
				chain_id = child.attrib["entityId"]	
				chain_mapping[chain_id] = {}				
				resi_map = {}
				
				for sub_elements in child.iter():
					if(sub_elements.tag.endswith("crossRefDb") and sub_elements.attrib["dbSource"]=="UniProt" and sub_elements.attrib["dbAccessionId"]==uniprot_id):
						uniprot_id = sub_elements.attrib["dbAccessionId"]
						uniprot_resinum = sub_elements.attrib["dbResNum"]
						uniprot_resiname = sub_elements.attrib["dbResName"]
						uniprot_resid = uniprot_resiname+"_"+str(uniprot_resinum)
						
						#mapped_resi = parent_map[sub_elements]
						parent = parent_map[sub_elements]
						mapped_resi = [ele for ele in parent.iter() if(ele.tag.endswith("crossRefDb") and ele.attrib["dbSource"]=="PDB")]
						pdbe_resinum = mapped_resi[0].attrib["dbResNum"]
						#null residue numbers mean the residue is unmodelled (missing) in PDB structure
						if(pdbe_resinum!="null"):
							pdbe_resiname = mapped_resi[0].attrib["dbResName"]
							if(pdbe_resiname in list(mod_res.keys())):  #XXX: Fixes error code 13
								pdbe_resiname = mod_res[pdbe_resiname]
								
							pdbe_resid = pdbe_resiname+"_"+str(pdbe_resinum)
						else:
							continue
						
						resi_map[uniprot_resid] = pdbe_resid
				
				chain_mapping[chain_id] = resi_map
				
		return chain_mapping[chain_target]  #XXX: A collections defaultdict object with SIFTS PDB -> UniProt residue numbering
	except:
		return {}
		
#Get SIFTS mapping between sequential residue numbering and UniProt residue numbering
def get_sifts_mapping_v2(sifts_chain_map):
	if(sifts_chain_map!={}):
		sequential_resimap = {}
		index = 1
		for key, value in sifts_chain_map.items():
			resiname = value.split("_")[0]
			sequential_resimap[resiname+"_"+str(index)] = key #UniProt residue mapped to PDB residue with sequential numbering
			index = index + 1
			
		return sequential_resimap
	else:
		return {}

#Find out if residue is mutated in PDB already - Ex: 5NIW: T30V - 'VAL_30': 'T_52'
def check_mutated_resi(wt_resi, mut_resi, mut_num, resi_map_inverted):
	pdbe_resi = "invalid"
	for pdb_res, unip_res in resi_map_inverted.items():
		res_num1 = pdb_res.split("_")[1]  #30
		res_name1 = pdb_res.split("_")[0]  #VAL
		
		res_num2 = unip_res.split("_")[1]  #52
		res_name2 = unip_res.split("_")[0]  #T
		
		try:
			if(res_num1==str(mut_num) and three_to_one[res_name1]!=wt_resi and res_name2==wt_resi):
				pdbe_resi = pdb_res
				break
		except:
			continue
	
	return pdbe_resi
	
#------------------------------------------------------------------------------------------------------------------------------------------
out = open(outfile, "w")
print("EC_number\tUniProt_ID\tOrganism_name\tStructure_type\tPDB_ID\tChain_ID\tBRENDA_mutation\tPDB_mutation\tMutation_type", file=out)

for i, row in pdb_mapping_df.iterrows():
	ec_num = row["EC_number"]
	uniprot_id = row["UniProt_ID"]
	orgn_name = row["Organism_name"]
	
	struct_type = ""
	if(row["Substrate+Cofactor"]!="-"):
		pdb_list = row["Substrate+Cofactor"].split(",")
		struct_type = "Substrate+Cofactor"
	elif(row["Substrate_only"]!="-"):
		pdb_list = row["Substrate_only"].split(",")
		struct_type = "Substrate only"
	elif(row["Cofactor_only"]!="-"):
		struct_type = "Cofactor only"
		
		pdb_list = []
		for j, r in cofactor_sub_df.iterrows():
			if(r["Cofactor_only"]==row["Cofactor_only"]):
				#pdb_id = r["Reference_PDB"].split("_")[0]
				#pdb_list.append(pdb_id)
				#pdb_list.append(r["Reference_PDB"])
				pdb_list = row["Cofactor_only"].split(",")
				break
		
	elif(row["Apo_structures"]!="-"):
		struct_type = "Apo only"
		
		pdb_list = []
		for j, r in apo_df.iterrows():
			if(r["Apo_structures"]==row["Apo_structures"]):
				#pdb_id = r["Holo_mapping"].split("_")[0]
				#pdb_list.append(pdb_id)
				pdb_list.append(r["Holo_mapping"])
	else:  #No structures mapped to PDB
		print(ec_num+" "+uniprot_id+" "+orgn_name+" cannot be resolved due to structures unavailable in PDB.")
		continue
		
	if(len(pdb_list)>0):
		#print(ec_num, uniprot_id, orgn_name, pdb_list, struct_type)
		mut_df = df[(df["EC_number"]==ec_num) & (df["Organism_name"]==orgn_name) & (df["Mutant"]=="yes")]
		mutations_unique = []
		prot_details = ""
		for k, mut in mut_df.iterrows():
			try:
				prot_id = literal_eval(mut["Protein_details"])['accessions']
			except:
				prot_id = mut["Protein_details"].split(",")
			
			if(uniprot_id in prot_id):
				if(not mut["Protein_details"].startswith("{")):
					prot_details_1 = mut["Protein_details"]
					prot_details_2 = "{'accessions': ['"+mut["Protein_details"]+"'], 'source': 'uniprot'}"
				else:
					prot_details_1 = mut["Protein_details"]
					prot_details_2 = uniprot_id

				mutations = mut["Mutation"].split("/")
				for m in mutations:
					if(m not in mutations_unique):
						try:
							wt_resi = m[0]
							mut_resi = m[-1]
							mutation_num = m[1:-1]
							if(wt_resi in list(one_to_three.keys()) and mut_resi in list(one_to_three.keys()) and mutation_num.isdigit()):  #Fixes error code 6
								mutations_unique.append(m)
						except:
							continue
		#print(mutations_unique)				
		
		bs_mutations = defaultdict(list)
		non_bs_mutations = defaultdict(list)
		
		for j, pdb_id in enumerate(pdb_list):
			pdb = pdb_id.split("_")[0]
			chain = pdb_id.split("_")[1]
			resi_map = get_sifts_mapping_v1(pdb, chain, uniprot_id, sifts_path)
			resi_map_inverted = dict((c,i) for i,c in resi_map.items())
			#print(pdb, chain, resi_map.keys())
			seq_resi_map = get_sifts_mapping_v2(resi_map)
			#print(seq_resi_map)
			
			binding_site = bs_df[bs_df[0]==str(pdb)+"_"+str(chain)].iloc[0][1]
			
			if(len(resi_map.keys())>0 and len(mutations_unique)>0 and len(binding_site)>0):
				#print(pdb, chain, resi_map)
				for mut in mutations_unique:					
					check_str1 = ec_num+" "+orgn_name+" "+prot_details_1+" "+pdb+" "+chain+" "+mut
					check_str2 = ec_num+" "+orgn_name+" "+prot_details_2+" "+pdb+" "+chain+" "+mut
						
					for ind, err in error_res_df.iterrows():
						val1 = err[0]
						val2 = int(err[1])
						val3 = err[2]
						if(val1.startswith(check_str1) and val2==2): #Wrongly curated errors corrected manually
							mut = val3  #Fixes error code 2
						elif(val1.startswith(check_str2) and val2==2):
							mut = val3  #Fixes error code 2
					
					wt_resi = mut[0]
					mut_resi = mut[-1]
					mut_num = int(mut[1:-1])
					
					if(wt_resi+"_"+str(mut_num) in resi_map):
						pdbe_resi = resi_map[wt_resi+"_"+str(mut_num)]
						check_res = pdbe_resi.replace("_", "")
					elif(one_to_three[wt_resi]+"_"+str(mut_num) in resi_map_inverted):  #Fixes error code 8
						pdbe_resi = resi_map_inverted[one_to_three[wt_resi]+"_"+str(mut_num)]
						resparts = pdbe_resi.split("_")
						check_res = one_to_three[resparts[0]]+str(resparts[1])
					elif(one_to_three[wt_resi]+"_"+str(mut_num) in seq_resi_map):  #Fixes error code 12
						uniprot_resi = seq_resi_map[one_to_three[wt_resi]+"_"+str(mut_num)]
						pdbe_resi = resi_map[uniprot_resi]
						check_res = pdbe_resi.replace("_", "")
					else:
						pdbe_resi = check_mutated_resi(wt_resi, mut_resi, mut_num, resi_map_inverted)
						#print(pdb, chain, mut, pdbe_resi)
						#pdbe_resi = "invalid"  #Can be invalid due to multiple reasons: Look at the 15 error codes
						if(pdbe_resi!="invalid"):
							#print(pdb, chain, mut, pdbe_resi)
							resparts = pdbe_resi.split("_")
							check_res = resparts[0]+str(resparts[1])
						
					#print(pdb, chain, mut, pdbe_resi)
					if(pdbe_resi!="invalid"):
						if(check_res in binding_site):
							bs_mutations[pdb+"_"+str(chain)].append(mut+":"+pdbe_resi)
							print(ec_num+"\t"+uniprot_id+"\t"+orgn_name+"\t"+struct_type+"\t"+str(pdb)+"\t"+str(chain)+"\t"+mut+"\t"+pdbe_resi+"\tBinding site", file=out)
						else:
							non_bs_mutations[pdb+"_"+str(chain)].append(mut+":"+pdbe_resi)
							print(ec_num+"\t"+uniprot_id+"\t"+orgn_name+"\t"+struct_type+"\t"+str(pdb)+"\t"+str(chain)+"\t"+mut+"\t"+pdbe_resi+"\tNon-binding site", file=out)
					else:
						print(pdb+" "+str(chain)+" "+mut+" could not be mapped to the PDB structure.")				
		
	else:  #Structures mapped for cofactor/apo cases could not be resolved
		print(ec_num+" "+uniprot_id+" "+orgn_name+" cannot be resolved for Cofactor/Apo only cases.")
		continue
		
	#break




























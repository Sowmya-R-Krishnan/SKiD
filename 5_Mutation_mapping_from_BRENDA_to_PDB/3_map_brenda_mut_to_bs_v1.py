#Program to map BRENDA mutations to UniProt sequences and PDB structure residue numbering
#PDB to UniProt residue numbering is taken from SIFTS database (EMBL EBI)

import os
import sys
import csv
import pandas as pd
import numpy as np
from collections import defaultdict
import xml.etree.ElementTree as ET

data = sys.argv[1]
bs_data = sys.argv[2]
sifts_path = sys.argv[3]
outfile = sys.argv[4]

df = pd.read_csv(data, sep="\t", header=0)
#print(df.columns)  #['EC_number', 'Substrate', 'Protein_details', 'PDB_IDs', 'Organism_name', 'kcat_value', 'Mutant', 'Mutation', 'pH', 'Temperature', 'References', 'Site_type', 'Substrate_SMILES', 'Mol_ID']

bs_df = pd.read_csv(bs_data, sep="\t", header=None)

def get_sifts_mapping(pdb_id, sifts_path):
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
					if(sub_elements.tag.endswith("crossRefDb") and sub_elements.attrib["dbSource"]=="UniProt"):
						uniprot_id = sub_elements.attrib["dbAccessionId"]
						uniprot_resinum = sub_elements.attrib["dbResNum"]
						uniprot_resiname = sub_elements.attrib["dbResName"]
						uniprot_resid = uniprot_resiname+"_"+str(uniprot_resinum)
						
						#mapped_resi = parent_map[sub_elements]
						parent = parent_map[sub_elements]
						mapped_resi = [ele for ele in parent.iter() if(ele.tag.endswith("crossRefDb") and ele.attrib["dbSource"]=="PDB")]
						pdbe_resinum = mapped_resi[0].attrib["dbResNum"]
						pdbe_resiname = mapped_resi[0].attrib["dbResName"]
						pdbe_resid = pdbe_resiname+"_"+str(pdbe_resinum)
						
						resi_map[uniprot_resid] = pdbe_resid
				
				chain_mapping[chain_id] = resi_map
				
		return chain_mapping  #XXX: A collections defaultdict object with SIFTS PDB -> UniProt residue numbering
	except:
		return {}
		
three_to_one = {"GLY":"G", "ALA":"A", "VAL":"V", "LEU":"L", "ILE":"I", "PRO":"P", "MET":"M", "SER":"S", "TYR":"Y", "TRP":"W", "PHE":"F", "THR":"T", "HIS":"H", "LYS":"K", "ARG":"R", "ASP":"D", "ASN":"N", "GLU":"E", "GLN":"Q", "CYS":"C"}
one_to_three = dict((c,i) for i, c in three_to_one.items())

bs_mutants = 0
non_bs_mutants = 0

out = open(outfile, "w")
print("EC_number\tSubstrate\tProtein_details\tPDB_IDs\tOrganism_name\tKm_value\tMutant\tMutation\tpH\tTemperature\tReferences\tSite_type\tSubstrate_SMILES\tMol_ID", file=out)

for i, row in df.iterrows():
	if(row["Mutant"]=="yes"):
		#print(row["EC_number"], row["Organism_name"])
		pdb_list = row["PDB_IDs"].split(",")
		mutation = row["Mutation"].split("/")  #XXX: Separate multi-mutation cases into individual mutations
		
		pdbe_resimap = defaultdict(list)
		mutant = 0
		for j, pdb in enumerate(pdb_list):
			resi_map = get_sifts_mapping(pdb, sifts_path)
			for mut in mutation:
				try:
					wt_resi = mut[0]
					mut_resi = mut[-1]
					mut_resinum = int(mut[1:-1])
					#print(mut_resi, wt_resi, mut_resinum)
				except:
					print(row["EC_number"], row["Organism_name"], row["Protein_details"], pdb, mut, row["References"])
					continue
				
				#TODO: First chain iterator; Check residue mapping with UniProt; Check UniProt mapping to PDB; Check presence of residue in binding site; Prepare stats
				#print(pdb, resi_map.keys())
				for chain, resnum_dict in resi_map.items():
					resnum_dict_inverted = dict((c,i) for i,c in resnum_dict.items())
					try:
						if(wt_resi+"_"+str(mut_resinum) in resnum_dict):
							pdbe_resi = resnum_dict[wt_resi+"_"+str(mut_resinum)]
							pdbe_resimap[pdb].append(pdbe_resi+"_"+str(chain))
							check_res = pdbe_resi.replace("_", "")  #BRENDA mutation comes from UniProt
						elif(one_to_three[wt_resi]+"_"+str(mut_resinum) in resnum_dict_inverted):
							pdbe_resi = resnum_dict_inverted[one_to_three[wt_resi]+"_"+str(mut_resinum)]
							pdbe_resimap[pdb].append(pdbe_resi+"_"+str(chain))
							resparts = pdbe_resi.split("_")
							check_res = one_to_three[resparts[0]]+str(resparts[1]) #BRENDA mutation comes from PDB
						else:
							print(row["EC_number"], row["Organism_name"], row["Protein_details"], pdb, chain, mut, row["References"])
							continue
						
							bs_resi = bs_df[bs_df[0]==str(pdb)+"_"+str(chain)].iloc[0][1]
							#print(check_res, bs_resi)
							#print(pdb, mut, pdbe_resi)
							if(check_res in bs_resi):
								mutant = mutant + 1
					except:
						print(row["EC_number"], row["Organism_name"], row["Protein_details"], pdb, chain, mut, row["References"])
						continue
		
		if(mutant > 0):
			bs_mutants = bs_mutants + 1
			contents = row.to_list()
			contents_final = [str(x) for x in contents]
			print("\t".join(contents_final), file=out)
		else:
			non_bs_mutants = non_bs_mutants + 1
		
		#if(i>500):	
			#break
			
	#break

print(bs_mutants, non_bs_mutants, bs_mutants+non_bs_mutants)






























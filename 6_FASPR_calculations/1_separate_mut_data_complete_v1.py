#Separate binding site mutations and complete PDB structures for final analysis

import sys
import csv
import pandas as pd
import numpy as np
import ast
from ast import literal_eval

mut_mapped_data = sys.argv[1]
data_all = sys.argv[2]
bs_data = sys.argv[3]
cofactor_resolve_df = sys.argv[4]
apo_resolve_df = sys.argv[5]
wt_mutant_data = sys.argv[6]
outfile1 = sys.argv[7]
outfile2 = sys.argv[8]
outfile3 = sys.argv[9]

mut_mapped_df = pd.read_csv(mut_mapped_data, sep="\t", header=0)
#print(mut_mapped_df.columns)  #['EC_number', 'UniProt_ID', 'Organism_name', 'Structure_type', 'PDB_ID', 'Chain_ID', 'BRENDA_mutation', 'PDB_mutation', 'Mutation_type']

brenda_raw_df = pd.read_csv(data_all, sep="\t", header=0)
#print(brenda_raw_df.columns)  #['EC_number', 'Substrate', 'Protein_details', 'PDB_IDs', 'Organism_name', 'kcat_value', 'Mutant', 'Mutation', 'pH', 'Temperature', 'References', 'Site_type', 'Substrate_SMILES', 'Mol_ID']

bs_df = pd.read_csv(bs_data, sep="\t", header=None)

cofactor_df = pd.read_csv(cofactor_resolve_df, sep="\t", header=0)  #XXX: Consider only proximal cofactor sites and not distal sites
#print(cofactor_df.columns)  #[EC_number, UniProt_ID, Organism_name, Cofactor_only, Substrate_site, Reference_PDB, Comments]
cofactor_sub_df = cofactor_df[cofactor_df["Substrate_site"]=="Proximal"]

apo_df = pd.read_csv(apo_resolve_df, sep="\t", header=0)
#print(apo_df.columns)  #[EC_number, UniProt_ID, Organism_name, Apo_structures, Holo_mapping]

wt_mutant_df = pd.read_csv(wt_mutant_data, sep="\t", header=0)
#print(wt_mutant_df.columns)  #['EC_number', 'UniProt_ID', 'Organism_name', 'WT_structures', 'Mutant_structures', 'Structure_type']

def extract_uniprot(pdb_comment):
	try:
		prot_dict = literal_eval(pdb_comment)
		return prot_dict["accessions"][0]
	except:
		return pdb_comment

#-------------------------------------------------------------------------------------------------------------------------------------------
out3 = open(outfile3, "w")
print('EC_number\tSubstrate\tUniProt_ID\tWT_PDB\tMutant_PDB\tMutant_PDB_only\tOrganism_name\tKm_value\tMutant\tMutation\tpH\tTemperature\tReferences\tSite_type\tSubstrate_SMILES\tMol_ID', file=out3)

out1 = open(outfile1, "w")
print('EC_number\tSubstrate\tUniProt_ID\tWT_PDB\tWT_mutmap\tMutant_PDB\tMut_mutmap\tMutant_PDB_only\tOrganism_name\tKm_value\tMutant\tMutation\tpH\tTemperature\tReferences\tSite_type\tSubstrate_SMILES\tMol_ID', file=out1)

out2 = open(outfile2, "w")
print('EC_number\tSubstrate\tUniProt_ID\tWT_PDB\tWT_mutmap\tMutant_PDB\tMut_mutmap\tMutant_PDB_only\tOrganism_name\tKm_value\tMutant\tMutation\tpH\tTemperature\tReferences\tSite_type\tSubstrate_SMILES\tMol_ID', file=out2)

for i, row in brenda_raw_df.iterrows():
	ec_num = row["EC_number"]
	orgn_name = row["Organism_name"]
	pdb_list = row["PDB_IDs"].split(",")
	uniprot_id = extract_uniprot(row["Protein_details"])
	
	if(row["Mutant"]=="yes"):
		sub_mut_df = mut_mapped_df[(mut_mapped_df["EC_number"]==ec_num) & (mut_mapped_df["UniProt_ID"]==uniprot_id) & (mut_mapped_df["Organism_name"]==orgn_name)]
		
		if(len(sub_mut_df.index)>0):
			mutations = row["Mutation"].split("/")  #To handle cases with multiple mutations in same entry
		
			for mut_target in mutations:
				bs_flag = 0  #Default 0: Non-binding site mutation
				mutant_flag = "no"
				best_wt_pdb = "-"
				best_mut_pdb = "-"
				wt_mutmap = "-"
				mut_mutmap = "-"
			
				wt_structs = []
				mut_structs = []
				
				mutant_info = wt_mutant_df[(wt_mutant_df["EC_number"]==ec_num) & (wt_mutant_df["UniProt_ID"]==uniprot_id) & (wt_mutant_df["Organism_name"]==orgn_name)]
				if(len(mutant_info.index)>0):
					wt_structs = mutant_info["WT_structures"].item().split(",")
					mut_structs = mutant_info["Mutant_structures"].item().split(",")
					
					wt_structs_processed = []
					mut_structs_processed = []
					if(mutant_info["Structure_type"].item()=="Cofactor only"):
						cof_row = cofactor_sub_df[(cofactor_sub_df["EC_number"]==ec_num) & (cofactor_sub_df["UniProt_ID"]==uniprot_id) & (cofactor_sub_df["Organism_name"]==orgn_name)]
						if(len(cof_row.index)>0):
							wt_structs_processed.extend(wt_structs)
							mut_structs_processed.extend(mut_structs)
						else:
							wt_structs_processed.append("-")
							mut_structs_processed.append("-")
					elif(mutant_info["Structure_type"].item()=="Apo only"):
						apo_row = apo_df[(apo_df["EC_number"]==ec_num) & (apo_df["UniProt_ID"]==uniprot_id) & (apo_df["Organism_name"]==orgn_name)]
						if(len(apo_row.index)>0):
							wt_structs_processed.extend(wt_structs)
							mut_structs_processed.extend(mut_structs)
						else:
							wt_structs_processed.append("-")
							mut_structs_processed.append("-")
					else:
						wt_structs_processed.extend(wt_structs)
						mut_structs_processed.extend(mut_structs)
							
					if(wt_structs_processed[0]!="-"):
						for l, wt in enumerate(wt_structs_processed):
							mut_match = sub_mut_df[(sub_mut_df["PDB_ID"]==wt) & (sub_mut_df["BRENDA_mutation"]==mut_target)]
							if(len(mut_match.index)>0):
								if(mut_match.iloc[0]["Mutation_type"]=="Binding site"):
									bs_flag = 1
									
								best_wt_pdb = wt
								wt_mutmap = mut_match.iloc[0]["PDB_mutation"]
								break
								
						if(mut_structs_processed[0]!="-"):
							for l, mut in enumerate(mut_structs_processed):
								mut_match = sub_mut_df[(sub_mut_df["PDB_ID"]==mut) & (sub_mut_df["BRENDA_mutation"]==mut_target)]
								if(len(mut_match.index)>0):
									if(mut_match.iloc[0]["Mutation_type"]=="Binding site"):
										bs_flag = 1

									best_mut_pdb = mut
									mut_mutmap = mut_match.iloc[0]["PDB_mutation"]
									break
					else:
						mutant_flag = "yes"
						for l, mut in enumerate(mut_structs_processed):
							mut_match = sub_mut_df[(sub_mut_df["PDB_ID"]==mut) & (sub_mut_df["BRENDA_mutation"]==mut_target)]
							if(len(mut_match.index)>0):
								if(mut_match.iloc[0]["Mutation_type"]=="Binding site"):
									bs_flag = 1

								best_mut_pdb = mut
								mut_mutmap = mut_match.iloc[0]["PDB_mutation"]
								break
								
					print(ec_num, uniprot_id, orgn_name, bs_flag, mutant_flag, mut_target, best_wt_pdb, wt_mutmap, best_mut_pdb, mut_mutmap)				
					if(bs_flag>0):
						print(ec_num+"\t"+row["Substrate"]+"\t"+uniprot_id+"\t"+str(best_wt_pdb)+"\t"+str(wt_mutmap)+"\t"+str(best_mut_pdb)+"\t"+str(mut_mutmap)+"\t"+mutant_flag+"\t"+orgn_name+"\t"+str(row["Km_value"])+"\t"+row["Mutant"]+"\t"+row["Mutation"]+"\t"+str(row["pH"])+"\t"+str(row["Temperature"])+"\t"+str(row["References"])+"\t"+row["Site_type"]+"\t"+row["Substrate_SMILES"]+"\t"+str(row["Mol_ID"]), file=out1)
					else:
						print(ec_num+"\t"+row["Substrate"]+"\t"+uniprot_id+"\t"+str(best_wt_pdb)+"\t"+str(wt_mutmap)+"\t"+str(best_mut_pdb)+"\t"+str(mut_mutmap)+"\t"+mutant_flag+"\t"+orgn_name+"\t"+str(row["Km_value"])+"\t"+row["Mutant"]+"\t"+row["Mutation"]+"\t"+str(row["pH"])+"\t"+str(row["Temperature"])+"\t"+str(row["References"])+"\t"+row["Site_type"]+"\t"+row["Substrate_SMILES"]+"\t"+str(row["Mol_ID"]), file=out2)
			
			#break
		else:  #Mutation could not be mapped to the PDB structure due to several errors identified in BRENDA
			continue
			
	else:
		mutant_flag = "no"
		best_wt_pdb = "-"
		best_mut_pdb = "-"
		
		mutant_info = wt_mutant_df[(wt_mutant_df["EC_number"]==ec_num) & (wt_mutant_df["UniProt_ID"]==uniprot_id) & (wt_mutant_df["Organism_name"]==orgn_name)]
		if(len(mutant_info.index)>0):
			wt_structs = mutant_info["WT_structures"].item().split(",")
			mut_structs = mutant_info["Mutant_structures"].item().split(",")
			
			wt_structs_processed = []
			mut_structs_processed = []
			if(mutant_info["Structure_type"].item()=="Cofactor only"):
				cof_row = cofactor_sub_df[(cofactor_sub_df["EC_number"]==ec_num) & (cofactor_sub_df["UniProt_ID"]==uniprot_id) & (cofactor_sub_df["Organism_name"]==orgn_name)]
				if(len(cof_row.index)>0):
					wt_structs_processed.extend(wt_structs)
					mut_structs_processed.extend(mut_structs)
				else:
					wt_structs_processed.append("-")
					mut_structs_processed.append("-")
			elif(mutant_info["Structure_type"].item()=="Apo only"):
				apo_row = apo_df[(apo_df["EC_number"]==ec_num) & (apo_df["UniProt_ID"]==uniprot_id) & (apo_df["Organism_name"]==orgn_name)]
				if(len(apo_row.index)>0):
					wt_structs_processed.extend(wt_structs)
					mut_structs_processed.extend(mut_structs)
				else:
					wt_structs_processed.append("-")
					mut_structs_processed.append("-")
			else:
				wt_structs_processed.extend(wt_structs)
				mut_structs_processed.extend(mut_structs)				
			
			#print(wt_structs, mut_structs)
			if(wt_structs_processed[0]!="-"):
				best_wt_pdb = wt_structs_processed[0]
				if(mut_structs_processed[0]!="-"):
					best_mut_pdb = mut_structs_processed[0]
			else:
				mutant_flag = "yes"
				best_mut_pdb = mut_structs_processed[0]
			
			#XXX: Any one structure should be available to use the data!
			if(best_wt_pdb!="-" or best_mut_pdb!="-"):	
				print(ec_num+"\t"+row["Substrate"]+"\t"+uniprot_id+"\t"+str(best_wt_pdb)+"\t"+str(best_mut_pdb)+"\t"+mutant_flag+"\t"+orgn_name+"\t"+str(row["Km_value"])+"\t"+row["Mutant"]+"\t"+row["Mutation"]+"\t"+str(row["pH"])+"\t"+str(row["Temperature"])+"\t"+str(row["References"])+"\t"+row["Site_type"]+"\t"+row["Substrate_SMILES"]+"\t"+str(row["Mol_ID"]), file=out3)	






































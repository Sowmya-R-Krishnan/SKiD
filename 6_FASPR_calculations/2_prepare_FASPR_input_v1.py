#Program to extract mutation datasets for running FASPR automation code

import sys
import csv
import pandas as pd
import numpy as np
from collections import defaultdict

infile = sys.argv[1]
mutmap_file = sys.argv[2]
bs_data = sys.argv[3]
outfile = sys.argv[4]

mut_mapped_df = pd.read_csv(infile, sep="\t", header=0)
#print(mut_mapped_df.columns)  #['EC_number', 'Substrate', 'UniProt_ID', 'WT_PDB', 'WT_mutmap', 'Mutant_PDB', 'Mut_mutmap', 'Mutant_PDB_only', 'Organism_name', 'kcat_value', 'Mutant', 'Mutation', 'pH', 'Temperature', 'References', 'Site_type', 'Substrate_SMILES', 'Mol_ID']

bs_df = pd.read_csv(bs_data, sep="\t", header=None)

three_to_one = {"ALA":"A", "GLY":"G", "VAL":"V", "LEU":"L", "ILE":"I", "PHE":"F", "TRP":"W", "TYR":"Y", "HIS":"H", "ARG":"R", "LYS":"K", "ASN":"N", "GLN":"Q", "ASP":"D", "GLU":"E", "SER":"S", "CYS":"C", "THR":"T", "MET":"M", "PRO":"P", "OCS":"C", "CSO":"C", "MSE":"M", "MHS":"H", "CME":"C", "CSD":"C", "TSY":"C", "CSS":"C", "SCY":"C", "PTR":"Y", "DDZ":"A", "JJJ":"C", "MHO":"M", "TPO":"T", "SEP":"S", "KPI":"K", "SLZ":"K", "KPF":"K", "KPY":"K", "74P":"K", "VPV":"K", "KGC":"K", "KYQ":"K", "LYF":"K", "CIR":"R", "SNC":"C", "CSP":"C", "LLP":"K", "KCX":"K"}

mutmap_df = pd.read_csv(mutmap_file, sep="\t", header=0)

out = open(outfile, "w")
print("affin.pdb_id\taffin.chain\tmutation", file=out)

mutmap = defaultdict(dict)

'''
for i, row in mut_mapped_df.iterrows():
	wt_pdb = row["WT_PDB"]
	mut_resi = row["WT_mutmap"]
	mutation = row["Mutation"]  #Mutation provided in BRENDA
	
	if(wt_pdb=="-"):
		wt_pdb = row["Mutant_PDB"]
		mut_resi = row["Mut_mutmap"]  #PDB-mapped mutation

	muts = []
	muts = mutation.split("/")
	final_mut = ""
	for mut in muts:
		#print(mut, mut_resi)
		try:
			mres = mut_resi.split("_")[0]
			resid = mut_resi.split("_")[1]
			
			if(len(mres)==1 and resid==mut[1:-1]):  #If single-letter code has been used and residue first letters match
				final_mut = mres+str(resid)+mut[-1]
				break
			elif(resid==mut[1:-1]):  #If three-letter code has been used and residue numbers match
				final_mut = three_to_one[mres]+str(resid)+mut[-1]
				break
			else:
				continue
		except:
			continue  #To handle cases where the mutation columns are empty => Wrong data curation
			
	#print(wt_pdb, final_mut)
	mutmap[wt_pdb][mut] = final_mut
	#break
'''

#EC_number, UniProt_ID, Organism_name, Structure_type, PDB_ID, Chain_ID, BRENDA_mutation, PDB_mutation, Mutation_type
for i, row in mutmap_df.iterrows():
	pdb = row["PDB_ID"]
	brenda_mut = row["BRENDA_mutation"]
	PDB_mut = row["PDB_mutation"]
	
	if(len(PDB_mut.split("_")[0])>1):
		pdb_mut_new = three_to_one[PDB_mut.split("_")[0]]+PDB_mut.split("_")[1]+brenda_mut[-1]
	else:
		pdb_mut_new = PDB_mut.split("_")[0]+PDB_mut.split("_")[1]+brenda_mut[-1]
		
	mutmap[pdb][brenda_mut] = pdb_mut_new
	
#print(mutmap["5AG0"])

done = []
for i, row in mut_mapped_df.iterrows():
	wt_pdb = row["WT_PDB"]
	mut_resi = row["WT_mutmap"]
	mutation = row["Mutation"]			
	
	if(wt_pdb=="-"):
		wt_pdb = row["Mutant_PDB"]
		mut_resi = row["Mut_mutmap"]
	
	try:	
		match_df = bs_df[bs_df[0].str.startswith(wt_pdb)]
		chain_id = match_df.iloc[0][0].split("_")[1]
			
		muts = mutation.split("/")
		final_muts = []
		mutstr = ""
		for mut in muts:
			try:
				final_mut = mutmap[wt_pdb][mut]
				if(final_mut not in final_muts):
					final_muts.append(final_mut)
			except:
				print(wt_pdb, chain_id, mut)
				continue
				
		mutstr = "/".join(final_muts)

		if((wt_pdb, chain_id, mutstr) not in done):
			print(str(wt_pdb)+"\t"+str(chain_id)+"\t"+mutstr, file=out)
			done.append((wt_pdb, chain_id, mutstr))
		else:
			continue
	except:
		continue  #IndexError: single positional indexer is out-of-bounds
	
	#break




























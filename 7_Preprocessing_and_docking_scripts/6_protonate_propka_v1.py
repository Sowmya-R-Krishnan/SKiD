#Protonate the protein based on pH from BRENDA using PDB2PQR and PROPKA tools

import os
import sys
import csv
import pandas as pd
import numpy as np
import pymol
from pymol import cmd

infile = sys.argv[1]
km_final_file = sys.argv[2]
pdbpath = sys.argv[3]
outpath = sys.argv[4]
logpath = sys.argv[5]

df = pd.read_csv(infile, sep="\t", header=0)
#print(df.columns)  #['EC_number', 'Substrate', 'UniProt_ID', 'Protein_file', 'Organism_name', 'Km_value', 'Mutant', 'Mutation', 'pH', 'Temperature', 'References', 'Site_type', 'Substrate_SMILES', 'Mol_file']

km_df = pd.read_csv(km_final_file, sep=" ", header=None)
fnames = os.listdir(pdbpath)

done = []
for i, row in df.iterrows():
	prot_file = row["Protein_file"]
	out_file = prot_file.replace(".pdb","").replace("_mut","")  #Covers both mutant and wild-type structures
	pdb_id = prot_file[0:4]
	if(pdb_id in list(km_df[0]) and out_file not in done):  #If the PDB ID is present among the list of valid structures in Km dataset, perform protonation
		try:
			ph_val = float(row["pH"])
			
			if(out_file+"_prepped.pdb" in fnames and not os.path.exists(outpath+out_file+"_"+str(ph_val)+".pdb")):
				prot_file = out_file+"_prepped.pdb"
				#Remove ligand atoms and save the protein file
				het_records = []
				prot_records = []
				with open(pdbpath+prot_file) as f:
					for line in f.readlines():
						line = line.strip()
						if(line.startswith("HETATM")):
							het_records.append(line)
						else:
							prot_records.append(line)
							
				temp_out = open(outpath+"tmp.pdb", "w")
				for line in prot_records:
					print(line, file=temp_out)
				temp_out.close()
				
				#Run PDB2PQR and PROPKA tools for protonating as per pH
				os.system("pdb2pqr "+outpath+"tmp.pdb "+outpath+"tmp.pqr --ff=CHARMM --pdb-output="+outpath+"tmp_1.pdb --with-ph="+str(ph_val)+" --titration-state-method=propka -o "+str(ph_val)+" -q")
				os.system("rm "+outpath+"tmp.pdb")
				os.system("rm "+outpath+"tmp.pqr")
				os.system("mv "+outpath+"tmp.log "+logpath+out_file+"_"+str(ph_val)+".log")
				
				#Append substrate and cofactor atoms back to the protonated PDB file
				mod_prot = open(outpath+out_file+"_"+str(ph_val)+".pdb", "w")
				with open(outpath+"tmp_1.pdb") as f:
					for line in f.readlines():
						line = line.strip()
						if(line!="END"):
							print(line, file=mod_prot)
						else:
							for l in het_records:
								print(l, file=mod_prot)
							print("END", file=mod_prot)
							
				mod_prot.close()
				os.system("rm "+outpath+"tmp_1.pdb")
				done.append(out_file)
							
			else:
				continue			
		except:
			ph_val = "nil"
			if(out_file+"_prepped.pdb" in fnames):
				prot_file = out_file+"_prepped.pdb"
				cmd.load(pdbpath+prot_file, "ref_chain")
				cmd.h_add("ref_chain")
				cmd.save(outpath+out_file+"_"+str(ph_val)+".pdb")
				cmd.reinitialize()
				done.append(out_file)
			else:
				continue
				
		print(out_file)
		
	else:
		continue





























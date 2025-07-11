#Resolve errors in FASPR calculations

import os
import sys
import csv
import pandas as pd
import numpy as np

infile = sys.argv[1]
pdbpath = sys.argv[2]
seq_path = sys.argv[3]
outpath = sys.argv[4]

df = pd.read_csv(infile, sep="\t", header=0)
#print(df.columns)  #['affin.pdb_id', 'mutation', 'affin.chain', 'Reason', 'Resolution', 'Type']

for i, row in df.iterrows():
	reasons = row["Reason"]
	
	if(row["Type"]=="Automated" and reasons.startswith("Chain ID")):
		mod_chain = row["Resolution"][-1]
		#print(mod_chain)
		
		v1 = row["affin.pdb_id"]
		if(not os.path.exists(pdbpath+v1+".pdb")):
			v1 = str.lower(v1)
			
		v2 = pdbpath
		v4 = seq_path
		mutation_list = str.split(row["mutation"], '/')
		v3 = mod_chain + ":" + mutation_list[0]
		for j in range(1, len(mutation_list)):
			v3 = v3 + "+" + mod_chain + ':' + mutation_list[j]
		
		os.system('python3 separate_chains_and_seqs_v2.py '+v1+' '+v2+' '+v3+' '+v4)
		os.system("pdbfixer "+seq_path+v1+"_"+mod_chain+".pdb --add-atoms=heavy --output="+seq_path+v1+"_"+mod_chain+".pdb ")
		
		os.system('./FASPR-master/FASPR -i ./FASPR_IP_kcat_mutant/'+v1+'_'+mod_chain+'.pdb -s ./FASPR_IP_kcat_mutant/'+v1+'_'+mod_chain+'_'+v3+'_mut.fasta -o ./FASPR_OP_kcat_mutant/'+v1+'_'+mod_chain+'_'+v3+'_mut.pdb')
		#print(v1, mod_chain)
		
	elif(row["Type"]=="Automated" and reasons=="Wrong mutation"):
		v1 = row["affin.pdb_id"]
		if(not os.path.exists(pdbpath+v1+".pdb")):
			v1 = str.lower(v1)
			
		v2 = pdbpath
		v4 = seq_path
		mutation_list = str.split(row["Resolution"], '/')
		v3 = row['affin.chain'] + ":" + mutation_list[0]
		for j in range(1, len(mutation_list)):
			v3 = v3 + "+" + row['affin.chain'] + ':' + mutation_list[j]
		
		os.system('python3 separate_chains_and_seqs_v2.py '+v1+' '+v2+' '+v3+' '+v4)
		os.system("pdbfixer "+seq_path+v1+"_"+row['affin.chain']+".pdb --add-atoms=heavy --output="+seq_path+v1+"_"+row['affin.chain']+".pdb ")
		
		os.system('./FASPR-master/FASPR -i ./FASPR_IP_kcat_mutant/'+v1+'_'+row['affin.chain']+'.pdb -s ./FASPR_IP_kcat_mutant/'+v1+'_'+row['affin.chain']+'_'+v3+'_mut.fasta -o ./FASPR_OP_kcat_mutant/'+v1+'_'+row['affin.chain']+'_'+v3+'_mut.pdb')
		
	elif(row["Type"]=="Automated" and reasons.startswith("Backbone incomplete")):		
		v1 = row["affin.pdb_id"]
		if(not os.path.exists(pdbpath+v1+".pdb")):
			v1 = str.lower(v1)
			
		v2 = pdbpath
		v4 = seq_path
		mutation_list = str.split(row["mutation"], '/')
		v3 = row['affin.chain'] + ":" + mutation_list[0]
		for j in range(1, len(mutation_list)):
			v3 = v3 + "+" + row['affin.chain'] + ':' + mutation_list[j]
		
		os.system('python3 separate_chains_and_seqs_v2.py '+v1+' '+v2+' '+v3+' '+v4)
		os.system("pdbfixer "+seq_path+v1+"_"+row['affin.chain']+".pdb --add-atoms=heavy --output="+seq_path+v1+"_"+row['affin.chain']+".pdb ")
		
		os.system('./FASPR-master/FASPR -i ./FASPR_IP_kcat_mutant/'+v1+'_'+row['affin.chain']+'.pdb -s ./FASPR_IP_kcat_mutant/'+v1+'_'+row['affin.chain']+'_'+v3+'_mut.fasta -o ./FASPR_OP_kcat_mutant/'+v1+'_'+row['affin.chain']+'_'+v3+'_mut.pdb')
		
	else:
		continue
			






























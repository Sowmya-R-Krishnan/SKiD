#Program to modify HETATM records of non-canonical residues in PDB files

import os
import sys
import csv
import re
import pandas as pd
import numpy as np

infile = sys.argv[1]
inpath = sys.argv[2]
outpath = sys.argv[3]

df = pd.read_csv(infile, sep="\t", header=0)

replace_dict = {"OCS":"CYS", "CSO":"CYS", "MSE":"MET", "MHS":"HIS", "CME":"CYS", "CSD":"CYS", "TSY":"CYS", "CSS":"CYS", "SCY":"CYS", "PTR":"TYR", "DDZ":"ALA", "JJJ":"CYS", "MHO":"MET", "TPO":"THR", "SEP":"SER", "KPI":"LYS", "SLZ":"LYS", "KPF":"LYS", "KPY":"LYS", "74P":"LYS", "VPV":"LYS", "KGC":"LYS", "KYQ":"LYS", "LYF":"LYS", "CIR":"ARG", "SNC":"CYS", "CSP":"CYS", "LLP":"LYS", "KCX":"LYS"}

for i, row in df.iterrows():
	if(row["Reason"]=="Residue is non-canonical"):
		mut = row["mutation"]
		chain = row["affin.chain"]
		pdb = row["affin.pdb_id"]
		
		muts = mut.split("/")
		
		contents = []
		with open(inpath+pdb+"_"+chain+".pdb") as pfile:
			for line in pfile.readlines():
				line = line.strip()
				if(line.startswith("HETATM")):
					for m in muts:
						if(int(line[22:26].strip())==int(m[1:-1])):
							if(line[12:16].strip() not in ["N", "CA", "C", "O"]):
								continue
							else:
								try:
									mod_res = line[17:20].strip()
									line = line.replace(mod_res, replace_dict[mod_res])
									line = line.replace("HETATM", "ATOM  ")
									contents.append(line)
								except:
									line = line.replace("HETATM", "ATOM  ")
									contents.append(line)
								#print(line)
								
				else:
					contents.append(line)
					
		out = open(outpath+pdb+"_"+chain+".pdb", "w")			
		for line in contents:
			print(line, file=out)	
			
		v3 = chain + ":" + muts[0]
		for j in range(1, len(muts)):
        		v3 = v3 + "+" + chain + ':' + muts[j]
        	
		os.system('python3 separate_chains_and_seqs_v2.py '+pdb+"_"+chain+' '+outpath+' '+v3+' '+outpath)
		os.system("pdbfixer "+outpath+pdb+"_"+chain+"_"+chain+".pdb --add-atoms=heavy --output="+outpath+pdb+"_"+chain+"_"+chain+".pdb ")
		os.system('./FASPR-master/FASPR -i ./FASPR_IP_Km_mutant/'+pdb+'_'+chain+'_'+chain+'.pdb -s ./FASPR_IP_Km_mutant/'+pdb+'_'+chain+'_'+chain+'_'+v3+'_mut.fasta -o ./FASPR_OP_Km_mutant/'+pdb+'_'+chain+'_'+chain+'_'+v3+'_mut.pdb')
			
		print(pdb)
		#break































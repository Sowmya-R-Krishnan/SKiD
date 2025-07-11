#Program to prepare all structures with a substrate and cofactor (WT and mutant) for docking pre-processing

import os
import sys
import csv
import pandas as pd
import numpy as np
import Bio
import ast
from ast import literal_eval
from Bio import PDB
from Bio.PDB import PDBIO, Select, PDBParser, NeighborSearch
import warnings
import pymol
from pymol import cmd

infile = sys.argv[1]
sub_cof_mapping = sys.argv[2]
refpath = sys.argv[3]
mutpath = sys.argv[4]
outpath = sys.argv[5]
outfile = sys.argv[6]

df = pd.read_csv(infile, sep="\t", header=0)
#print(df.columns)  #['Protein_file', 'Ligand_file', 'Data_type']

sub_cof_df = pd.read_csv(sub_cof_mapping, sep=" ", header=None)

wt_fnames = os.listdir(refpath)
mut_fnames = os.listdir(mutpath)

out = open(outfile, "w")
print("Protein_file\tReference_ligand\tLigand_file\tData_type", file=out)

for i, row in df.iterrows():
	#print(row)
	pdb_file = row["Protein_file"]
	pdb_id = pdb_file[0:4]
	
	final_fname = ""
	ref_sub = ""
	
	if(row["Data_type"]=="WT"):
		found = 0
		for n in wt_fnames:
			if(n.startswith(pdb_id)):
				found = 1
				final_fname = n
				break
				
		if(found==1):
			try:
				ref_df = sub_cof_df[sub_cof_df[0]==pdb_id]
				ref_sub = ref_df.iloc[0][2]
				os.system("cp "+refpath+final_fname+" "+outpath)
				print(final_fname+"\t"+str(ref_sub)+"\t"+row["Ligand_file"]+"\tWT", file=out)
				print(final_fname+"\t"+str(ref_sub)+"\t"+row["Ligand_file"]+"\tWT")
			except:
				continue
	else:
		found = 0
		for n in wt_fnames:
			if(n.startswith(pdb_id)):
				found = 1
				final_fname = n
				break
				
		if(found==1):
			try:
				ref_df = sub_cof_df[sub_cof_df[0]==pdb_id]
				ref_sub = ref_df.iloc[0][2]
				
				het_list = []
				with warnings.catch_warnings():
					warnings.simplefilter('ignore')
					pdb = PDBParser().get_structure(pdb_id, refpath+final_fname)
					io = PDBIO()
					io.set_structure(pdb)
					
					for model in pdb:
						for chain in model:
							for resi in chain:
								if(resi.id[0].startswith("H_")):
									het_list.append(resi.resname)

				#---------------------------------------PYMOL ALIGNMENT AND SAVING COORDS---------------------------------#
				cmd.load(refpath+final_fname, "ref_chain")
				cmd.load(mutpath+pdb_file, "align_chain")
				
				cmd.align("align_chain", "ref_chain")
				
				for j, idx in enumerate(het_list):
					sel_cmd = ["ref_chain"]
					sel_cmd.append("resn "+str(idx))
					sel_txt = " and ".join(sel_cmd)
					cmd.select(sel_txt)  #Extract cofactor atoms post alignment
					cmd.extract("het_atoms_"+str(j), "sele")
				
				cmd.delete("ref_chain")  #Delete aligned chain after extraction of cofactor atoms
				cmd.save(outpath+pdb_file)
				cmd.reinitialize()
				#--------------------------------------------------------------------------------------------------------#

				print(pdb_file+"\t"+str(ref_sub)+"\t"+row["Ligand_file"]+"\tMUT", file=out)
				print(pdb_file+"\t"+str(ref_sub)+"\t"+row["Ligand_file"]+"\tMUT")
			except:
				continue
			
	#break





























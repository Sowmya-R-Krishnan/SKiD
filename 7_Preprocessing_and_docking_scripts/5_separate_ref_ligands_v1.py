#Program to extract the substrates into a separate file for autobox option

import os
import sys
import csv
import pandas as pd
import numpy as np
import Bio
import ast
from ast import literal_eval
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, Select
import warnings

infile = sys.argv[1]
inpath = sys.argv[2]
outpath = sys.argv[3]
cof_sub_data = sys.argv[4]
apo_only = sys.argv[5]
cof_only = sys.argv[6]

df = pd.read_csv(infile, sep="\t", header=0)
#Protein_file	Reference_ligand	Ligand_file	Data_type

cof_sub_df = pd.read_csv(cof_sub_data, sep="\t", header=0)
apo_df = pd.read_csv(apo_only, sep="\t", header=0)
cof_df = pd.read_csv(cof_only, sep="\t", header=0)

#Residue selection class
class ResidueSelect(Select):
	def __init__(self, accepted_resi):
		self.accepted_resi = accepted_resi

	def accept_residue(self, resi):
		if(resi in self.accepted_resi):
			return 1
		else:
			return 0

done = []
for i, row in df.iterrows():
	if(row["Protein_file"] not in done):
		#ref_ligand = row["Reference_ligand"].replace("[\'","").replace("\']","")
		ref_ligand = literal_eval(row["Reference_ligand"])
		
		prot_resi = []
		ref_lig = []

		with warnings.catch_warnings():
			warnings.simplefilter('ignore')
			pdb = PDBParser().get_structure(row["Protein_file"], inpath+row["Protein_file"])
			io = PDBIO()
			io.set_structure(pdb)
			
			for model in pdb:
				for chain in model:
					for resi in chain:
						#if(resi.id[0]==ref_ligand):
						if(resi.id[0] in ref_ligand):
							ref_lig.append(resi)
						else:
							prot_resi.append(resi)

		if(len(ref_lig)>0):
			ref_lig_final = []
			
			if(len(ref_lig)>1):
				print(row["Protein_file"])
				cof_data = cof_sub_df[cof_sub_df["PDB_ID"]==row["Protein_file"][0:4]]
				if(cof_data.iloc[0]["Structure_type"]=="Substrate+Cofactor"):
					cof_id = literal_eval(cof_data.iloc[0]["Cofactors"])[0][1]
				elif(cof_data.iloc[0]["Structure_type"]=="Substrate only"):
					cof_id = "-"
				elif(cof_data.iloc[0]["Structure_type"]=="Apo structure"):
					cof_data2 = apo_df[apo_df["PDB_ID"]==row["Protein_file"][0:4]]
					cof_id = cof_data2.iloc[0]["Cofactor_ID"]
				else:
					cof_id = literal_eval(cof_data.iloc[0]["Cofactors"])[0][1]
					
				if(cof_id!="-" and cof_id!=""):
					for resi in prot_resi:
						if(resi.id[0]=="H_"+cof_id):
							com_ref = resi.center_of_mass()
							break
					
					dists = []		
					for resi in ref_lig:
						com_lig = resi.center_of_mass()
						dist = np.linalg.norm(com_ref - com_lig)
						dists.append(dist)
					
					ref_lig_final.append(ref_lig[np.argmin(dists)])  #Closest to cofactor is considered for autobox
				
				else:
					ref_lig_final.append(ref_lig[0])
			else:
				ref_lig_final.append(ref_lig[0])			
		
			outfile = row["Protein_file"].replace("_prep_v1.pdb","").replace("_mut.pdb","")		
			with warnings.catch_warnings():
				warnings.simplefilter('ignore')
				pdb = PDBParser().get_structure(row["Protein_file"], inpath+row["Protein_file"])
				io = PDBIO()
				io.set_structure(pdb)	
				io.save(outpath+outfile+"_orig.pdb", ResidueSelect(ref_lig_final))
				io.save(outpath+outfile+"_prepped.pdb", ResidueSelect(prot_resi))
				
			os.system("obabel -ipdb "+outpath+outfile+"_orig.pdb -osdf -O "+outpath+outfile+"_orig.sdf")
			os.system("rm "+outpath+outfile+"_orig.pdb")
			
			done.append(row["Protein_file"])
			#print(outfile)
		#break






























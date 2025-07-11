#Program to extract binding pockets from the Km and kcat datasets

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

input_data = sys.argv[1]
cofactor = sys.argv[2]
apo = sys.argv[3]
pdbpath = sys.argv[4]  #WT structures
apopath = sys.argv[5]
cofpath = sys.argv[6]
outpath = sys.argv[7]

pdb_df = pd.read_csv(input_data, sep="\t", header=0)
#print(pdb_df.columns)  #['EC_number', 'PDB_ID', 'Structure_type', 'Substrates', 'Cofactors']

cof_df = pd.read_csv(cofactor, sep="\t", header=0)
#print(cof_df.columns)  #['EC_number', 'UniProt_ID', 'Organism_name', 'Cofactor_only', 'Substrate_site', 'Reference_PDB', 'Comments']

apo_df = pd.read_csv(apo, sep="\t", header=0)
#print(apo_df.columns)  #['EC_number', 'UniProt_ID', 'Organism_name', 'Apo_structures', 'Holo_mapping']

#Residue selection class
class ResidueSelect(Select):
	def __init__(self, accepted_resi):
		self.accepted_resi = accepted_resi
		#self.omit_resi = omit_resi
		
		self.ions = ['NCO', 'NA', 'RHD', 'FE', 'SR', 'PO4', 'ZN', 'MG', 'AG', 'K', 'AU3', 'NI', 'SE4', 'ACT', 'IRI', 'IR', 'CL', 'LU', 'TB', 'CA', 'CS', 'TL', 'SO4', 'CU', 'CAC', 'AU', 'BR', 'SM', 'CD', 'MN', 'NH4', 'IR3', 'CO', 'HG', '3CO', 'OS', 'CON', 'BA', 'FE2', 'PB', 'SIN', 'MES', 'GOL', 'MPD', 'BME', '1PE', 'EPE', 'OXY', 'O', 'NO3', 'FES', 'SF4', 'MOS', 'TRS', 'YT3', 'LI', 'DIO', 'MRD', 'CO3', 'F', 'LI']
		self.nc_resi = ["OCS", "CSO", "MSE", "MHS", "CME", "CSD", "TSY", "CSS", "SCY", "PTR", "DDZ", "JJJ", "MHO", "TPO", "SEP", "KPI", "SLZ", "KPF", "KPY", "74P", "VPV", "KGC", "KYQ", "LYF", "CIR", "SNC", "CSP", "LLP", "KCX", "MLY", "MLZ", "OHI", "HTR", "MTY", "OAS", "CSX", "OMT", "TRQ", "CXM", "SCH", "ALY", "TYX", "CAS", "IAS", "SNN", "NEP", "D4B", "2DA", "YCM", "CSR", "HS8", "A0A", "ALS", "TYS", "C3Y", "ASL", "ACE", "DTY", "OMY", "OMZ", "PCA"]

	def accept_residue(self, resi):
		#if (resi in self.accepted_resi and resi.id[0] not in self.omit_resi and resi.resname!="HOH" and resi.resname not in self.ions and resi.resname not in self.nc_resi):
		if(resi in self.accepted_resi):
			return 1
		else:
			return 0

#Function to perform Neighbor Search and extract binding site residues
def extract_bs_resi(pdb_id, chain_id, lig_ids, pdb_path):
	bs_resi = []
	resi_done = []
	
	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		pdb = PDBParser().get_structure(pdb_id, pdb_path+pdb_id+".pdb")
		io = PDBIO()
		io.set_structure(pdb)
		
		for model in pdb:
			for chain in model:
				if(chain.id==chain_id):
					for resi in chain:
						if(resi.resname in lig_ids and resi.resname not in resi_done):
							center_atoms = Bio.PDB.Selection.unfold_entities(resi, 'A')
							alist = Bio.PDB.Selection.unfold_entities(pdb, 'A')
							
							search_center = NeighborSearch(alist)
							selected_resi = {res for center_atom in center_atoms for res in search_center.search(center_atom.coord, 6.5, 'R')}
							bs_resi.extend(selected_resi)
							resi_done.append(resi.resname)  #Only first copy of the substrate is considered
					
			break  #Only first model considered in NMR structures

		bs_resi = list(set(bs_resi))

	return bs_resi

#Function to extract all residues of a chain in PDB structure
def extract_target_chain(pdb_id, chain_id, retain_resi, pdb_path):
	resi_list = []
	het_list = []
	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		pdb = PDBParser().get_structure(pdb_id, pdb_path+pdb_id+".pdb")
		io = PDBIO()
		io.set_structure(pdb)
		
		for model in pdb:
			for chain in model:
				if(chain.id==chain_id):
					for resi in chain:
						if(resi.id[0].startswith("H_") and resi.id[0] in retain_resi):
							resi_list.append(resi)
						elif(not resi.id[0].startswith("H_") and resi.resname!="HOH"):
							resi_list.append(resi)
						else:
							continue			
			break
			
	return resi_list

for i, row in pdb_df.iterrows():
	sub_ids = []
	cof_ids = []

	if(row["Structure_type"]=="Substrate+Cofactor"):
		pdb_id = str(row["PDB_ID"])
		subs = literal_eval(row["Substrates"])
		cofs = literal_eval(row["Cofactors"])
		
		sub_chains = [x for x,y in subs]
		cof_chains = [x for x,y in cofs]
		common_chains = list(set(sub_chains).intersection(set(cof_chains)))
		
		if(len(common_chains)>0): #Cofactor and substrate present in same chain of structure
			sub_ids = []
			cof_ids = []
			
			chain_id = common_chains[0]  #Always first chain considered for binding site extraction
			sub_ids = [y for x,y in subs if x==chain_id]  #Accounts for presence of multiple substrates
			cof_ids = [y for x,y in cofs if x==chain_id]  #Accounts for presence of multiple cofactors
			
			#sub_pocket = extract_bs_resi(pdb_id, chain_id, sub_ids, pdbpath)
			#cof_pocket = extract_bs_resi(pdb_id, chain_id, cof_ids, pdbpath)
			#final_pocket = []
			#final_pocket.extend(sub_pocket)
			#final_pocket.extend(cof_pocket)
			
			#omit_resi = ["H_"+str(y) for y in sub_ids]  #Delete substrate copy in binding site, but retain cofactor copy
			
			sub_resi = ["H_"+str(y) for y in sub_ids]
			cof_resi = ["H_"+str(y) for y in cof_ids]
			retain_resi = []
			retain_resi.extend(sub_resi)
			retain_resi.extend(cof_resi)
			final_chain = extract_target_chain(pdb_id, chain_id, retain_resi, pdbpath)
			
			'''
			with warnings.catch_warnings():
				warnings.simplefilter('ignore')
				pdb = PDBParser().get_structure(pdb_id, pdbpath+pdb_id+".pdb")
				io = PDBIO()
				io.set_structure(pdb)
				io.save(outpath+pdb_id+"_pocket.pdb", ResidueSelect(final_pocket, omit_resi))
			'''
			
			with warnings.catch_warnings():
				warnings.simplefilter('ignore')
				pdb = PDBParser().get_structure(pdb_id, pdbpath+pdb_id+".pdb")
				io = PDBIO()
				io.set_structure(pdb)
				io.save(outpath+pdb_id+"_prep_v1.pdb", ResidueSelect(final_chain))
					
			print(pdb_id, chain_id, sub_resi, cof_resi)
		
		else:  #Cofactor and substrate are present in different chains of the same structure - NEEDS ALIGNMENT
			sub_ids = []
			cof_ids = []
			cof_chain_id = cof_chains[0]  #To be aligned to substrate chain
			sub_chain_id = sub_chains[0]  #Reference chain for alignment
			sub_ids = [y for x,y in subs if x==sub_chain_id]  #Accounts for presence of multiple substrates
			cof_ids = [y for x,y in cofs if x==cof_chain_id]  #Accounts for presence of multiple cofactors
			
			#---------------------------------------PYMOL ALIGNMENT AND SAVING COORDS--------------------------------------#
			cmd.load(pdbpath+pdb_id+".pdb", "mol1")
			#cmd.select("ref_chain", "chain "+str(sub_chain_id))
			#cmd.select("align_chain", "chain "+str(cof_chain_id))
			cmd.extract("ref_chain", "chain "+str(sub_chain_id))
			cmd.extract("align_chain", "chain "+str(cof_chain_id))
			cmd.align("align_chain", "ref_chain")
			
			sel_cmd = ["align_chain"]
			
			id_list = "resn "
			for j, idx in enumerate(cof_ids):
				if(j==0):
					id_list = id_list+str(idx)
				else:
					id_list = id_list+"+"+str(idx)
					
			sel_cmd.append(id_list)
			sel_txt = " and ".join(sel_cmd)		
			cmd.select(sel_txt)  #Extract cofactor atoms post alignment
			cmd.extract("cof_atoms", "sele")
			cmd.delete("align_chain")  #Delete aligned chain after extraction of cofactor atoms
			cmd.save(outpath+"tmp.pdb")
			cmd.reinitialize()
			#--------------------------------------------------------------------------------------------------------------#
			
			#sub_pocket = extract_bs_resi("tmp", sub_chain_id, sub_ids, outpath)
			#cof_pocket = extract_bs_resi("tmp", sub_chain_id, cof_ids, outpath)
			#final_pocket = []
			#final_pocket.extend(sub_pocket)
			#final_pocket.extend(cof_pocket)
			
			#omit_resi = ["H_"+str(y) for y in sub_ids]  #Delete substrate copy in binding site, but retain cofactor copy
			'''
			with warnings.catch_warnings():
				warnings.simplefilter('ignore')
				pdb = PDBParser().get_structure(pdb_id, outpath+"tmp.pdb")
				io = PDBIO()
				io.set_structure(pdb)
				io.save(outpath+pdb_id+"_pocket.pdb", ResidueSelect(final_pocket, omit_resi))
			'''
			
			sub_resi = ["H_"+str(y) for y in sub_ids]
			cof_resi = ["H_"+str(y) for y in cof_ids]
			retain_resi = []
			retain_resi.extend(sub_resi)
			retain_resi.extend(cof_resi)
			final_chain = extract_target_chain("tmp", sub_chain_id, retain_resi, outpath)
			
			with warnings.catch_warnings():
				warnings.simplefilter('ignore')
				pdb = PDBParser().get_structure(pdb_id, outpath+"tmp.pdb")
				io = PDBIO()
				io.set_structure(pdb)
				io.save(outpath+pdb_id+"_prep_v1.pdb", ResidueSelect(final_chain))
						
			os.remove(outpath+"tmp.pdb")
			print(pdb_id, sub_chain_id, sub_resi, cof_resi)
			
	elif(row["Structure_type"]=="Substrate only"):
		sub_ids = []
		pdb_id = str(row["PDB_ID"])
		subs = literal_eval(row["Substrates"])
		sub_chains = [x for x,y in subs]
		chain_id = sub_chains[0]  #Always first chain considered for binding site extraction
		sub_ids = [y for x,y in subs if x==chain_id]  #Accounts for presence of multiple substrates
		
		'''
		sub_pocket = extract_bs_resi(pdb_id, chain_id, sub_ids, pdbpath)
		final_pocket = []
		final_pocket.extend(sub_pocket)
		
		omit_resi = ["H_"+str(y) for y in sub_ids]  #Delete substrate copy in binding site
		with warnings.catch_warnings():
			warnings.simplefilter('ignore')
			pdb = PDBParser().get_structure(pdb_id, pdbpath+pdb_id+".pdb")
			io = PDBIO()
			io.set_structure(pdb)
			io.save(outpath+pdb_id+"_pocket.pdb", ResidueSelect(final_pocket, omit_resi))
		'''
		
		sub_resi = ["H_"+str(y) for y in sub_ids]
		retain_resi = []
		retain_resi.extend(sub_resi)
		final_chain = extract_target_chain(pdb_id, chain_id, retain_resi, pdbpath)
		
		with warnings.catch_warnings():
			warnings.simplefilter('ignore')
			pdb = PDBParser().get_structure(pdb_id, pdbpath+pdb_id+".pdb")
			io = PDBIO()
			io.set_structure(pdb)
			io.save(outpath+pdb_id+"_prep_v1.pdb", ResidueSelect(final_chain))
			
		print(pdb_id, chain_id, sub_resi, "[]")
		
	elif(row["Structure_type"]=="Apo structure"):
		sub_ids = []
		cof_ids = []
		
		pdb_id = str(row["PDB_ID"])
		chains = []
		with warnings.catch_warnings():
			warnings.simplefilter('ignore')
			pdb = PDBParser().get_structure(pdb_id, pdbpath+pdb_id+".pdb")
			io = PDBIO()
			io.set_structure(pdb)
			for model in pdb:
				for chain in model:
					chains.append(chain.id)
				break
		
		sub_apo_df = apo_df[apo_df["PDB_ID"]==pdb_id]
		holo_id = sub_apo_df.iloc[0]["Holo_mapping"]
		if(holo_id!="-"):
			sub_ids = []
			cof_ids = []
			
			holo_chain = str(holo_id.split("_")[1])
			pymol_cmd = sub_apo_df.iloc[0]["pymol"]
			sub_ids = sub_apo_df.iloc[0]["Substrate_ID"].split(",")
			#sub_pocket = extract_bs_resi(str(holo_id.split("_")[0]).lower(), holo_chain, sub_ids, apopath)
			#print(sub_pocket)
			
			cof_ids = []
			if(sub_apo_df.iloc[0]["Cofactor_ID"]!="-"):
				cof_ids = sub_apo_df.iloc[0]["Cofactor_ID"].split(",")
			
			#print(pdb_id, chains[0], holo_id, holo_chain, sub_ids, cof_ids)
			
			#---------------------------------------PYMOL ALIGNMENT AND SAVING COORDS--------------------------------------#
			cmd.load(pdbpath+pdb_id+".pdb", "mol1")
			cmd.select("mol1 and chain "+str(chains[0]))
			cmd.extract("align_chain", "sele")
			cmd.delete("mol1")
			cmd.load(apopath+str(holo_id.split("_")[0]).lower()+".pdb", "mol2")
			cmd.select("mol2 and chain "+str(holo_chain))
			cmd.extract("ref_chain", "sele")
			cmd.delete("mol2")
			
			if(pymol_cmd=="align"):
				cmd.align("align_chain", "ref_chain")
			else:
				cmd.super("align_chain", "ref_chain")
			
			sel_cmd = ["ref_chain"]
			id_list = "resn "
			for j, idx in enumerate(sub_ids):
				if(j==0):
					id_list = id_list+str(idx)
				else:
					id_list = id_list+"+"+str(idx)
					
			sel_cmd.append(id_list)	
			sel_txt = " and ".join(sel_cmd)
			
			if(len(cof_ids)>0):
				for idx in cof_ids:
					sel_txt = sel_txt+"+"+str(idx)
				
			cmd.select(sel_txt)  #Extract cofactor atoms post alignment
			cmd.extract("sub_cof_atoms", "sele")
			cmd.delete("ref_chain")  #Delete aligned chain after extraction of cofactor atoms
			cmd.save(outpath+"tmp.pdb")
			cmd.reinitialize()
			#--------------------------------------------------------------------------------------------------------------#
			
			'''
			sub_pocket = extract_bs_resi("tmp", holo_chain, sub_ids, outpath)
			final_pocket = []
			final_pocket.extend(sub_pocket)
			
			omit_resi = ["H_"+str(y) for y in sub_ids]  #Delete substrate copy in binding site, but retain cofactor copy
			with warnings.catch_warnings():
				warnings.simplefilter('ignore')
				pdb = PDBParser().get_structure(pdb_id, outpath+"tmp.pdb")
				io = PDBIO()
				io.set_structure(pdb)
				io.save(outpath+pdb_id+"_pocket.pdb", ResidueSelect(final_pocket, omit_resi))
			'''
			
			sub_resi = ["H_"+str(y) for y in sub_ids]
			if(len(cof_ids)>0):
				cof_resi = ["H_"+str(y) for y in cof_ids]
			else:
				cof_resi = "[]"
				
			retain_resi = []
			retain_resi.extend(sub_resi)
			if(len(cof_ids)>0):
				retain_resi.extend(cof_resi)
				
			final_chain = extract_target_chain("tmp", chains[0], retain_resi, outpath)
			
			with warnings.catch_warnings():
				warnings.simplefilter('ignore')
				pdb = PDBParser().get_structure(pdb_id, outpath+"tmp.pdb")
				io = PDBIO()
				io.set_structure(pdb)
				io.save(outpath+pdb_id+"_prep_v1.pdb", ResidueSelect(final_chain))
			
			os.remove(outpath+"tmp.pdb")
			print(pdb_id, chains[0], sub_resi, cof_resi)
		else:
			continue
		
	else:
		#print(row)
		sub_ids = []
		cof_ids = []
		
		pdb_id = str(row["PDB_ID"])
		chains = []
		with warnings.catch_warnings():
			warnings.simplefilter('ignore')
			pdb = PDBParser().get_structure(pdb_id, pdbpath+pdb_id+".pdb")
			io = PDBIO()
			io.set_structure(pdb)
			for model in pdb:
				for chain in model:
					chains.append(chain.id)
				break
		
		sub_cof_df = cof_df[cof_df["PDB_ID"]==pdb_id]
		holo_id = sub_cof_df.iloc[0]["Holo_mapping"]
		if(holo_id!="-"):
			sub_ids = []
			cof_ids = []
			holo_chain = str(holo_id.split("_")[1])
			pymol_cmd = sub_cof_df.iloc[0]["pymol"]
			sub_ids = sub_cof_df.iloc[0]["Substrate_ID"].split(",")
			cofs = literal_eval(row["Cofactors"])
			cof_ids = [y for x,y in cofs if x==holo_chain]  #TODO: Check if any cofactor IDs are same as substrate IDs
			
			#---------------------------------------PYMOL ALIGNMENT AND SAVING COORDS--------------------------------------#
			cmd.load(pdbpath+pdb_id+".pdb", "mol1")
			cmd.select("mol1 and chain "+str(chains[0]))
			cmd.extract("align_chain", "sele")
			cmd.delete("mol1")
			cmd.load(cofpath+str(holo_id.split("_")[0]).lower()+".pdb", "mol2")
			cmd.select("mol2 and chain "+str(holo_chain))
			cmd.extract("ref_chain", "sele")
			cmd.delete("mol2")
			
			if(pymol_cmd=="align"):
				cmd.align("align_chain", "ref_chain")
			else:
				cmd.super("align_chain", "ref_chain")
			
			sel_cmd = ["ref_chain"]
			id_list = "resn "
			for j, idx in enumerate(sub_ids):
				if(j==0):
					id_list = id_list+str(idx)
				else:
					id_list = id_list+"+"+str(idx)
					
			sel_cmd.append(id_list)
			sel_txt = " and ".join(sel_cmd)
		
			cmd.select(sel_txt)  #Extract cofactor atoms post alignment
			cmd.extract("sub_atoms", "sele")
			cmd.delete("ref_chain")  #Delete aligned chain after extraction of cofactor atoms
			cmd.save(outpath+"tmp.pdb")
			cmd.reinitialize()
			#--------------------------------------------------------------------------------------------------------------#
			
			'''
			sub_pocket = extract_bs_resi("tmp", holo_chain, sub_ids, outpath)
			final_pocket = []
			final_pocket.extend(sub_pocket)
			
			omit_resi = ["H_"+str(y) for y in sub_ids]  #Delete substrate copy in binding site, but retain cofactor copy
			with warnings.catch_warnings():
				warnings.simplefilter('ignore')
				pdb = PDBParser().get_structure(pdb_id, outpath+"tmp.pdb")
				io = PDBIO()
				io.set_structure(pdb)
				io.save(outpath+pdb_id+"_pocket.pdb", ResidueSelect(final_pocket, omit_resi))
			'''
			
			sub_resi = ["H_"+str(y) for y in sub_ids]
			if(len(cof_ids)>0):
				cof_resi = ["H_"+str(y) for y in cof_ids]
			else:
				cof_resi = "[]"
				
			retain_resi = []
			retain_resi.extend(sub_resi)
			if(len(cof_ids)>0):
				retain_resi.extend(cof_resi)
				
			final_chain = extract_target_chain("tmp", chains[0], retain_resi, outpath)
			
			with warnings.catch_warnings():
				warnings.simplefilter('ignore')
				pdb = PDBParser().get_structure(pdb_id, outpath+"tmp.pdb")
				io = PDBIO()
				io.set_structure(pdb)
				io.save(outpath+pdb_id+"_prep_v1.pdb", ResidueSelect(final_chain))
			
			os.remove(outpath+"tmp.pdb")
			
			print(pdb_id, chains[0], sub_resi, cof_resi)
		
	#break


































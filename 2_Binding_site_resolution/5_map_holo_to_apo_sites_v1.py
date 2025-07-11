#Program to map the binding site residues from holo structure to apo structure automatically

import os
import sys
import csv
import pandas as pd
import numpy as np
import Bio
from Bio import PDB, SeqIO, AlignIO
from Bio.PDB import PDBIO, PDBParser
import warnings

apo_data = sys.argv[1]
pdb_path = sys.argv[2]
msa_path = sys.argv[3]
site_path = sys.argv[4]
outfile = sys.argv[5]

apo_df = pd.read_csv(apo_data, sep="\t", header=0)
#print(apo_df.columns)  #['EC_number', 'UniProt_ID', 'Organism_name', 'Apo_structures', 'Holo_mapping']

out = open(outfile, "w")
print("PDB_ID\tBinding_site_residues\tSite_type", file=out)

three_to_one = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}
one_to_three = dict((c,i) for i,c in three_to_one.items())

for i, row in apo_df.iterrows():
	ec = row["EC_number"]
	apo_id = row["Apo_structures"].split(",")[0].split("_")[0]
	apo_chain = row["Apo_structures"].split(",")[0].split("_")[1]
	holo_id = row["Holo_mapping"].split("_")[0]
	holo_chain = row["Holo_mapping"].split("_")[1]

	holo_bs_resi = []
	site_type = ""
	sub_cof_df = pd.read_csv(site_path+"BRENDA_kcat_sub_and_cof_sites_v1.csv", sep="\t", header=0)
	sub_df = pd.read_csv(site_path+"BRENDA_kcat_substrate_only_sites_v1.csv", sep="\t", header=0)

	print(ec, apo_id, apo_chain, holo_id, holo_chain)
	test_df1 = sub_cof_df[sub_cof_df["PDB_ID"]==holo_id+"_"+holo_chain]
	if(len(test_df1.index) > 0):
		holo_bs_resi = test_df1.iloc[0]["Binding_site_residues"].split(",")
		site_type = "Substrate+Cofactor"
	else:
		test_df2 = sub_df[sub_df["PDB_ID"]==holo_id+"_"+holo_chain]
		holo_bs_resi = test_df2.iloc[0]["Binding_site_residues"].split(",")
		site_type = "Substrate_only"

	apo_fname = ""
	for fname in os.listdir(pdb_path):
		if(fname.startswith(str.lower(apo_id))):
			apo_fname = fname
			break

	apo_residict = {}
	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		pdb = PDBParser().get_structure(apo_id, pdb_path+apo_fname)
		io = PDBIO()
		io.set_structure(pdb)

		for model in pdb:
			for chain in model:
				if(chain.id==apo_chain):
					for resi in chain:
						try:
							apo_residict[resi.id[1]] = three_to_one[str(resi.resname)]
						except:
							continue
				else:
					continue

	holo_residict = {}
	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		pdb = PDBParser().get_structure(holo_id, pdb_path+str.lower(holo_id)+".pdb")
		io = PDBIO()
		io.set_structure(pdb)

		for model in pdb:
			for chain in model:
				if(chain.id==holo_chain):
					for resi in chain:
						try:
							holo_residict[resi.id[1]] = three_to_one[str(resi.resname)]
						except:
							continue
				else:
					continue

	align = AlignIO.read(msa_path+ec+"_msa.clu", "clustal")

	apo_record = ""
	holo_record = ""
	for rec in align._records:
		if(rec.id=="Apo_"+apo_id+"_"+str(apo_chain)):
			apo_record = rec
		elif(rec.id=="Gen_"+holo_id+"_"+str(holo_chain)):
			holo_record = rec

	bs_resi_final = []
	for resi in holo_bs_resi:
		bs_resi_final.append(int(resi[3:]))

	bs_resi_final = list(set(bs_resi_final))
	bs_resi_final.sort()
	apo_resi_final = []

	holo_resids = list(holo_residict.keys())
	apo_resids = list(apo_residict.keys())

	apo_index = 0
	holo_index = 0
	apo_counter = 0
	holo_counter = 0

	prev_holo = 0
	prev_apo = 0
	for hch, ach in zip(holo_record.seq, apo_record.seq):
		try:
			if(hch!="-"):
				holo_index = holo_resids[holo_counter]
				holo_counter = holo_counter + 1
			if(ach!="-"):
				apo_index = apo_resids[apo_counter]
				apo_counter = apo_counter + 1

			if(holo_index!=prev_holo and apo_index!=prev_apo and holo_index in bs_resi_final):
				apo_resi_final.append(one_to_three[apo_residict[apo_index]]+str(apo_index))

			prev_holo = holo_index
			prev_apo = apo_index
		except:
			continue

	print(str(row["Apo_structures"].split(",")[0])+"\t"+str(",".join(apo_resi_final))+"\t"+site_type, file=out)





























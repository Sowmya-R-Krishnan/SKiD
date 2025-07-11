#Program to get the mapping between PDB IDs and Chain IDs

import os
import sys
import csv
import pandas as pd
import numpy as np
import Bio
from Bio import PDB
from Bio.PDB import PDBIO, PDBParser
import warnings

inpath = sys.argv[1]
outfile = sys.argv[2]

out = open(outfile, "w")
print("PDB_ID\tChains", file=out)

for prot in os.listdir(inpath):
	chain_list = []
	pdb_id = prot.replace(".pdb","")

	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		pdb = PDBParser().get_structure(pdb_id, inpath+prot)
		io = PDBIO()
		io.set_structure(pdb)

		for model in pdb:
			for chain in model:
				chain_list.append(str(chain.id))
			break

	print(str(pdb_id)+"\t"+",".join(chain_list), file=out)

print("Done!")






























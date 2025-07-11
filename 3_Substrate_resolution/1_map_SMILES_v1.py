#Program to map IUPAC names to SMILES using OPSIN and pubchempy

import sys
import csv
import pandas as pd
import numpy as np
from py2opsin import py2opsin
import pubchempy as pcp

import warnings
warnings.filterwarnings('ignore')

infile = sys.argv[1]
mapped_outfile = sys.argv[2]
unmapped_outfile = sys.argv[3]

def iupac_2_smiles(iupac_name):
	smiles = py2opsin(chemical_name = iupac_name, output_format = "SMILES",)
	if len(smiles) > 0:
		#print('Using OPSIN')
		return smiles
	else:
		#print('Using PubChem')
		results = pcp.get_compounds(iupac_name, 'name')
		if len(results) > 0:
			return results[0].canonical_smiles
		else:
			return None

df = pd.read_csv(infile, sep="\t", header=0)
df.fillna("-", inplace=True)

out1 = open(mapped_outfile, "w") 
out2 = open(unmapped_outfile, "w")
print("SMILES\tSubstrate\tMol_ID", file=out1)
print("Substrate\tMol_ID", file=out2)

for i, row in df.iterrows():
	smiles = iupac_2_smiles(row["Substrate"])
	print(smiles)
	if(smiles):
		print(smiles+"\t"+row["Substrate"]+"\t"+str(row["Mol_ID"]), file=out1)
	else:
		print(row["Substrate"]+"\t"+str(row["Mol_ID"]), file=out2)






























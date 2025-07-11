#Program to use RDKit to obtain the InChIKey mapping for SMILES

import sys
import csv
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import inchi

infile = sys.argv[1]
outfile = sys.argv[2]

df = pd.read_csv(infile, sep="\t", header=0)
#print(df.columns)  #XXX: Atleast one column with SMILES header should be present

inchikeys = []
inchis = []

for i, row in df.iterrows():
	try:
		mol = Chem.MolFromSmiles(row["SMILES"])
		inch = inchi.MolToInchi(mol)
		inchikey = inchi.MolToInchiKey(mol)
		
		inchis.append(inch)
		inchikeys.append(inchikey)
	except:
		inchis.append("-")
		inchikeys.append("-")
		print(row["Ligand_name"])

df["InChI"] = inchis
df["InChIKey"] = inchikeys
df.to_csv(outfile, sep="\t", header=True, index=False)
print(df)

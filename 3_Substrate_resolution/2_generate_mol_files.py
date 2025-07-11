#Python script to write the mol files using RDKit ETKDG method of generating 3D coordinates

import os
import sys
import csv

from rdkit import Chem
from rdkit.Chem import AllChem

smilesfile = sys.argv[1]
savepath = sys.argv[2]

with open(smilesfile) as molecules:
	molreader = csv.reader(molecules, delimiter='\t', quotechar='"')
	for molecule in molreader:
		try:
			smiles = ''.join(molecule[0]).strip()
			mol_id = ''.join(molecule[1]).strip()
			mol_name = ''.join(molecule[2]).strip()
			rdkitmol = Chem.MolFromSmiles(smiles)
			mol2 = Chem.AddHs(rdkitmol)
			AllChem.EmbedMolecule(mol2, randomSeed=0xf00d)  #Random seed is for reproducibility
			AllChem.MMFFOptimizeMolecule(mol2)

			fname = savepath+"mol_"+str(mol_id)+".sdf"
			writer = Chem.SDWriter(fname)
			writer.write(mol2)
		except:
			fname = savepath+"mol_"+str(mol_id)+".sdf"
			print(smiles+"\t"+str(mol_id)+"\t"+mol_name)
			os.system("obabel -:\""+smiles+"\" -osdf -O "+fname+" --gen3d -h")
			continue

		


































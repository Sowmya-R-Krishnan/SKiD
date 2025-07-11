#Program to extract SIFTs mapping

import os
import sys
import csv
import pandas as pd
import numpy as np

inpath1 = sys.argv[1]
inpath2 = sys.argv[2]
outpath = sys.argv[3]

pdb_list = []
for fname in os.listdir(inpath1):
	pdb_id = str.lower(fname).replace(".pdb","")
	pdb_list.append(pdb_id)
	
for fname in os.listdir(inpath2):
	pdb_id = str.lower(fname).replace(".pdb","")
	pdb_list.append(pdb_id)

pdb_list = list(set(pdb_list))
print(len(pdb_list))
	
for i, fname in enumerate(pdb_list):
	if(str(fname)+".xml" not in os.listdir(outpath)):  #Download only if PDB SIFTS file is not already present in the output path
		os.system("wget \"https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/"+str(fname)+".xml.gz\"")
		os.system("mv "+str(fname)+".xml.gz "+outpath)
			
print("Done")































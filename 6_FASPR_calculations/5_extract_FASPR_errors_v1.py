#Program to extract errors in FASPR modelling of mutations

import os
import sys
import csv
import pandas as pd
import numpy as np

infile = sys.argv[1]
inpath = sys.argv[2]
outfile = sys.argv[3]

df = pd.read_csv(infile, sep="\t", header=0)
#print(df.columns)

out = open(outfile, "w")
print("affin.pdb_id\tmutation\taffin.chain", file=out)

for i, row in df.iterrows():
	mutations = row["mutation"].split("/")
	chain_id = row["affin.chain"]
	mutstr = ""
	for mut in mutations:
		mutstr = mutstr + chain_id + ":" + mut + "+"
	mutstr = mutstr[:-1]
	
	fname = row["affin.pdb_id"]+"_"+chain_id+"_"+mutstr+"_mut.pdb"
	if(os.path.exists(inpath+fname)):
		continue
	else:
		print(row["affin.pdb_id"]+"\t"+row["mutation"]+"\t"+chain_id, file=out)




























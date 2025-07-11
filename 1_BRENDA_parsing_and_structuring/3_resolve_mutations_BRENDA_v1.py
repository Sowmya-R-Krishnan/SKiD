#Program to extract mutation information and conditions from comments in BRENDA data

import os
import sys
import csv
import re
import pandas as pd
import numpy as np

inpath = sys.argv[1]
outpath = sys.argv[2]
multi_organism_outfile = sys.argv[3]

mout = open(multi_organism_outfile, "w")
print("EC_number\tSubstrate\tProtein_details\tPDB_IDs\tOrganism_name\tKinetic_value\tConditions\tReferences\tKinetic_type", file=mout)

for fname in os.listdir(inpath):
	data = pd.read_csv(inpath+fname, sep="\t", header=0)
	#print(data.columns)  #['Substrate', 'Protein_details', 'PDB_IDs', 'Organism_name', 'kcat_value', 'Conditions', 'References']

	cols_list = list(data.columns)[:-2]
	cols_list.extend(["Mutant", "Mutation", "pH", "Temperature", "References"])

	out = open(outpath+fname, "w")
	print("\t".join(cols_list), file=out)

	for i, row in data.iterrows():
		comments = row["Conditions"]
		org_matches = re.findall("(#{1}\d+#{1})", comments)
		
		try:
			if(len(org_matches)>1):
				row_list = [str(x) for x in row.values.flatten().tolist()]
				if("kcat_value" in cols_list):
					print(fname.replace(".csv", "").split("_")[0]+"\t"+"\t".join(row_list)+"\tkcat", file=mout)
				else:
					print(fname.replace(".csv", "").split("_")[0]+"\t"+"\t".join(row_list)+"\tKm", file=mout)

			else:
				contents = comments.split(" ")

				mutant = "no"
				mutation = "-----"
				ph_val = "-----"
				temp_val = "-----"

				mutations = []
				for j, val in enumerate(contents):
					if(re.search("[A-Z]{1}\d{1,}[A-Z]{1}", val)):
						mutations.append(val.replace(",", ""))
						mutant = "yes"
					elif(re.search("[pH]{2}", val) and contents[j+1]!="not"):  #To account for cases with (pH = XX), pH not specified
						ph_val = contents[j+1].replace(",", "").replace(")", "")
					elif(re.search("[°C]{2}", val) and contents[j+1]!="and" and contents[j+1]!="not"):  #To account for cases with temperature and pH not specified
						temp_val = val.replace(",", "").replace(")", "")
				
				if(len(mutations)==0):
					mutation = "-----"
					mutant = "no"
					if("kcat_value" in cols_list):
						print(row["Substrate"]+"\t"+row["Protein_details"]+"\t"+str(row["PDB_IDs"])+"\t"+row["Organism_name"]+"\t"+str(row["kcat_value"])+"\t"+mutant+"\t"+mutation+"\t"+str(ph_val)+"\t"+str(temp_val)+"\t"+str(row["References"]), file=out)
					else:
						print(row["Substrate"]+"\t"+row["Protein_details"]+"\t"+str(row["PDB_IDs"])+"\t"+row["Organism_name"]+"\t"+str(row["Km_value"])+"\t"+mutant+"\t"+mutation+"\t"+str(ph_val)+"\t"+str(temp_val)+"\t"+str(row["References"]), file=out)

				elif(len(mutations)==1):
					mutation = mutations[0]
					mutant = "yes"

					if("kcat_value" in cols_list):
						print(row["Substrate"]+"\t"+row["Protein_details"]+"\t"+str(row["PDB_IDs"])+"\t"+row["Organism_name"]+"\t"+str(row["kcat_value"])+"\t"+mutant+"\t"+mutation+"\t"+str(ph_val)+"\t"+str(temp_val)+"\t"+str(row["References"]), file=out)
					else:
						print(row["Substrate"]+"\t"+row["Protein_details"]+"\t"+str(row["PDB_IDs"])+"\t"+row["Organism_name"]+"\t"+str(row["Km_value"])+"\t"+mutant+"\t"+mutation+"\t"+str(ph_val)+"\t"+str(temp_val)+"\t"+str(row["References"]), file=out)

				elif(len(mutations)>1):
					for mut in mutations:
						mutation = mut
						mutant = "yes"

						if("kcat_value" in cols_list):
							print(row["Substrate"]+"\t"+row["Protein_details"]+"\t"+str(row["PDB_IDs"])+"\t"+row["Organism_name"]+"\t"+str(row["kcat_value"])+"\t"+mutant+"\t"+mutation+"\t"+str(ph_val)+"\t"+str(temp_val)+"\t"+str(row["References"]), file=out)
						else:
							print(row["Substrate"]+"\t"+row["Protein_details"]+"\t"+str(row["PDB_IDs"])+"\t"+row["Organism_name"]+"\t"+str(row["Km_value"])+"\t"+mutant+"\t"+mutation+"\t"+str(ph_val)+"\t"+str(temp_val)+"\t"+str(row["References"]), file=out)

		except:
			row_list = [str(x) for x in row.values.flatten().tolist()]
			if("kcat_value" in cols_list):
				print(fname.replace(".csv", "").split("_")[0]+"\t"+"\t".join(row_list)+"\tkcat", file=mout)
			else:
				print(fname.replace(".csv", "").split("_")[0]+"\t"+"\t".join(row_list)+"\tKm", file=mout)
			

	print(fname.replace(".csv", ""))






























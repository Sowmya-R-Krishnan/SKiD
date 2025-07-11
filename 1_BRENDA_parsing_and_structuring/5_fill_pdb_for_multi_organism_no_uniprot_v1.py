#Program to automate UniProt -> PDB mapping for multi-organism data with no UniProt IDs available

import sys
import csv
import os
import pandas as pd
import numpy as np
import json

data = sys.argv[1]
taxonomy_map = sys.argv[2]
outfile1 = sys.argv[3]  #Km file
outfile2 = sys.argv[4]  #kcat file
outfile3 = sys.argv[5]  #Manual check file

input_df = pd.read_csv(data, sep="\t", header=0)
tax_df = pd.read_csv(taxonomy_map, sep="\t", header=0)
#print(input_df.columns)  #['EC_number', 'Substrate', 'Organism_name', 'Kinetic_value', 'Conditions', 'References', 'Kinetic_type']
#print(tax_df.columns)  #['Organism', 'Taxonomy_ID']

out1 = open(outfile1, "w")
print("EC_number\tSubstrate\tProtein_details\tPDB_IDs\tOrganism_name\tKm_value\tConditions\tReferences", file=out1)

out2 = open(outfile2, "w")
print("EC_number\tSubstrate\tProtein_details\tPDB_IDs\tOrganism_name\tkcat_value\tConditions\tReferences", file=out2)

out3 = open(outfile3, "w")
print('EC_number\tSubstrate\tOrganism_name\tKinetic_value\tConditions\tReferences\tKinetic_type', file=out3)

#curl "https://rest.uniprot.org/uniprotkb/search?&query=taxonomy_id:833+AND+ec:3.2.1.73" > test.json
for i, row in input_df.iterrows():
	ec_num = row["EC_number"]
	orgn_name = row["Organism_name"]
	if(orgn_name not in ["Fungi", "Mammalia"] and not orgn_name.endswith("sp.")):
		tax_sub = tax_df[tax_df["Organism"]==orgn_name]
		tax_id = tax_sub["Taxonomy_ID"].item()

		if(tax_id!="-"):
			os.system("curl \"https://rest.uniprot.org/uniprotkb/search?&query=taxonomy_id:"+str(tax_id)+"+AND+ec:"+ec_num+"\" > test.json")
			with open("test.json", 'r') as f:
				uniprot_db = json.load(f)

			try:
				results = uniprot_db["results"][0]
				uniprot_id = results["primaryAccession"]

				pdb_list = []
				for crossref in results["uniProtKBCrossReferences"]:
					if(crossref["database"]=="PDB"):
						pdb_list.append(str(crossref["id"]))

				pdb_ids = ",".join(pdb_list)
				if(row["Kinetic_type"]=="Km"):
					print(str(ec_num)+"\t"+row["Substrate"]+"\t"+str(uniprot_id)+"\t"+pdb_ids+"\t"+orgn_name+"\t"+str(row["Kinetic_value"])+"\t"+row["Conditions"]+"\t"+str(row["References"]), file=out1)
				else:
					print(str(ec_num)+"\t"+row["Substrate"]+"\t"+str(uniprot_id)+"\t"+pdb_ids+"\t"+orgn_name+"\t"+str(row["Kinetic_value"])+"\t"+row["Conditions"]+"\t"+str(row["References"]), file=out2)

				os.remove("test.json")
			except:
				row_list = [str(x) for x in row.values.flatten().tolist()]
				print("\t".join(row_list), file=out3)
				os.remove("test.json")

		else:
			row_list = [str(x) for x in row.values.flatten().tolist()]
			print("\t".join(row_list), file=out3)

	else:
		row_list = [str(x) for x in row.values.flatten().tolist()]
		print("\t".join(row_list), file=out3)

































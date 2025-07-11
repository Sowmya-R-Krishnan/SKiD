#Program to extract the organism-specific kinetic and structure data from BRENDA JSON file

import sys
import csv
import json
import pandas as pd
import numpy as np

json_input = sys.argv[1]
outpath1 = sys.argv[2]  #Organisms
outpath2 = sys.argv[3]  #Proteins
outpath3 = sys.argv[4]  #References
outpath4 = sys.argv[5]  #kcat values
outpath5 = sys.argv[6]  #Km values

with open(json_input, 'r') as f:
	brenda_db = json.load(f)

print("BRENDA database loaded...")

ec_list = list(brenda_db['data'].keys())[1:]

for i, ec_num in enumerate(ec_list):
	try:
		organisms = brenda_db['data'][ec_num]['organisms']
		references = brenda_db['data'][ec_num]['references']
		proteins = brenda_db['data'][ec_num]['proteins']
		kcat = brenda_db['data'][ec_num]['turnover_number']
		km = brenda_db['data'][ec_num]['km_value']
		
		#Tabulate organism data
		with open(outpath1+ec_num+".csv", "w") as out1:
			print("Organism_ID\tOrganism_name", file=out1)
			for key, value in organisms:
				print(str(key)+"\t"+str(value['value']), file=out1)
		out1.close()

		#Tabulate protein data
		with open(outpath2+ec_num+".csv", "w") as out2:
			print("Protein_ID\tAccession\tSource\tComments", file=out2)
			for key, value in proteins:
				try:
					print(str(key)+"\t"+str(",".join(value['accessions']))+"\t"+str(value['source'])+"\t-----", file=out2)
				except:
					print(str(key)+"\t-----\t-----\t"+str(value['comment']), file=out2)
		out2.close()
		
		#Tabulate references
		with open(outpath3+ec_num+".csv", "w") as out3:
			print("ID\tTitle\tAuthors\tJournal\tYear\tPages\tVolume\tPMID", file=out3)
			for key, value in references:
				print(str(key)+"\t"+value['title']+"\t"+value['authors']+"\t"+value['journal']+"\t"+str(value['year'])+"\t"+str(value['pages'])+"\t"+str(value['vol'])+"\t"+str(value['pmid']), file=out3)
				
		out3.close()
		
		#Tabulate kcat values
		with open(outpath4+ec_num+".csv", "w") as out4:
			print("Substrate\tProtein_details\tOrganism_ID\tOrganism_name\tkcat_value\tConditions\tReferences", file=out4)
			for key, value in kcat:
				print(value['value']+"\t"+proteins[value['proteins'][0]]+"\t"+value['organisms'][0]+"\t"+organisms[value['organisms'][0]]['value']+"\t"+str(value['num_value'])+"\t"+value['comment']+"\t"+value['references'], file=out4)
				
		out4.close()
		
		#Tabulate Km values
		with open(outpath5+ec_num+".csv", "w") as out5:
			print("Substrate\tProtein_details\tOrganism_ID\tOrganism_name\tKm_value\tConditions\tReferences", file=out5)
			for key, value in km:
				print(value['value']+"\t"+proteins[value['proteins'][0]]+"\t"+value['organisms'][0]+"\t"+organisms[value['organisms'][0]]['value']+"\t"+str(value['num_value'])+"\t"+value['comment']+"\t"+value['references'], file=out5)
				
		out5.close()			
	except:
		print(ec_num+" failed.")
		continue

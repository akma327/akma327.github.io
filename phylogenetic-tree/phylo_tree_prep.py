# Author: Anthony Kai Kwang Ma
# Email: akma327@stanford.edu
# phylo_tree_prep.py

#!/share/PI/rondror/software/miniconda/bin/python

from __future__ import division 

import requests
import urllib2
import json

USAGE_STR = """

# Purpose 
# Generate Json File representing the phylogenetic tree for Class A GPCR Receptors

# USAGE:
# python phylo_tree_prep.py

"""


def Map_uniprot_to_gpcr_family(url):

	"""
		Retrieve GPCR data for given url requests
	"""

	uniprot_to_gpcr_family = {}
	response = urllib2.urlopen(url)
	gpcr_family_data = json.loads(response.read().decode('utf-8'))
	for entry in gpcr_family_data:
		key = entry['entry_name']
		if("human" in key):
			uniprot_to_gpcr_family[str(key.split("_human")[0])] = str(entry["family"])
	return uniprot_to_gpcr_family

def gen_phylo_tree_json(uniprot_to_gpcr_family):
	f = open("/scratch/PI/rondror/akma327/DynamicNetworks/src/analysis/phylogenetic-tree/phylo_data.txt", 'w')
	json_dict = {}
	for k in uniprot_to_gpcr_family:
		v = uniprot_to_gpcr_family[k]
		l1, l2, l3, l4 = v.strip().split("_")
		# json_dict["name"] = l1
		# json_dict["children"] = []
		f.write(l1 + "\t" + l2 + "\t" + l3 + "\t" + l4 + "\t" + k + "\n")
		# if(l1 not in json_dict):
		# 	json_dict[l1] = [{l2: [{l3: [{l4: [k]}]}]}]
		# else:
		# 	if(l2 not in json_dict):
		# 		json_dict[l1][l2] = [{l3: [{l4: [k]}]}]
		# 	else:
		# 		if(l3 not in json_dict):
		# 			json_dict[l1][l2][l3] = [{l4 : [k]}]
		# 		else:
		# 			if(l4 not in json_dict):
		# 				json_dict[l1][l2][l3][l4] = [k]
		# 			else:
		# 				json_dict[l1][l2][l3][l4].append(k)

	return json_dict

def get_full_gpcr_family():
	full_gpcr_family = {}
	for x in range(1,12):
		url = "http://gpcrdb.org:80/services/proteinfamily/proteins/001_" + "%03d" %(x,)
		uniprot_to_gpcr_family = Map_uniprot_to_gpcr_family(url)
		full_gpcr_family.update(uniprot_to_gpcr_family)
	return full_gpcr_family


def phylo_tree_driver():
	full_gpcr_family = get_full_gpcr_family()
	json_dict = gen_phylo_tree_json(full_gpcr_family)
	print json_dict

if __name__ == "__main__":
	phylo_tree_driver()



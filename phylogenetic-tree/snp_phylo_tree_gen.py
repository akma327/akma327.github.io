# Author: Anthony Kai Kwang Ma
# Email: akma327@stanford.edu
# snp_phylo_tree_gen.py

from __future__ import division 

import requests
import urllib2
import json
import sys
import math
import random

"""
# Purpose 
# Generate Json File representing the phylogenetic tree for Class A GPCR Receptors

# USAGE:
# python snp_phylo_tree_gen.py <DATA> <OUT>

# Arguments:
<DATA> List format data for class A GPCR phylogeny tree 
<OUT> Output json file for d3 js radial tree graph 

# Example

DATA="/scratch/PI/rondror/akma327/DynamicNetworks/data/phylogenetic-tree/phylo_data.txt"
OUT="/scratch/PI/rondror/akma327/DynamicNetworks/data/phylogenetic-tree/snp/flare-snp-gradient.json"
python snp_phylo_tree_gen.py $DATA $OUT

"""

red_gradient = ["#FFFFFF", "#FFFCFC", "#FFF9F9", "#FFF7F7", "#FFF4F4", "#FFF2F2", "#FFEFEF", "#FFECEC", "#FFEAEA", "#FFE7E7", "#FFE5E5", "#FFE2E2", "#FFE0E0", "#FFDDDD", "#FFDADA", "#FFD8D8", "#FFD5D5", "#FFD3D3", "#FFD0D0", "#FFCECE", "#FFCBCB", "#FFC8C8", "#FFC6C6", "#FFC3C3", "#FFC1C1", "#FFBEBE", "#FFBCBC", "#FFB9B9", "#FFB6B6", "#FFB4B4", "#FFB1B1", "#FFAFAF", "#FFACAC", "#FFAAAA", "#FFA7A7", "#FFA4A4", "#FFA2A2", "#FF9F9F", "#FF9D9D", "#FF9A9A", "#FF9797", "#FF9595", "#FF9292", "#FF9090", "#FF8D8D", "#FF8B8B", "#FF8888", "#FF8585", "#FF8383", "#FF8080", "#FF7E7E", "#FF7B7B", "#FF7979", "#FF7676", "#FF7373", "#FF7171", "#FF6E6E", "#FF6C6C", "#FF6969", "#FF6767", "#FF6464", "#FF6161", "#FF5F5F", "#FF5C5C", "#FF5A5A", "#FF5757", "#FF5555", "#FF5252", "#FF4F4F", "#FF4D4D", "#FF4A4A", "#FF4848", "#FF4545", "#FF4242", "#FF4040", "#FF3D3D", "#FF3B3B", "#FF3838", "#FF3636", "#FF3333", "#FF3030", "#FF2E2E", "#FF2B2B", "#FF2929", "#FF2626", "#FF2424", "#FF2121", "#FF1E1E", "#FF1C1C", "#FF1919", "#FF1717", "#FF1414", "#FF1212", "#FF0F0F", "#FF0C0C", "#FF0A0A", "#FF0707", "#FF0505", "#FF0202", "#FF0000"]
SNP_PATH="/scratch/PI/rondror/akma327/DynamicNetworks/data/phylogenetic-tree/snp/classA_missense_freq.txt"

def phylo_tree_list(DATA):
	phylo_list = []
	f = open(DATA, 'r')
	for line in f:
		leaves = line.strip().split("\t")
		phylo_list.append(leaves)
	return phylo_list


def leaf_key_exists(key, array):
	"""
		For an array of format [{"name": key1, "children": ___ },
								{"name": key2, "children": ___}]
		Check whether key is part of set {key1, key2} already
	"""
	for inner_dict in array:
		if(inner_dict["name"] == key): return True
	return False

def find_leaf_dict(key, array):
	"""
		For an array of format [{"name": key1, "children": ___ },
								{"name": key2, "children": ___}]
		Return the sub dictionary that contains desired key
	"""
	for inner_dict in array:
		if(inner_dict["name"] == key): return inner_dict
	return None

def uniprot_to_missense():
	"""	
		Load in snp count file for each human GPCR
	"""

	uniprot_snp_dict = {}
	f = open(SNP_PATH, 'r')
	for line in f:
		linfo = line.strip().split("\t")
		uniprot = linfo[0].split("_")[0]
		snp_missense = int(linfo[1])
		uniprot_snp_dict[uniprot] = snp_missense
	return uniprot_snp_dict, min(uniprot_snp_dict.values()), max(uniprot_snp_dict.values())

def getColor(name, uniprot_snp_dict, min_missense, max_missense):
	"""
		Scale from white to red based on number of missense in relation
		to min and max.
	"""
	m = 0
	key = name.upper()
	if (key in uniprot_snp_dict):
		m = uniprot_snp_dict[key]

	inc = float(max_missense - min_missense + 2)/100
	grad_index = int(math.ceil((m - min_missense)/inc))
	if(grad_index == 100): grad_index -= 1
	return red_gradient[grad_index]


def phylo_tree_json_arr(phylo_list, leaf_dict, uniprot_snp_dict, min_missense, max_missense):
	"""
		Generate the overall phylo_tree json_arr in levels manually
	"""
	json_arr = []
	for entry in phylo_list:
		l1_key = leaf_dict[entry[0]]

		if not leaf_key_exists(l1_key, json_arr):
			json_arr.append({"name": l1_key, "children": []})

	for entry in phylo_list:
		l1_key = leaf_dict[entry[0]]
		l2_key = leaf_dict[entry[0] + "_" + entry[1]]
		child_array = find_leaf_dict(l1_key, json_arr)["children"]
		if not leaf_key_exists(l2_key, child_array):
			child_array.append({"name": l2_key, "children": []})

	for entry in phylo_list:
		l1_key = leaf_dict[entry[0]]
		l2_key = leaf_dict[entry[0] + "_" + entry[1]]
		# l3_key = leaf_dict[entry[0] + "_" + entry[1] + "_" + entry[2]]
		l3_key = entry[2]
		child_array = find_leaf_dict(l2_key, find_leaf_dict(l1_key, json_arr)["children"])["children"]
		if not leaf_key_exists(l3_key, child_array):
			child_array.append({"name": l3_key, "children": []})

	for entry in phylo_list:
		# l1_key, l2_key, l3_key = entry[0], entry[1], entry[2]
		final_key = entry[4] # Fourth level
		l1_key = leaf_dict[entry[0]]
		l2_key = leaf_dict[entry[0] + "_" + entry[1]]
		# l3_key = leaf_dict[entry[0] + "_" + entry[1] + "_" + entry[2]]
		l3_key = entry[2]
		child_array = find_leaf_dict(l3_key, find_leaf_dict(l2_key, find_leaf_dict(l1_key, json_arr)["children"])["children"])["children"]
		if not leaf_key_exists(final_key, child_array):
			child_array.append({"name": final_key, "size": 2000, "colorid": getColor(final_key, uniprot_snp_dict, min_missense, max_missense)})

	return json_arr

def get_slug_array():
	url = "http://gpcrdb.org/services/proteinfamily/"
	response = urllib2.urlopen(url)
	slug_arr = json.loads(response.read().decode('utf-8'))
	leaf_dict = {}
	for slug_dict in slug_arr:
		slug, name = str(slug_dict['slug']).strip(), str(slug_dict['name']).strip()
		leaf_dict[slug] = name
	return leaf_dict


def phylo_tree_gen_driver(DATA, OUT):
	phylo_list = phylo_tree_list(DATA)
	leaf_dict_arr = get_slug_array()
	uniprot_snp_dict, min_missense, max_missense = uniprot_to_missense()

	json_arr = phylo_tree_json_arr(phylo_list, leaf_dict_arr, uniprot_snp_dict, min_missense, max_missense)
	json_dict = json_arr[0]
	with open(OUT, 'w') as fp:
		json.dump(json_dict, fp)


if __name__ == "__main__":
	(DATA, OUT) = (sys.argv[1], sys.argv[2])
	phylo_tree_gen_driver(DATA, OUT)

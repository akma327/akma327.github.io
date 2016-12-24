# Author: Anthony Kai Kwang Ma
# Email: akma327@stanford.edu

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
# python phylo_tree_gen.py <DATA> <OUT>

# Arguments:
<DATA> List format data for class A GPCR phylogeny tree 
<OUT> Output json file for d3 js radial tree graph 

# Example

DATA="/scratch/PI/rondror/akma327/DynamicNetworks/data/phylogenetic-tree/phylo_data.txt"
OUT="/scratch/PI/rondror/akma327/DynamicNetworks/data/phylogenetic-tree/snp/flare-snp-randint.json"
python phylo_tree_gen.py $DATA $OUT

"""

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

def genRandomColor():
	letters = '0123456789ABCDEF'
	color = "0x"
	for i in range(6):
		color += letters[int(math.floor(random.random()*16))]
	return int(color, 16)
	return color


# def phylo_tree_json_arr(phylo_list, leaf_dict_arr):
# 	"""
# 		Generate the overall phylo_tree json_arr in levels manually
# 	"""
# 	json_arr = []
# 	for entry in phylo_list:
# 		# l1_key = entry[0] #First level
# 		l1_key = leaf_dict_arr[0][entry[0]]
# 		if not leaf_key_exists(l1_key, json_arr):
# 			json_arr.append({"name": l1_key, "children": []})

# 	for entry in phylo_list:
# 		# l1_key = entry[0]
# 		# l2_key = entry[1] # Second level
# 		l1_key = leaf_dict_arr[0][entry[0]]
# 		l2_key = leaf_dict_arr[1][entry[1]]
# 		child_array = find_leaf_dict(l1_key, json_arr)["children"]
# 		if not leaf_key_exists(l2_key, child_array):
# 			child_array.append({"name": l2_key, "children": []})

# 	for entry in phylo_list:
# 		# l1_key, l2_key = entry[0], entry[1]
# 		# l3_key = entry[2] # Third level
# 		l1_key = leaf_dict_arr[0][entry[0]]
# 		l2_key = leaf_dict_arr[1][entry[1]]
# 		l3_key = entry[2]
# 		child_array = find_leaf_dict(l2_key, find_leaf_dict(l1_key, json_arr)["children"])["children"]
# 		if not leaf_key_exists(l3_key, child_array):
# 			child_array.append({"name": l3_key, "children": []})

# 	for entry in phylo_list:
# 		# l1_key, l2_key, l3_key = entry[0], entry[1], entry[2]
# 		final_key = entry[4] # Fourth level
# 		l1_key = leaf_dict_arr[0][entry[0]]
# 		l2_key = leaf_dict_arr[1][entry[1]]
# 		l3_key = entry[2]
# 		child_array = find_leaf_dict(l3_key, find_leaf_dict(l2_key, find_leaf_dict(l1_key, json_arr)["children"])["children"])["children"]
# 		if not leaf_key_exists(final_key, child_array):
# 			child_array.append({"name": final_key, "size": 2000})

# 	return json_arr

# def get_slug_array():
# 	url = "http://gpcrdb.org/services/proteinfamily/"
# 	response = urllib2.urlopen(url)
# 	slug_arr = json.loads(response.read().decode('utf-8'))
# 	leaf_dict_arr = [{}, {}, {}, {}]
# 	for slug_dict in slug_arr:
# 		slug, name = str(slug_dict['slug']).strip(), str(slug_dict['name']).strip()
# 		slug_parts = slug.split("_")
# 		leaf_dict_arr[len(slug_parts) - 1][slug_parts[-1]] = name
# 	return leaf_dict_arr


def phylo_tree_json_arr(phylo_list, leaf_dict):
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
			child_array.append({"name": final_key, "size": 2000, "colorid": "#F0BCEE"})

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
	json_arr = phylo_tree_json_arr(phylo_list, leaf_dict_arr)
	json_dict = json_arr[0]
	with open(OUT, 'w') as fp:
		json.dump(json_dict, fp)


if __name__ == "__main__":
	(DATA, OUT) = (sys.argv[1], sys.argv[2])
	phylo_tree_gen_driver(DATA, OUT)

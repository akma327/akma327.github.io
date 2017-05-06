# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# uniprot_to_crystal_pdb.py

import json 
import sys 
import collections

USAGE_STR = """
# Purpose 
# Generate mapping from subset of uniprot to solved crystal pdb structures 

# Usage 
# python uniprot_to_crystal_pdb.py <OUTPUT_FILE>

# Arguments
# <OUTPUT_FILE> Output file path 

# Example 
OUTPUT_FILE="/scratch/PI/rondror/akma327/GeneticVariations/data/phylogenetic-tree/uniprot_to_pdb/uniprot_to_crystal_pdb.json"
cd /scratch/PI/rondror/akma327/GeneticVariations/src/phylogenetic-tree
python uniprot_to_crystal_pdb.py $OUTPUT_FILE

"""

K_MIN_ARG = 2

class_A_pdb_list = "/scratch/PI/rondror/akma327/GeneticVariations/data/phylogenetic-tree/uniprot_to_pdb/pdb_list.txt"
gpcr_pdb_list = "/scratch/PI/rondror/akma327/DynamicNetworks/data/crystal-analysis/simulation-analysis/gpcrdb-freq-config/GPCR_PDB_List.txt"


def calc_uniprot_to_snp():
	"""
		For every uniprot generate a dictionary from gpcrdb position to 
		number of missense mutations at that position. 
	"""

	f = open(missense_mutations, 'r')
	uniprot_to_snp = {}
	counter = 0
	for line in f:
		# if(counter > 500): break
		linfo = line.strip().split("\t")
		uniprot, resid, gpcrdb, allele_freq = linfo[0].lower(), int(linfo[1]), linfo[2], float(linfo[10])
		# print(uniprot, resid, gpcrdb, allele_freq)
		if(uniprot not in uniprot_to_snp):
			uniprot_to_snp[uniprot] = {resid: 1}
		else:
			if(resid not in uniprot_to_snp[uniprot]):
				uniprot_to_snp[uniprot][resid] = 1
			else:
				uniprot_to_snp[uniprot][resid] += 1
		counter +=1

	# for key in uniprot_to_snp:
	# 	uniprot_to_snp[key] = collections.Counter(uniprot_to_snp[key])
	return uniprot_to_snp


def gen_uniprot_to_pdb_mapping(OUTPUT_FILE):
	"""
		Generate a mapping between uniprot and list of pdbs 
	"""

	f1 = open(class_A_pdb_list, 'r')
	f2 = open(gpcr_pdb_list, 'r')

	pdb_list = []
	for line in f1:
		pdb_list.append(line.strip())

	pdb_to_uniprot = {}
	for line in f2:
		if(line == "\n"): continue
		linfo = line.strip().split("\t")
		uniprot, pdb = linfo[0], linfo[2]
		pdb_to_uniprot[pdb] = uniprot

	uniprot_to_pdbs = {}
	for pdb in pdb_list:
		if(pdb in pdb_to_uniprot):
			uniprot = pdb_to_uniprot[pdb].lower()
			if(uniprot not in uniprot_to_pdbs):
				uniprot_to_pdbs[uniprot] = [pdb.upper()]
			else:
				uniprot_to_pdbs[uniprot].append(pdb.upper())

	json_output = []
	for uniprot in uniprot_to_pdbs:
		pdbs = ",".join(uniprot_to_pdbs[uniprot])
		json_output.append({"uniprot": uniprot, "pdbs" : pdbs})

	with open(OUTPUT_FILE, 'w') as f: 
		json.dump(json_output, f)


def driver(OUTPUT_FILE):
	gen_uniprot_to_pdb_mapping(OUTPUT_FILE)


if __name__ == "__main__":
	if(len(sys.argv) < K_MIN_ARG):
		print(USAGE_STR)
		exit(0)
	(OUTPUT_FILE) = (sys.argv[1])
	driver(OUTPUT_FILE)

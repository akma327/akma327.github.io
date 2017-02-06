# Author: Anthony Kai Kwang Ma
# Email: akma327@stanford.edu
# uniprot_to_snp.py

import json 
import sys 
import collections

USAGE_STR = """
# Purpose 
# Generate dictionary from uniprot to snp data per residue

# Usage 
# python uniprot_to_snp.py <OUTPUT_FILE>

# Arguments
# <OUTPUT_FILE> Output file path

# Example 
OUTPUT_FILE="/scratch/PI/rondror/akma327/GeneticVariations/data/phylogenetic-tree/uniprot_to_pdb/uniprot_to_snp.json"
cd /scratch/PI/rondror/akma327/GeneticVariations/src/phylogenetic-tree
python uniprot_to_snp.py $OUTPUT_FILE

"""

K_MIN_ARG = 2

missense_mutations = "/scratch/PI/rondror/akma327/GeneticVariations/data/gpcrdb-mutations/classA_human_missense_mutations.txt"


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

	return uniprot_to_snp


def driver(OUTPUT_FILE):
	uniprot_to_snp = calc_uniprot_to_snp()
	with open(OUTPUT_FILE, 'w') as f:
		json.dump(uniprot_to_snp, f)


if __name__ == "__main__":
	if(len(sys.argv) < K_MIN_ARG):
		print(USAGE_STR)
		exit(0)
	(OUTPUT_FILE) = (sys.argv[1])
	driver(OUTPUT_FILE)

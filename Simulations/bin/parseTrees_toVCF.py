import msprime, allel, tskit, pyslim, gzip ## Special packages
import argparse, glob, random ## Standard packages
import pandas as pd ## Pandas 
import numpy as np ## Numpy


## A script to take the tree sequence files from either MSprime or SLiM and calculate summary statistics

def read_slim(slim_file, keepMuts = False, mutation_rate = 2.5e-6, sampleSize = 10, vcf = "none"):
## Use this function to sample individuals from the SLiM trees
	
## Get a list of the individuals from each sub-pop and their corresponding nodes
## There's no need for recapitation as I ran the sims until coalescence

## Read in the slim trees
	orig_ts = pyslim.load(slim_file)

## Make a dictionary containing the data on all individuals alive at time 0 (i.e. the end of the simulation)
	individuals_to_test = [ j for j in orig_ts.individuals_alive_at(0)]

## Make a dict containing all the alive individuals' nodes in the tree
	pop_dict = {}
	for i in orig_ts.individuals():
		if i.id not in individuals_to_test: continue
		try:
			pop_dict[i.population].append(i.nodes)
		except KeyError:
			pop_dict[i.population] = [i.nodes]


## For each of the 10 populations, grab a sample of 10 individuals (20 nodes) from the node dicts
	nodes_to_keep_list = []
	for p in range(1):
		nodes_to_keep_list.append( random.sample( pop_dict[p+1], sampleSize ) )
	nodes_to_keep = np.sort( np.concatenate(nodes_to_keep_list, axis = None) )

## Simplify the riginal sequence with the 
	simplified_tree_seq = orig_ts.simplify( nodes_to_keep )

	if keepMuts:
		return simplified_tree_seq
	elif keepMuts == False:
		ts_added = pyslim.SlimTreeSequence(msprime.mutate(simplified_tree_seq,
		 rate=mutation_rate, 
		 keep=False))
		return ts_added

	
def main():
	parser = argparse.ArgumentParser(description="")

	parser.add_argument("--tree", 
		required = True,
		dest = "tree",
		type = str, 
		help = "The  tree sequences")

	args = parser.parse_args()
		
## Do this for SLiM trees		
	for tree in glob.glob(args.tree + "/*trees"):
		nameString = tree.split("/")[-1]
		print(nameString)
		if len(nameString.split(".")) == 3:
			gen = 30000
			rep = int(nameString.split(".")[1])
		else:
			gen = int(nameString.split(".")[1][3:])
			rep = int(nameString.split(".")[2][3:])
		if gen > 33000 or rep > 10: continue	
		print(tree)	
		mutated_tree = read_slim( tree , keepMuts = False)
		
		vcf_name = tree.split("/")[-1]
		with gzip.open(vcf_name+".vcf.gz", "wt") as f:
		    mutated_tree.write_vcf(f)
		    


main()

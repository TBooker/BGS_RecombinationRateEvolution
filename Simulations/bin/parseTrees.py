import msprime, allel, tskit, pyslim ## Special packages
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
	parser.add_argument("--rec", 
		required = True,
		dest = "rec",
		type = str, 
		help = "The recombiantion rate map")

	parser.add_argument("--output", 
		required = True,
		dest = "output",
		type = str,
		help = "Give the name of the output file")
 		

	args = parser.parse_args()
	
	output = []

	r_regime_1 = [1e-8, 2e-8, 3e-8, 4e-7, 5e-7, 6e-7, 7e-7, 8e-7, 9e-7, 10e-7]
	r_regime_2 = r_regime_1[::-1]
	
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
		if gen > 33000: continue	
		mutated_tree = read_slim( tree , keepMuts = False)	

		pos_vector = [int(v.position)+1 for v in mutated_tree.variants()]
		print(str(len(pos_vector)) + "\n")
		msprime_genotype_matrix = mutated_tree.genotype_matrix()
# Convert msprime's haplotype matrix into genotypes by randomly merging chromosomes
		haplotype_array = allel.HaplotypeArray( msprime_genotype_matrix )
		genotype_array = haplotype_array.to_genotypes(ploidy=2)
		ac = genotype_array.count_alleles()
		windowWidth = 10000

		pi, windows, n_bases, counts = allel.windowed_diversity(pos_vector, ac, size=windowWidth, start=1, stop=10000000)

		recom_length = 1000000
		
		
		window_locations = windows.sum(axis = 1)/2

		if gen >30000:
			r = r_regime_2
		else:
			r = r_regime_1
	
		
		window_r = sum(([r[i]]*int(recom_length/windowWidth) for i in range(len(r))),[]) 

		dat =  pd.DataFrame( [pi, window_r, window_locations] ).transpose()
		dat["gen"] = gen
		dat["rep"] = rep
		output.append(dat)
	pd.concat(output).to_csv(args.output, index = False)


main()

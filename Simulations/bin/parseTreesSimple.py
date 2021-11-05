import msprime, allel, tskit, pyslim ## Special packages
import argparse, glob, random ## Standard packages
import pandas as pd ## Pandas
import numpy as np ## Numpy
from multiprocessing import Pool


## A script to take the tree sequence files from either MSprime or SLiM and calculate summary statistics

def read_slim(slim_file, keepMuts = False, mutation_rate = 0.5e-6, sampleSize = 50):
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

def analyse_tree(argsFromPool):

	treeSeq = argsFromPool[0]

	r_regime_1 = argsFromPool[1][0]

	print(treeSeq)

	nameString = treeSeq.split("/")[-1]

	gen = int( nameString.split('.')[1][3:] )

	if len( nameString.split('.') ) == 7:
		R = float(nameString.split('.')[4] + "." + nameString.split('.')[5] )
	elif len( nameString.split('.') ) == 6:
		R = float(nameString.split('.')[4] )
	rep = nameString.split('.')[2][3:]

	r_regime_2 = [r_regime_1[0]*R]

	mutated_tree = read_slim( treeSeq , keepMuts = False)

	pos_vector = [int(v.position)+1 for v in mutated_tree.variants()]

	msprime_genotype_matrix = mutated_tree.genotype_matrix()

# Convert msprime's haplotype matrix into genotypes by randomly merging chromosomes
	haplotype_array = allel.HaplotypeArray( msprime_genotype_matrix )
	genotype_array = haplotype_array.to_genotypes(ploidy=2)
	ac = genotype_array.count_alleles()
	windowWidth = 5000

	pi, windows, n_bases, counts = allel.windowed_diversity(pos_vector, ac, size=windowWidth, start=1, stop=25000)

	window_locations = windows.sum(axis = 1)/2

	if gen >100001:
		r = r_regime_2
	else:
		r = r_regime_1

	dat =  pd.DataFrame( [pi, window_locations] ).transpose()
	dat["gen"] = gen
	dat["rep"] = rep
	dat["R"] = R

	return(dat)

def main():
	parser = argparse.ArgumentParser(description="")

	parser.add_argument("--tree",
		required = True,
		dest = "tree",
		type = str,
		help = "The  tree sequences directory")

	parser.add_argument("--output",
		required = True,
		dest = "output",
		type = str,
		help = "Give the name of the output file")

	parser.add_argument("--nProc",
		required = False,
		dest = "nProc",
		type = int,
		help = "Give the number of threads for this program [5]")
	args = parser.parse_args()

	output = []
	rRegime1 = [5e-7]

	arg_set = [ [tree, [rRegime1]] for tree in glob.glob(args.tree + "/*trees") ]

## Do this for SLiM trees

	with Pool(args.nProc) as p:
		output= p.map(analyse_tree, arg_set)
	pd.concat(output).to_csv(args.output, index = False)


main()

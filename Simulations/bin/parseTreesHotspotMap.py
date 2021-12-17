import msprime, allel, tskit, pyslim ## Special packages
import argparse, glob, random ## Standard packages
import pandas as pd ## Pandas
import numpy as np ## Numpy


## A script to take the tree sequence files from either MSprime or SLiM and calculate summary statistics

def read_slim(slim_file, keepMuts = False, mutation_rate = 2.5e-7, sampleSize = 50, refGen = 0):
## Use this function to sample individuals from the SLiM trees

## Get a list of the individuals from each sub-pop and their corresponding nodes
## There's no need for recapitation as I ran the sims until coalescence

## Read in the slim trees
	orig_ts = pyslim.load(slim_file)
## Make a dictionary containing the data on all individuals alive at time 0 (i.e. the end of the simulation)
	individuals_to_test = [ j for j in orig_ts.individuals_alive_at(refGen)]

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


def parseRecRates( rate_file , hotspot_r = 6e-6):

	hotspots_raw = [ i.strip() for i in open(rate_file).readlines() ]

	num_hotspots_r1 = int(hotspots_raw.index("regime_2") /2)

	num_hotspots_r2 = int((len(hotspots_raw) - ((num_hotspots_r1*2) + 2))/2)

	regime_1_locs = [1]+[int(k)+1 for  k in  hotspots_raw[1:int(num_hotspots_r1*2)]]
	regime_2_locs = [1]+[int(k)+1 for  k in  hotspots_raw[int(num_hotspots_r1*2+1):]]
#	print(regime_2_locs)

	r_1 = [hotspot_r/100] + [hotspot_r/100, hotspot_r]*num_hotspots_r1 + [hotspot_r/100]
	r_2 = [hotspot_r/100] + [hotspot_r/100, hotspot_r]*num_hotspots_r2 + [hotspot_r/100]


	rec_DF_r1 = pd.DataFrame([regime_1_locs,r_1]).transpose()
	rec_DF_r1.columns = ['breaks', 'rates']
	rec_DF_r2 = pd.DataFrame([regime_2_locs,r_2]).transpose()
	rec_DF_r2.columns = ['breaks', 'rates']

	return {"r1":rec_DF_r1, "r2":rec_DF_r2}


def get_r_in_windows( windows, rates ):

	outRates = []

	for w in windows:

		if w[0] < 1 or w[1] > 1e7:
			outRates.append(np.nan)

			continue
		mini_DF = rates[(rates.breaks < w[1])&(rates.breaks >= w[0])].copy()
		mini_DF.reset_index()
		if len (mini_DF) == 0 or len (mini_DF) == 1 and mini_DF.index[0] ==0:
			midpoint = (w[0] + w[1]) /2
			temp = rates.copy()
			temp[ "dist_to_mid" ] = temp[ "breaks" ] - midpoint
#			print(temp[temp["dist_to_mid"] > 0] )
			outRates.append( temp[temp["dist_to_mid"] > 0].iloc[0]['rates'])
#			print(temp[temp["dist_to_mid"] > 0].iloc[0]['rates'])
		else:

			if mini_DF.index.min() ==0:
				new_mini_DF = rates[mini_DF.index.min():mini_DF.index.max()+2].copy()
			else:
				new_mini_DF = rates[mini_DF.index.min()-1:mini_DF.index.max()+2].copy()

			new_mini_DF.reset_index()
			new_mini_DF.iloc[-1]["breaks"] = w[1]
			new_mini_DF.iloc[0]["breaks"] = w[0]-1

			new_mini_DF["diff"] = new_mini_DF.breaks.diff()
#		mini_DF.iloc[0]["diff"] = list(mini_DF["breaks"])[0] - w[0]
			new_mini_DF.loc[0, 'diff'] = list(new_mini_DF["breaks"])[0] - w[0]
			new_mini_DF["CumRec"] = (new_mini_DF["diff"] * new_mini_DF["rates"])
#			print(new_mini_DF)
			outRates.append(new_mini_DF.CumRec.sum()/( 1 + w[1] - w[0]))
#			print(new_mini_DF.CumRec.sum()/( 1 + w[1] - w[0]))

	return outRates

def main():
	parser = argparse.ArgumentParser(description="")

	parser.add_argument("--tree",
		required = True,
		dest = "tree",
		type = str,
		help = "The  tree sequences")
#	parser.add_argument("--rec",
#		required = True,
#		dest = "rec",
#		type = str,
#		help = "The recombiantion rate map")
	parser.add_argument("--win",
		required = True,
		dest = "win",
		type = int,
		help = "The width of the analysis window")
	parser.add_argument("--recWin",
		required = False,
		dest = "recWin",
		type = int,
		help = "The width of the recombination rate windows. Default is to use args.win",
		default = 1)
	parser.add_argument("--output",
		required = True,
		dest = "output",
		type = str,
		help = "Give the name of the output file")


	args = parser.parse_args()

	output = []


#	r_regime_1 = [1e-8, 2e-8, 3e-8, 4e-7, 5e-7, 6e-7, 7e-7, 8e-7, 9e-7, 10e-7]
#	r_regime_2 = r_regime_1[::-1]
	print(args.win, args.recWin)
	if args.recWin == 1:
		recWin = args.win
	elif args.recWin < args.win:
		return
	else:
		recWin = args.recWin

## Do this for SLiM trees
	for tree in glob.glob(args.tree + "/*trees"):


		nameString = tree.split("/")[-1]
		print( nameString.split(".") )
		if len(nameString.split(".")) == 3:
			gen = 80000
			rep = int(nameString.split(".")[1])
			recRates = parseRecRates(args.tree + "/replicate_" + str(rep) + ".recombinationHotspots.txt")
			mutated_tree = read_slim( tree , keepMuts = False, refGen = 1)
		else:
			gen = int(nameString.split(".")[1][3:])
			rep = int(nameString.split(".")[2][3:])
			recRates = parseRecRates(args.tree + "/replicate_" + str(rep) + ".recombinationHotspots.txt")
			mutated_tree = read_slim( tree , keepMuts = False, refGen = 0)
		if gen > 160001: continue
#		print(tree)
#		if rep != 2 and gen != 135000: continue


		pos_vector = [int(v.position)+1 for v in mutated_tree.variants()]

		msprime_genotype_matrix = mutated_tree.genotype_matrix()

# Convert msprime's haplotype matrix into genotypes by randomly merging chromosomes
		haplotype_array = allel.HaplotypeArray( msprime_genotype_matrix )
		genotype_array = haplotype_array.to_genotypes(ploidy=2)
		ac = genotype_array.count_alleles()
		windowWidth = args.win

		pi, windows, n_bases, counts = allel.windowed_diversity(pos_vector, ac, size=windowWidth, start=1, stop=10000000)

#		temp_pi, rec_windows, temp_n_bases, counts = allel.windowed_diversity(pos_vector, ac, size=recWin, start=1, stop=10000000)

#		recom_length = 1000000

		window_locations = windows.sum(axis = 1)/2

		if gen >80000:
			regime = "r2"
		else:
			regime = "r1"

		if args.recWin == 1 or args.recWin/args.win==1.:
			window_r = get_r_in_windows(windows, recRates[regime])
			print(window_r)
		elif args.recWin > args.win:
			## do something cool
			print(windows)
			w_to_add = int((args.recWin - args.win)/2)
			rec_windows = np.array([[dw[0]-w_to_add,dw[1]+w_to_add]  for dw in windows])
			window_r = get_r_in_windows(rec_windows, recRates[regime])
		else:
			print("sorry, does not compute")
			return


		dat =  pd.DataFrame( {"pi":pi, "window_r":window_r, "window_locations":window_locations} )
		dat["gen"] = gen
		dat["rep"] = rep

		output.append(dat)

	pd.concat(output).to_csv(args.output, index = False)


main()

import argparse, glob
import pandas as pd
## Read in recMaps inferred by PyRho and make a nice CSV for plotting

	
def main():
	parser = argparse.ArgumentParser(description="")

	parser.add_argument("--rmaps", 
		required = True,
		dest = "rmaps",
		type = str, 
		help = "The dir with recombination maps")

	parser.add_argument("--output", 
		required = True,
		dest = "output",
		type = str,
		help = "Give the name of the output file")
 		

	args = parser.parse_args()
	
	output = []

## Do this for all rmaps	
	for rmap in glob.glob(args.rmaps + "/*rmap"):
	
		nameString = rmap.split("/")[-1]
		print(nameString)
		if len(nameString.split(".")) == 3:
			gen = 30000
			rep = int(nameString.split(".")[1])
		else:
			gen = int(nameString.split(".")[1][3:])
			rep = int(nameString.split(".")[2][3:])
		R = float(nameString.split(".")[4])
		if R == 0:
			R = 0.1
		
		
		csv = pd.read_csv(rmap, header = None, sep = "\t")
		span  = csv[1].max() - csv[0].min()
		cum_rate = ( ( csv[1] - csv[0] ) * csv[2] ).sum()  / span

		output.append( [rep, R, gen, cum_rate] )		
	pd.DataFrame(output, columns = ['rep','R','generation','r_rate']).to_csv(args.output, index = False)
main()
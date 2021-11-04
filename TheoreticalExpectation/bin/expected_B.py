## A script to calculate the expected B for a neutral site located at a specified position on a given chromosome map
import scipy.stats
import pandas as pd
from scipy.integrate import quad
import numpy as np
import sys
## The simulated chromosome is 15,000bp long with a 5,000bp functional region in the centre
## All mutations that occur in the functional region have an equal chance of mutating to a deleterious allele
## Deleterious mutational effects are drawn from a gamma distribution with mean -0.1 and shape 0.1
## Recombination occurs at a constant rate of 2.5e-6
## At generation X it changes from 2.5e-6 to (2.5e-6)*R
## The population size is N=1000 diploids
def B_gamma_integrand( s, u_x, r_xv, g_shape, g_mean):
	t = s*0.5
	f_t = scipy.stats.gamma.pdf(t, g_shape, loc = 0, scale = g_mean*0.5)
	return (u_x * f_t)/(t*(1 + ((1-t)*r_xv)/t)**2 ) 

def B_fixed_s( s, u_x, r_xv):
	t = s*0.5
	return (u_x * 1.)/(t*(1 + ((1-t)*r_xv)/t)**2 ) 
	
def B_all_sites( focalSite, 
				funcSites,
				recRate,
				mutRate,
				gamma_shape,
				gamma_mean):
	B_summer = 0
	for j in funcSites:
		recDist = abs(focalSite - j) * recRate	
		if recDist == 0:
			recDist = 1.*recRate
#		B_summer +=  quad(B_gamma_integrand, 0.005/100, 1, args=( mutRate, recDist, gamma_shape, gamma_mean))[0]
		B_summer +=  B_fixed_s( 0.05/1000, mutRate, recDist)

	return np.exp(-1.*B_summer)

print(sys.argv[1])		
functionalSites = [x for x in range(10000, 15000) ]

output = []

for R in [0.1, 1.0, 10.0]:
#for R in [1.0]:
	for i in range(1, 25000, 1000):
#	for i in [1,  7500, 12500, 17500, 22500]:
#	print( "Physical dist to selected target:", i -12500 )
#	print( "Genetic dist to selected target:", (i -12500)*2.5e-6 )
		print(i, R)
		output.append( [i, B_all_sites(i, functionalSites, R*2.5e-6/1000, 2.5e-6/1000, 0.1, 0.1), R])
		
pd.DataFrame(output, columns = ["pos","B","R"]).to_csv(sys.argv[1], index = False)

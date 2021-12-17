## Script to run the analysis of the recombination rate map evolution simulations

#I ran the simulations on a server that uses the SLURM scheduler, the following script runs an array job, for 30 simulations:


# The data from the server end up zipped (actually tszipped), so I can download them more quickly

# the first thing to do is to unzip the data:
parallel "tsunzip {}" ::: $(ls broadScale_3x/*tsz)
parallel "tsunzip {}" ::: $(ls hotspots/*tsz)


parallel -j2 "python bin/parseTreesRecRateMap.py --tree broadScale_3x/ --win {1} --recWin {2} --output broadScale_3x_analysis/broadScale_3x_w{1}_r{2}.csv " ::: 10000 100000 1000000 ::: 10000 100000 1000000

parallel -j4 "python bin/parseTreesHotspotMap.py --tree hotspots/ --win {1} --recWin {2} --output hotspots_map_analysis/hotspot_map_w{1}_r{2}.csv" ::: 10000 100000 1000000 ::: 10000 100000 1000000

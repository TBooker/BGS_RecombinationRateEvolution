rm(list = ls())

## Make a plot showing the recombination rate landscapes before and after the change for the broad-scale and fine scale simulations

## Hotspots 

# Read in the file of starts and stops from a single simulation rep
hotspots <- read.table("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations/hotspotMap/replicate_1.recombinationHotspots.txt")

# Extract the vector corresponding to the first set of hotspots
initial_hotspots <- c(1, as.numeric( hotspots$V1[2:(nrow(hotspots)/2)] ) )
# Add an integer corresponding to the position immediately before each hotspot
initial_hotspots
# Extract the vector corresponding to the second set of hotspots
final_hotspots <- c(1, as.numeric( hotspots$V1[(2 + nrow(hotspots)/2):nrow(hotspots)] ) )

# Make a vector for the recombination rates within each hotspot
rates <- c( rep(rep(c(2.08e-7, 2.08e-5),each = 2), times = length(initial_hotspots)/4), 2.08e-7, 2.08e-7)

hotspot_maps <- data.frame( initial_hotspots = initial_hotspots, final_hotspots = final_hotspots, rates = rates )

library(ggplot2)
library(reshape2)
melted_hotspots <- melt(hotspot_maps, id = "rates")

ggplot(data = melted_hotspots, aes(x = value, y = rates))+
  geom_line()+
  facet_wrap(~variable)+
  scale_y_log10(limits = c(1e-8, 1e-3))

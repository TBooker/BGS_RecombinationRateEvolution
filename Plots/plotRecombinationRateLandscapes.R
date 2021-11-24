rm(list = ls())

## Make a plot showing the recombination rate landscapes before and after the change for the broad-scale and fine scale simulations
setwd("~/UBC/Mouse_BGS_RecombinationRateEvolution/Plots")

## Hotspots 
# Read in the file of starts and stops from a single simulation rep
hotspots <- read.table("../Simulations/hotspotMap/hs_10000_randomNumber/replicate_7.recombinationHotspots.txt")
# Get the index for where the regime changes in the file
initial_break <- (1:nrow(hotspots))[hotspots$V1 == "regime_2"]
# Extract the vector corresponding to the first set of hotspots
initial_hotspots <- c(1, as.numeric(as.character(hotspots$V1[2:(initial_break-1)] ) ) )
# Add an integer corresponding to the position immediately before each hotspot

initial_hotspots <- sort( c( initial_hotspots, 
                             initial_hotspots[seq(2,length(initial_hotspots)-1,2)]-1, 
                             initial_hotspots[seq(3,length(initial_hotspots)-1,2)]+1) )
# Make a vector for the recombination rates within each hotspot
intitial_rates <- c( rep(rep(c( 4.16e-8,  4.16e-6),each = 2), times = length(initial_hotspots)/4), 4.16e-8, 4.16e-8)


# Extract the vector corresponding to the second set of hotspots
final_hotspots <- c(1, as.numeric(as.character(hotspots$V1[(initial_break+1):nrow(hotspots)] )))
final_hotspots <- sort( c( final_hotspots, 
                           final_hotspots[seq(2,length(final_hotspots)-1,2)]-1, 
                           final_hotspots[seq(3,length(final_hotspots)-1,2)]+1) )
# Make a vector for the recombination rates within each hotspot
final_rates <- c( rep(rep(c( 4.16e-8,  4.16e-6),each = 2), times = length(final_hotspots)/4), 4.16e-8, 4.16e-8)


initial_hotspot_map <- data.frame( hotspots = initial_hotspots,  rates = intitial_rates )
final_hotspot_map <- data.frame( hotspots = final_hotspots,  rates = final_rates )

initial_hotspot_map$source <- "Before"
final_hotspot_map$source <- "After"

library(ggplot2)
library(reshape2)

hotspot_maps <- rbind(initial_hotspot_map, final_hotspot_map)
melted_hotspots <- melt(hotspot_maps, id = c("rates","source"))

melted_hotspots$source <- factor(melted_hotspots$source, 
                                 levels = c("Before","After"))

hotspot_plot <- ggplot(data = melted_hotspots, aes(x = value/1e6, y = rates*4*5000))+
  geom_line()+
  facet_wrap(source~., ncol =1 )+
  scale_x_continuous("Position in Genome (Mbp)",
                     breaks = 0:10)+
  scale_y_log10(expression("Recombination Rate (4"*italic(N[e]*"r")*")"),
                limit = c(0.0005,0.5))+
  theme_bw()+
  theme(
    strip.background.x = element_blank(),
    strip.text.x = element_text( size = 10)
  )

# There, plotted the hotspot map nicely in ggplot2

# Now let's do the same for the broad scale map - this is a lot less involved
recRates <- c(2.5e-7, 3e-7, 4e-7, 4.5e-7, 5e-7, 5e-7, 5.5e-7, 6e-7, 7e-7, 7.5e-7)
recRates <- rep(recRates, each = 2)
starts <- seq(1,10e6,1e6)
ends <- starts+999999
broadScale_before <- data.frame(position = sort(c(starts, ends)),
                                rates = recRates,
                                source = "Before")
broadScale_after <- data.frame(position = sort(c(starts, ends)),
                                rates = rev(recRates),
                                source = "After")

broadScale_maps <- rbind(broadScale_before, broadScale_after)

broadScale_maps$source <- factor(broadScale_maps$source, 
                                 levels = c("Before","After"))


broadScale_plot <- ggplot(data = broadScale_maps, aes(x = position/1e6, y = rates*4*5000))+
  geom_line()+
  facet_wrap(source~., ncol = 1)+
  scale_x_continuous("Position in Genome (Mbp)",
                     breaks = 0:10)+
  scale_y_log10(expression("Recombination Rate (4"*italic(N[e]*"r")*")"),
                limit = c(0.0005,0.5))+
  theme_bw()+
  theme(
    strip.background.x = element_blank(),
    strip.text.x = element_text( size = 10)
  )



library(ggpubr)

map_schematic <- ggarrange( broadScale_plot,
           hotspot_plot,
           ncol = 2,
           labels = "AUTO", 
           align = "h")
pdf("recombinationMapDiagram.pdf", width = 10, height = 6)
print(map_schematic)
dev.off()  

recRates <- c(2.5e-7, 3e-7, 4e-7, 4.5e-7, 5e-7, 5e-7, 5.5e-7, 6e-7, 7e-7, 7.5e-7) *(4.164/5)
sum(recRates*(10^6))

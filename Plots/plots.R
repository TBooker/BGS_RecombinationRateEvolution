rm (list  = ls())

library(ggplot2)
options(scipen=999)

getWindowDataFrame <- function(data_raw, win, label){
  cor_df<-list()
  for (g in seq(30000,33000,500)){
    temp_data <- data_raw[data_raw$gen == g,]
    gen_cor <- with(temp_data, cor.test( pi, r, method = "kendall"))
    gen_out = c(win, g, label, mean(temp_data$pi),
                sd(temp_data$pi), 
                gen_cor$estimate, 
                gen_cor$p.value)
    temp_df <- rbind.data.frame(gen_out)
    names(temp_df)<-c("Window",
                      "Generation",
                      "Label",
                      "Ave.Diversity",
                      "SD.Diversity",
                      "Rho.est",
                      "Rho.p")
    cor_df <- append(cor_df, list(temp_df))
  }
  
  cor_df_merged <- do.call(rbind, cor_df)
  return(cor_df_merged)
}  

r0.3_10000 <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations/broadScale_range0.3.10000.csv")
names (r0.3_10000 ) <- c("pi","POS","gen","rep","r")
win_0.3_10000 <- getWindowDataFrame(r0.3_10000, "10,000 bp", "0.3")

r0.3_100000 <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations/broadScale_range0.3.100000.csv")
names (r0.3_10000 ) <- c("pi","POS","gen","rep","r")
win_0.3_100000 <- getWindowDataFrame(r0.3_10000, "100,000 bp", "0.3")

r0.3_1000000 <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations/broadScale_range0.3.1000000.csv")
names (r0.3_1000000 ) <- c("pi","POS","gen","rep","r")
win_0.3_1000000 <- getWindowDataFrame(r0.3_1000000, "1,000,000 bp", "0.3")
# 
# r0.111_10000 <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations/broadScale_range0.111.10000.csv")
# names (r0.111_10000 ) <- c("pi","POS","gen","rep","r")
# win_0.111_10000 <- getWindowDataFrame(r0.111_10000, "10 Kbp", "0.111")
# 
# r0.111_10000 <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations/broadScale_range0.111.100000.csv")
# names (r0.111_10000 ) <- c("pi","POS","gen","rep","r")
# win_0.111_100000 <- getWindowDataFrame(r0.111_10000, "100 Kbp", "0.111")
# 
# r0.111_1000000 <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations/broadScale_range0.111.1000000.csv")
# names (r0.111_1000000 ) <- c("pi","POS","gen","rep","r")
# win_0.111_1000000 <- getWindowDataFrame(r0.111_1000000, "1 Mbp", "0.111")


# allDat_0.111 <- rbind(win_0.111_10000,win_0.111_100000,win_0.111_1000000)

allDat_0.3 <- rbind(win_0.3_10000,win_0.3_100000, win_0.3_1000000)
allDat_0.3$Window <- factor(allDat_0.3$Window,
                            levels = c("10,000 bp",
                                       "100,000 bp",
                                       "1,000,000 bp"))
corPlot_broadScale <- ggplot(data = allDat_0.3, aes(x = (as.numeric(as.character(Generation)) - 30000)/2000, y =   as.numeric(as.character(Rho.est)), col =   as.numeric(as.character(Rho.p)) < 0.0001))+
  geom_point()+
  facet_grid(Window~.)+
  scale_y_continuous(expression("Spearman's "*italic(rho) *" (Recombination v. Diversity)"))+
  scale_x_continuous("Time Since Recombination Rate Changed\n(2N Generations)")+
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey")+
  scale_color_manual(expression(p<0.0001), values = c("red","black"))+
  theme_bw()

pdf("~/UBC/Mouse_BGS_RecombinationRateEvolution/Plots/simpleMap03_pi_r_correlationOverTime.pdf", height = 8, width = 5)
print(corPlot_0.3)
dev.off()
# 
# corPlot_0.111 <- ggplot(data = allDat_0.111, aes(x = (as.numeric(as.character(Generation)) - 30000)/2000, y =   as.numeric(as.character(Rho.est)), col =   as.numeric(as.character(Rho.p)) < 0.0001))+
#   geom_point()+
#   facet_grid(Window~.)+
#   scale_y_continuous(expression("Spearman's "*italic(rho) *" (Recombination v. Diversity)"))+
#   scale_x_continuous("Time Since Recombination Rate Changed (2N Generations)")+
#   geom_hline(aes(yintercept = 0), lty = 2, col = "grey")+
#   scale_color_manual(expression(p<0.0001), values = c("red","black"))+
#   theme_bw()
# 
# pdf("~/UBC/Mouse_BGS_RecombinationRateEvolution/Plots/simpleMap0111_pi_r_correlationOverTime.pdf", height = 8, width = 5)
# print(corPlot_0.111)
# dev.off()
# 



corPlot_0.3 <- ggplot(data = allDat_0.3, aes(x = (as.numeric(as.character(Generation)) - 30000)/2000, y =   as.numeric(as.character(Rho.est)), col =   as.numeric(as.character(Rho.p)) < 0.0001))+
  geom_point()+
  facet_grid(Window~.)+
  scale_y_continuous(expression("Spearman's "*italic(rho) *" (Recombination v. Diversity)"))+
  scale_x_continuous("Time Since Recombination Rate Changed (2N Generations)")+
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey")+
  scale_color_manual(expression(p<0.0001), values = c("red","black"))+
  theme_bw()

pdf("~/UBC/Mouse_BGS_RecombinationRateEvolution/Plots/for_talk_simpleMap03_pi_r_correlationOverTime.pdf", height = 8, width = 5)
print(corPlot_0.3)
dev.off()

corPlot_0.3 <- ggplot(data = win_0.3_100000, aes(x = (as.numeric(as.character(Generation)) - 30000)/2000, y =   as.numeric(as.character(Rho.est)), col =   as.numeric(as.character(Rho.p)) < 0.0001))+
  geom_point()+
  facet_grid(~Window)+
  scale_y_continuous(expression("Spearman's "*italic(rho) *" (Recombination v. Diversity)"))+
  scale_x_continuous("Time Since Recombination Rate Changed (2N Generations)")+
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey")+
  scale_color_manual(expression(p<0.0001), values = c("red","black"))+
  theme_bw()

getWindowDataFrame <- function(data_raw, window){
  cor_df<-list()
  for (g in c(15000, 15500, 17000, 16000, 17500, 18000, 16500)){
    temp_data <- data_raw[data_raw$gen == g,]
    gen_cor <- with(temp_data, cor.test( pi, r, method = "spearman"))
    gen_out = c(window, g, "Gough Island", mean(temp_data$pi),
                sd(temp_data$pi), 
                gen_cor$estimate, 
                gen_cor$p.value)
    temp_df <- rbind.data.frame(gen_out)
    names(temp_df)<-c("Window",
                      "Generation",
                      "Population",
                      "Ave.Diversity",
                      "SD.Diversity",
                      "Rho.est",
                      "Rho.p")
    cor_df <- append(cor_df, list(temp_df))
  }
  
  cor_df_merged <- do.call(rbind, cor_df)
  return(cor_df_merged)
}  

raw_10000 <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations//hotspotMap/hs_10000_randomNumber/win10Kbp.csv")
names (raw_10000 ) <- c("pi","r","POS","gen","rep")
win_10000 <- getWindowDataFrame(raw_10000, "10,000 bp")
# 
# raw_50000 <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations//hotspotMap/hs_10000_randomNumber/win50Kbp.csv")
# names (raw_50000 ) <- c("pi","r","POS","gen","rep")
# win_50000 <- getWindowDataFrame(raw_50000, "50 Kbp")

raw_100000 <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations//hotspotMap/hs_10000_randomNumber/win100Kbp.csv")
names (raw_100000 ) <- c("pi","r","POS","gen","rep")
win_100000 <- getWindowDataFrame(raw_100000, "100,000 bp")
# 
# raw_500000 <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations/hotspotMap/hs_10000_randomNumber/win500Kbp.csv")
# names (raw_500000 ) <- c("pi","r","POS","gen","rep")
# win_500000 <- getWindowDataFrame(raw_500000, "500 Kbp")

raw_1000000 <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations/hotspotMap/hs_10000_randomNumber/win1Mbp.csv")
names (raw_1000000 ) <- c("pi","r","POS","gen","rep")
win_1000000 <- getWindowDataFrame(raw_1000000, "1,000,000 bp")


allDat <- rbind(win_10000, win_100000, win_1000000)
allDat$Window <- factor(allDat$Window,
                        levels = c("10,000 bp",
                                   "100,000 bp",
                                   "1,000,000 bp"))

## For sup figure
corPlot_hotSpots <- ggplot(data = allDat, aes(x = (as.numeric(as.character(Generation)) - 15000)/2000, y =   as.numeric(as.character(Rho.est)), col =   as.numeric(as.character(Rho.p)) < 0.0001))+
  geom_point()+
  facet_grid(Window~.)+
  scale_y_continuous(expression("Spearman's "*italic(rho) *" (Recombination v. Diversity)"), limits = c(-0.3,0.5))+
  scale_x_continuous("Time Since Recombination Rate Changed\n(2N Generations)")+
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey")+
  scale_color_manual(expression(p<0.0001), values = c("red","black"))+
  theme_bw()

### For talk
corPlot <- ggplot(data = allDat, aes(x = (as.numeric(as.character(Generation)) - 15000)/2000, y =   as.numeric(as.character(Rho.est)), col =   as.numeric(as.character(Rho.p)) < 0.0001))+
  geom_point()+
  facet_grid(~Window)+
  scale_y_continuous(expression("Spearman's "*italic(rho) *" (Recombination v. Diversity)"), limits = c(-0.3,0.5))+
  scale_x_continuous("Time Since Recombination Rate Changed\n(2Ne Generations)")+
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey")+
  scale_color_manual(expression(p<0.0001), values = c("red","black"))+
  theme_bw()+
  theme(
    strip.background.x =  element_blank(),
    strip.text.x =  element_blank()
  )


pdf("~/UBC/Mouse_BGS_RecombinationRateEvolution/Plots/pi_r_correlationOverTime.pdf", height = 4, width = 4)
print(corPlot)
dev.off()


library(ggpubr)


corPlot_simple <- ggplot(data = win_0.3_100000, aes(x = (as.numeric(as.character(Generation)) - 30000)/2000, y =   as.numeric(as.character(Rho.est)), col =   as.numeric(as.character(Rho.p)) < 0.0001))+
  geom_point()+
  facet_grid(~Window)+
  scale_y_continuous(expression("Spearman's "*italic(rho) *" (Recombination v. Diversity)"), limits = c(-0.3,0.5))+
  scale_x_continuous("Time Since Recombination Rate Changed\n(2N Generations)")+
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey")+
  scale_color_manual(expression(p<0.0001), values = c("red","black"))+
  theme_bw()+
  theme(
    strip.background.x =  element_blank(),
    strip.text.x =  element_blank()
  )


corPlot_hotspot <- ggplot(data =win_100000, aes(x = (as.numeric(as.character(Generation)) - 30000)/2000, y =   as.numeric(as.character(Rho.est)), col =   as.numeric(as.character(Rho.p)) < 0.0001))+
  geom_point()+
  facet_grid(~Window)+
  scale_y_continuous(expression("Spearman's "*italic(rho) *" (Recombination v. Diversity)"), limits = c(-0.3,0.5))+
  scale_x_continuous("Time Since Recombination Rate Changed\n(2N Generations)")+
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey")+
  scale_color_manual(expression(p<0.0001), values = c("red","black"))+
  theme_bw()+
  theme(
    strip.background.x =  element_blank(),
    strip.text.x =  element_blank()
  )


Figure2 <- ggarrange(corPlot_simple,
          corPlot_hotspot,
          nrow=1,
          ncol=2,
          labels = "AUTO",
         common.legend = T,
         legend = "right")
 

pdf("~/UBC/Mouse_BGS_RecombinationRateEvolution/Plots/pi_r_correlationOverTime_bothMaps.pdf", height = 4, width = 8)
print(Figure2)
dev.off()




combined_suppFig5 <- ggarrange(corPlot_broadScale, corPlot_hotSpots,
          nrow= 1,
          ncol = 2,
          common.legend = T,
          labels = "AUTO",
          legend = "right")


pdf("~/UBC/Mouse_BGS_RecombinationRateEvolution/Plots/pi_r_correlationOverTime_bothMaps_allWindows.pdf", height = 8, width = 8)
print(combined_suppFig5)
dev.off()

rm(list= ls())
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
win_0.3_10000 <- getWindowDataFrame(r0.3_10000, "10 Kbp", "0.3")

r0.3_10000 <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations/broadScale_range0.3.100000.csv")
names (r0.3_10000 ) <- c("pi","POS","gen","rep","r")
win_0.3_100000 <- getWindowDataFrame(r0.3_10000, "100 Kbp", "0.3")

r0.3_1000000 <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations/broadScale_range0.3.1000000.csv")
names (r0.3_1000000 ) <- c("pi","POS","gen","rep","r")
win_0.3_1000000 <- getWindowDataFrame(r0.3_1000000, "1 Mbp", "0.3")

r0.111_10000 <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations/broadScale_range0.111.10000.csv")
names (r0.111_10000 ) <- c("pi","POS","gen","rep","r")
win_0.111_10000 <- getWindowDataFrame(r0.111_10000, "10 Kbp", "0.111")

r0.111_10000 <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations/broadScale_range0.111.100000.csv")
names (r0.111_10000 ) <- c("pi","POS","gen","rep","r")
win_0.111_100000 <- getWindowDataFrame(r0.111_10000, "100 Kbp", "0.111")

r0.111_1000000 <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations/broadScale_range0.111.1000000.csv")
names (r0.111_1000000 ) <- c("pi","POS","gen","rep","r")
win_0.111_1000000 <- getWindowDataFrame(r0.111_1000000, "1 Mbp", "0.111")


allDat_0.111 <- rbind(win_0.111_10000,win_0.111_100000,win_0.111_1000000)
allDat_0.3 <- rbind(win_0.3_10000,win_0.3_100000,win_0.3_1000000)

corPlot_0.3 <- ggplot(data = allDat_0.3, aes(x = (as.numeric(as.character(Generation)) - 30000)/2000, y =   as.numeric(as.character(Rho.est)), col =   as.numeric(as.character(Rho.p)) < 0.0001))+
  geom_point()+
  facet_grid(~Window)+
  scale_y_continuous(expression("Spearman's "*italic(rho) *" (Recombination v. Diversity)"))+
  scale_x_continuous("Time Since Recombination Rate Changed (2N Generations)")+
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey")+
  scale_color_manual(expression(p<0.0001), values = c("red","black"))+
  theme_bw()

pdf("~/UBC/Mouse_BGS_RecombinationRateEvolution/Plots/simpleMap03_pi_r_correlationOverTime.pdf", height = 8, width = 5)
print(corPlot_0.3)
dev.off()

corPlot_0.111 <- ggplot(data = allDat_0.111, aes(x = (as.numeric(as.character(Generation)) - 30000)/2000, y =   as.numeric(as.character(Rho.est)), col =   as.numeric(as.character(Rho.p)) < 0.0001))+
  geom_point()+
  facet_grid(Window~.)+
  scale_y_continuous(expression("Spearman's "*italic(rho) *" (Recombination v. Diversity)"))+
  scale_x_continuous("Time Since Recombination Rate Changed (2N Generations)")+
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey")+
  scale_color_manual(expression(p<0.0001), values = c("red","black"))+
  theme_bw()

pdf("~/UBC/Mouse_BGS_RecombinationRateEvolution/Plots/simpleMap0111_pi_r_correlationOverTime.pdf", height = 8, width = 5)
print(corPlot_0.111)
dev.off()



allDat_0.3 <- rbind(win_0.3_10000,win_0.3_100000)

corPlot_0.3 <- ggplot(data = allDat_0.3, aes(x = (as.numeric(as.character(Generation)) - 30000)/2000, y =   as.numeric(as.character(Rho.est)), col =   as.numeric(as.character(Rho.p)) < 0.0001))+
  geom_point()+
  facet_grid(~Window)+
  scale_y_continuous(expression("Spearman's "*italic(rho) *" (Recombination v. Diversity)"))+
  scale_x_continuous("Time Since Recombination Rate Changed (2N Generations)")+
  geom_hline(aes(yintercept = 0), lty = 2, col = "grey")+
  scale_color_manual(expression(p<0.0001), values = c("red","black"))+
  theme_bw()

pdf("~/UBC/Mouse_BGS_RecombinationRateEvolution/Plots/for_talk_simpleMap03_pi_r_correlationOverTime.pdf", height = 8, width = 5)
print(corPlot_0.3)
dev.off()



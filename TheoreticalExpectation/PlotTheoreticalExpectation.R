rm(list = ls())
library(ggplot2)
## Using the formula for the coalescence time for a step change in B
## (1 + (R1 - 1)*exp(-1*T_0)) from Johri et al 2020

R1 = 0.9/0.92
R2 = 0.9/0.99
R3 = 1/R1
R4 = 1/R2
T_0 = 1:600/100

B1 = 0.67
B2 = 0.81
plot( T_0, B1*(1 + (B2/B1 - 1)*exp(-1*T_0)), type = "l", xlab = "Time since the recombination rate changed (2Ne)", ylab = "Net Coalesence Time", ylim = c(0.8,1.1))
lines( T_0, (1 + (R2 - 1)*exp(-1*T_0)), type = "l")
lines( T_0, (1 + (R3 - 1)*exp(-1*T_0)), type = "l")
lines( T_0, (1 + (R4 - 1)*exp(-1*T_0)), type = "l")
abline(h = 1, lty = 2)

##### Gamma DFE

BGS_DF_part1 <- data.frame( B = 0.67*(1 + (0.81/0.67 - 1)*exp(-1*T_0)) , gen = T_0, R = 0.1, pos = 12500.5)
BGS_DF_part2 <- data.frame( B = 0.81*(1 + (0.81/0.81 - 1)*exp(-1*T_0)) , gen = T_0, R = 1, pos = 12500.5)
BGS_DF_part3 <- data.frame( B = 0.95*(1 + (0.81/0.95 - 1)*exp(-1*T_0)) , gen = T_0, R = 10, pos = 12500.5)

BGS_DF_part4 <- data.frame( B = 0.74*(1 + (0.93/0.74 - 1)*exp(-1*T_0)) , gen = T_0, R = 0.1, pos = 17500.5)
BGS_DF_part5 <- data.frame( B = 0.92*(1 + (0.93/0.93 - 1)*exp(-1*T_0)) , gen = T_0, R = 1, pos = 17500.5)
BGS_DF_part6 <- data.frame( B = 0.98*(1 + (0.93/0.98 - 1)*exp(-1*T_0)) , gen = T_0, R = 10, pos = 17500.5)

BGS_DF_part7 <- data.frame( B = 0.8*(1 + (0.95/0.8 - 1)*exp(-1*T_0)) , gen = T_0, R = 0.1, pos = 22500.5)
BGS_DF_part8 <- data.frame( B = 0.97*(1 + (0.97/0.97 - 1)*exp(-1*T_0)) , gen = T_0, R = 1, pos = 22500.5)
BGS_DF_part9 <- data.frame( B = 1.0*(1 + (0.96/1.0 - 1)*exp(-1*T_0)) , gen = T_0, R = 10, pos = 22500.5)

sims <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations/summary_simple_map_constant_s0.01.csv")
names(sims) <- c("pi","pos","gen","rep","R")
unique(sims$pos)
library(ggplot2)

BGS <- rbind(BGS_DF_part1,BGS_DF_part2,BGS_DF_part3,
      BGS_DF_part4,BGS_DF_part5,BGS_DF_part6,
      BGS_DF_part7,BGS_DF_part8,BGS_DF_part9)

BGS$pos_label <- abs(BGS$pos-12500.5)
BGS$pos_label <- factor( BGS$pos_label, c(0,5000,10000), labels = c("Within the functional element", "5Kbp away","10Kbp away"))

sims$pos_label <- abs(sims$pos-12500.5)
sims$pos_label <- factor( sims$pos_label, c(0,5000,10000), labels = c("Within the functional element", "5Kbp away","10Kbp away"))

ggplot(data = sims, aes( x= (gen-100000)/(10000), y = pi, group = R))+
  facet_grid(.~pos_label)+
  geom_line(data= BGS, aes(x = gen, y = B*0.01, col = as.factor(R)), lty =2)+
  scale_y_continuous(expression(pi))+
  scale_x_continuous("Time since recombination rate changed (2N generations)")+
  scale_color_brewer(expression(lambda), palette= "Dark2")+
  stat_summary(fun = "mean", aes(colour = as.factor(R)), size = 2, geom = "point")+
  theme_bw()+
  theme(
    axis.title.y =element_text(size = 20, angle = 0, vjust = 0.5),
    legend.title =element_text(size = 20),
    legend.text =element_text(size = 14)
  )


### Fixed fitness effect

B <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/TheoreticalExpectation/B_fixed_s_0.01.csv")
sims_fixed <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations/summary_simple_map_constant_s0.01.csv")
names(sims_fixed) <- c("pi","pos","gen","rep","R")

B_fixed_plot <- ggplot(data = sims_fixed[sims_fixed$gen>150000,], aes( x = pos/1000, y = pi/0.01, col = as.factor(R)))+
  stat_summary(fun = "mean", aes(colour = as.factor(R)), size = 2, geom = "point")+
  stat_summary(fun.data = mean_se, aes(colour = as.factor(R)), geom = "errorbar", width = 0)+
  geom_line(data = B, aes( x = pos/1000, y = B, col = as.factor(R)))+
  scale_y_continuous(expression(italic("B")))+
  scale_x_continuous("Position in chromosome (Kbp)")+
  coord_cartesian(ylim=c(0.5,1)) +
  scale_color_brewer(expression(lambda), palette= "Dark2")+
  theme_bw()+
    theme(
      axis.title.y =element_text(size = 20, angle = 0, vjust = 0.5),
      axis.title.x =element_text(size = 17, angle = 0, vjust = 0.5),
      legend.title =element_text(size = 20),
      strip.background = element_blank(),
      legend.text =element_text(size = 14)
    )

pdf("~/UBC/Mouse_BGS_RecombinationRateEvolution/TheoreticalExpectation/B_fixed_plot.pdf", width = 6, height = 5)
print(B_fixed_plot)
dev.off()

B$pos_label <- abs(B$pos-12500.5)

B1_mid <- mean(B[(B$R==1)&(B$pos_label<2500),]$B)
B1_near <- mean(B[(B$R==1)&(B$pos_label>2500)&(B$pos_label<=7500),]$B)
B1_far <- mean(B[(B$R==1)&(B$pos_label>7500)&(B$pos_label<=12500),]$B)

B2_mid_R01 <- mean(B[(B$R==0.1)&(B$pos_label<2500),]$B)
B2_near_R01 <- mean(B[(B$R==0.1)&(B$pos_label>2500)&(B$pos_label<=7500),]$B)
B2_far_R01 <- mean(B[(B$R==0.1)&(B$pos_label>7500)&(B$pos_label<=12500),]$B)

B2_mid_R10 <- mean(B[(B$R==10)&(B$pos_label<2500),]$B)
B2_near_R10 <- mean(B[(B$R==10)&(B$pos_label>2500)&(B$pos_label<=7500),]$B)
B2_far_R10 <- mean(B[(B$R==10)&(B$pos_label>7500)&(B$pos_label<=12500),]$B)

BGS_DF_part1 <- data.frame( B = B2_mid_R01*(1 + (B1_mid/B2_mid_R01 - 1)*exp(-1*T_0)) , gen = T_0, R = 0.1, pos = 12500.5)
BGS_DF_part2 <- data.frame( B = B1_mid*(1 + (B1_mid/B1_mid - 1)*exp(-1*T_0)) , gen = T_0, R = 1, pos = 12500.5)
BGS_DF_part3 <- data.frame( B = B2_mid_R10*(1 + (B1_mid/B2_mid_R10 - 1)*exp(-1*T_0)) , gen = T_0, R = 10, pos = 12500.5)

BGS_DF_part4 <- data.frame( B = B2_near_R01*(1 + (B1_near/B2_near_R01 - 1)*exp(-1*T_0)) , gen = T_0, R = 0.1, pos = 17500.5)
BGS_DF_part5 <- data.frame( B = B1_near*(1 + (B1_near/B1_near - 1)*exp(-1*T_0)) , gen = T_0, R = 1, pos = 17500.5)
BGS_DF_part6 <- data.frame( B = B2_near_R10*(1 + (B1_near/B2_near_R10 - 1)*exp(-1*T_0)) , gen = T_0, R = 10, pos = 17500.5)

BGS_DF_part7 <- data.frame( B = B2_far_R01*(1 + (B1_far/B2_far_R01 - 1)*exp(-1*T_0)) , gen = T_0, R = 0.1, pos = 22500.5)
BGS_DF_part8 <- data.frame( B = B1_far*(1 + (B1_far/B1_far - 1)*exp(-1*T_0)) , gen = T_0, R = 1, pos = 22500.5)
BGS_DF_part9 <- data.frame( B = B2_far_R10*(1 + (B1_far/B2_far_R10 - 1)*exp(-1*T_0)) , gen = T_0, R = 10, pos = 22500.5)


BGS <- rbind(BGS_DF_part1,BGS_DF_part2,BGS_DF_part3,
             BGS_DF_part4,BGS_DF_part5,BGS_DF_part6,
             BGS_DF_part7,BGS_DF_part8,BGS_DF_part9)
BGS$pos_label <- abs(BGS$pos-12500.5)
BGS$pos_label <- factor( BGS$pos_label, c(0,5000,10000), labels = c("Within the functional element", "5Kbp away","10Kbp away"))

sims_fixed <- sims_fixed[sims_fixed$pos%in%c(12500.5, 17500.5, 22500.5),]
sims_fixed$pos_label <- abs(sims_fixed$pos-12500.5)
sims_fixed$pos_label <- factor( sims_fixed$pos_label, c(0,5000,10000), labels = c("Within the functional element", "5Kbp away","10Kbp away"))

BGS_fixed_s_time_plot <- ggplot(data = sims_fixed, aes( x= (gen-100000)/(10000), y = pi, group = R))+
  facet_grid(.~pos_label)+
  geom_line(data= BGS, aes(x = gen, y = B*0.01, col = as.factor(R)), lty =2)+
  scale_y_continuous(expression(pi))+
  scale_x_continuous("Time since recombination rate changed (2N generations)")+
  scale_color_brewer(expression(lambda), palette= "Dark2")+
  stat_summary(fun = "mean", aes(colour = as.factor(R)), size = 2, geom = "point")+
  theme_bw()+
  theme(
    axis.title.y =element_text(size = 20, angle = 0, vjust = 0.5),
    axis.title.x =element_text(size = 20, angle = 0, vjust = 0.5),
    legend.title =element_text(size = 20),
    legend.text =element_text(size = 14),
    strip.text = element_text(size = 13),
    strip.background = element_blank(),
    
  )

pdf("~/UBC/Mouse_BGS_RecombinationRateEvolution/TheoreticalExpectation/B_over_time_fixed_s_plot.pdf", width = 10, height = 5)
print(BGS_fixed_s_time_plot)
dev.off()

### Single panel

BGS_fixed_s_time_plot_singlePanel <- ggplot(data = sims_fixed[sims_fixed$pos_label == "10Kbp away",], aes( x= (gen-100000)/(10000), y = pi, group = R))+
  geom_line(data= BGS[BGS$pos_label == "10Kbp away",], aes(x = gen, y = B*0.01, col = as.factor(R)), lty =2)+
  scale_y_continuous(expression(pi))+
  scale_x_continuous(expression("Time since recombination rate changed (2N"[e]*" generations)"))+
  scale_color_brewer(expression(lambda), palette= "Dark2")+
  stat_summary(fun = "mean", aes(colour = as.factor(R)), size = 2, geom = "point")+
  theme_bw()+
  theme(
    axis.title.y =element_text(size = 15, angle = 0, vjust = 0.5),
    axis.title.x =element_text(size = 12, angle = 0),
    axis.text = element_text(color = "black"),
    legend.title =element_text(size = 15),
    legend.text =element_text(size = 13),
    strip.text = element_text(size = 13),
    strip.background = element_blank(),
    
  )

pdf("~/UBC/Mouse_BGS_RecombinationRateEvolution/TheoreticalExpectation/B_over_time_fixed_s_plot_singlePanel.pdf", width = 6, height = 5)
print(BGS_fixed_s_time_plot_singlePanel)
dev.off()

recRates <- read.csv("~/UBC/Mouse_BGS_RecombinationRateEvolution/Simulations/simple_simResults_N1000_fixedDFE_recMaps/simple_simResults_N1000_fixedDFE_recMaps.csv")

R = c(0.1, 1.0, 10.0)
rec_rates = 2.5e-6*R*4000

recRateDF <- data.frame(R = R, rec_rates = rec_rates)

recRates_plot <- ggplot(data = recRates, aes( x= (generation-15000)/(2000), y = r_rate*4000, col = as.factor(R)))+
  geom_jitter(width = 0.01, height = 0)+
  geom_smooth(span=0.2)+
  scale_y_sqrt(expression("Population scaled recombination rate ("*italic(rho*" = 4N"[e]*"r")*")"),
               breaks = c(0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5))+
  scale_x_continuous("Time since recombination rate changed (2N generations)")+
  geom_hline(data = recRateDF, aes(yintercept = rec_rates, col = as.factor(R)),lty = 2)+
  scale_color_brewer(expression(lambda), palette= "Dark2")+
  theme_bw()+
  theme(
    axis.title.y =element_text(size = 20, vjust = 0.5),
    axis.title.x =element_text(size = 20, angle = 0, vjust = 0.5),
    legend.title =element_text(size = 20),
    legend.text =element_text(size = 14),
    strip.text = element_text(size = 13),
    strip.background = element_blank()
  )

pdf("~/UBC/Mouse_BGS_RecombinationRateEvolution/TheoreticalExpectation/recRatesPlot.pdf", width = 6, height = 5)
print(recRates_plot)
dev.off()  

### For presentation
B_fixed_plot_1 <- ggplot(data = sims_fixed[sims_fixed$gen>26000,], aes( x = pos/5000, y = pi, col = as.factor(R)))+
#  stat_summary(fun = "mean", aes(colour = as.factor(R)), size = 2, geom = "point")+
#  stat_summary(fun.data = mean_se, aes(colour = as.factor(R)), geom = "errorbar", width = 0)+
  geom_line(data = B[B$R == 1,], aes( x = pos/1000, y = B*0.01, col = as.factor(R)), lwd = 2)+
  scale_y_continuous(expression(italic("Nucleotide Diversity ("*pi*")")))+
  scale_x_continuous("Position in chromosome (Kbp)")+
  coord_cartesian(ylim=c(0.005,0.01)) +
  scale_color_manual(expression(frac(italic(r["Before"]), italic(r["After"]))), values = c("#D95F02"))+
  theme_bw()+
  theme(
    axis.title.y =element_text(size = 20, angle = 90, vjust = 0.5),
    axis.title.x =element_text(size = 17, angle = 0, vjust = 0.5),
    legend.title =element_text(size = 20),
    strip.background = element_blank(),
    legend.text =element_text(size = 14)
  )

B_fixed_plot_2 <- ggplot(data = sims_fixed[sims_fixed$gen>26000,], aes( x = pos/5000, y = pi, col = as.factor(R)))+
  #  stat_summary(fun = "mean", aes(colour = as.factor(R)), size = 2, geom = "point")+
  #  stat_summary(fun.data = mean_se, aes(colour = as.factor(R)), geom = "errorbar", width = 0)+
  geom_line(data = B, aes( x = pos/1000, y = B*0.01, col = as.factor(R)), lwd = 2)+
  scale_y_continuous(expression(italic("Nucleotide Diversity ("*pi*")")))+
  scale_x_continuous("Position in chromosome (Kbp)")+
  coord_cartesian(ylim=c(0.005,0.01)) +
  scale_color_brewer(expression(frac(italic(r["Before"]), italic(r["After"]))), palette= "Dark2")+
  theme_bw()+
  theme(
    axis.title.y =element_text(size = 20, angle = 90, vjust = 0.5),
    axis.title.x =element_text(size = 17, angle = 0, vjust = 0.5),
    legend.title =element_text(size = 20),
    strip.background = element_blank(),
    legend.text =element_text(size = 14)
  )

pdf("~/UBC/Mouse_BGS_RecombinationRateEvolution/TheoreticalExpectation/forTalk_B1_fixed_plot.pdf", width = 6, height = 5)
print(B_fixed_plot_1)
dev.off()
pdf("~/UBC/Mouse_BGS_RecombinationRateEvolution/TheoreticalExpectation/forTalk_B2_fixed_plot.pdf", width = 6, height = 5)
print(B_fixed_plot_2)
dev.off()

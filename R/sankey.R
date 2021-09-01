setwd("~/Data/Mouse_K9ac/3.0macs2/Peak_Sort")
options(stringsAsFactors = F)

library(ggplot2)
library(reshape2)
library(ggalluvial)
library(tidyr)

#all signal
signal_binary <- read.table("./all_peak.tab",header =T,sep = "\t")
colnames(signal_binary) <- c("chr","start","end","GV","M2","sperm","PN5","e2C","l2C","C4","C8","morula","ICM","TE","ESC","TSC")
signal_binary_ICM <- signal_binary[,c(4:5,7:13,15)]
signal_binary_TE <- signal_binary[,c(4:5,7:12,14,16)]

signal_binary_ICM$id=1:nrow(signal_binary_ICM)
plot_ICM <- gather(signal_binary_ICM,period,group,-id)
plot_ICM$period <- factor(plot_ICM$period,levels = c("GV","M2","PN5","e2C","l2C","C4","C8","morula","ICM","ESC"))
plot_ICM$group <- factor(plot_ICM$group,levels = c(0,1))
ggplot(plot_ICM,
       aes(x = period, stratum = group, alluvium = id,fill = group)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_flow(stat = "flow", lode.guidance = "rightleft",color = 'lightgrey',aes.flow = "forward") +
  geom_stratum()+
  theme_classic()+
  ylab('peak number')+
  theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5))
ggsave("~/Data/Mouse_K9ac/R/plot/peak_process/alluvila_ICM.pdf")

signal_binary_TE$id=1:nrow(signal_binary_TE)
plot_TE <- gather(signal_binary_TE,period,group,-id)
plot_TE$period <- factor(plot_TE$period,levels = c("GV","M2","PN5","e2C","l2C","C4","C8","morula","TE","TSC"))
plot_TE$group <- factor(plot_TE$group,levels = c(0,1))
ggplot(plot_TE,
       aes(x = period, stratum = group, alluvium = id,fill = group)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_flow(stat = "flow", lode.guidance = "rightleft",color = 'lightgrey',aes.flow = "forward") +
  geom_stratum()+
  theme_classic()+
  ylab('peak number')+
  theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5))
ggsave("~/Data/Mouse_K9ac/R/plot/peak_process/alluvila_TE.pdf")

save.image("~/Data/Mouse_K9ac/R/peak_process/alluvial.RData")

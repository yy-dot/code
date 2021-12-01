library(ggplot2)
library(tidyr)

barplot <- function(anno, xvalue, yvalue, fillvalue, ylable) {
  p <- ggplot(data=anno) +
    geom_bar(mapping=aes(x=xvalue, y=yvalue, fill=fillvalue), stat="identity",position=position_dodge(), width = 0.7) +
    #theme_classic() +
    theme(panel.background = element_rect(fill = NA),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, vjust = 1),
          axis.text.y = element_text(size = 12, margin=margin(0,0,0,10)),
          axis.title.x=element_blank(), 
          axis.title.y = element_text(size=12),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(y= ylable) +
    scale_fill_manual(values=c("#EB6101","#94ADDA", "#FCD575", "#A69425", "#9078B6", "#B4766B")) +
    geom_hline(aes(yintercept=0))
  return(p)
}


anno <- gather(raw.data,sample,values,-annotation)
anno$sample <-factor(anno$sample,level=anno$sample)
p <- barplot(anno=anno, xvalue=anno$sample, yvalue=anno$values, fillvalue=anno$annotation, ylable='Ratio of observed to expected peaks(log2)')

library(ggplot2)
library(tidyr)

barplot <- function(anno, xvalue, yvalue, fillvalue, ylable) {
  p <- ggplot(data=anno) +
    geom_bar(mapping=aes(x=xvalue, y=yvalue, fill=fillvalue), stat="identity",position=position_dodge(), width = 0.7) + #position_stack(),position_fill()
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




#需要使bar图y轴的最小值大于0
#可以使用geom_col和coord_cartesian强制让ymin=特定值,geom_col=geom_bar(stat = "identity")
ggplot(data = qc,mapping = aes(x = sample, y = `5mCpG/all_CpG(%)`)) +
  geom_col(fill='#71A0C6',width =0.7)+
  coord_cartesian(ylim = c(73,76))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   axis.line = element_line(colour = "black"),
                   legend.key.size=unit(10, "pt"),
                   axis.text.x = element_text(angle = 90, hjust = 1))

#stat = "count"：对于geom_bar()，默认是计算每个x值的行数，不需要提供y值。此时条形高度为x值的行数，即默认值stat = "count"。
#stat = "identity"：如果设置identity，说明需要提供统计好的y值。此时条形高度为提供的数值。
#1.position_identity()：不进行任何位置调整，元素按照原始数据的位置进行排列。
#2.position_stack()：在柱状图和面积图中使用，将元素按照数值堆叠在一起。
#3.position_fill()：在柱状图和面积图中使用，将元素按照数值进行归一化，填充整个绘图区域。
#4.position_dodge()：在柱状图和面积图中使用，将元素按照分组进行分开排列，避免重叠。
#5.position_jitter()：在散点图中使用，对元素的位置进行随机抖动，避免重叠。
#6.position_nudge()：在散点图和线图中使用，对元素的位置进行微调，使其偏离原始位置。

library(ggplot2)
library(tidyr)
anno <- gather(raw.data,sample,values,-annotation)
anno$sample <-factor(anno$sample,level=anno$sample)

barplot.p <- ggplot(data=anno) +
  geom_bar(mapping=aes(x=xvalue, y=yvalue, fill=fillvalue), stat="identity",position=position_dodge(), width = 0.7) + 
  theme(panel.background = element_rect(fill = NA),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1,size = 12,color="black"),
        axis.text.y = element_text(size = 12, margin=margin(0,0,0,10),color="black"),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=14)) +
  labs(x= xlable, y=ylable, title=title) +
  scale_fill_manual(values=c("#EB6101","#94ADDA", "#FCD575", "#A69425", "#9078B6", "#B4766B")) +
  scale_y_continuous(expand = c(0, 0))#使bar图紧贴x轴 #geom_hline(aes(yintercept=0))#加入y=0的横线




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

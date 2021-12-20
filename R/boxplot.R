ggplot(data=phase, mapping = aes(x = motification, y = value,fill = factor(class)))+
          geom_boxplot(notch = T,outlier.shape = NA,linetype="dashed",position = position_dodge(width = 0.8)) +
          stat_boxplot(aes(ymin=..lower..,ymax=..upper..),size=1, outlier.shape = NA, position = position_dodge(width = 0.8))+
          stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.5,size=1, position = position_dodge(width = 0.8)) + 
          stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.5,size=1 , position = position_dodge(width = 0.8)) +
           theme(panel.background = element_rect(fill = NA),
                 axis.line = element_line(colour = "black"),
                 axis.title.x=element_blank(), 
                 axis.line.x = element_blank(),
                 axis.ticks.x = element_blank()) +
           labs(y='log2(H3K9ac RPKM at promoter + 1)') +
           labs(title = colnames(k9ac_bin.rowdata)[i])+ scale_color_manual(values = c('#ADD8E6','#0000FF'))

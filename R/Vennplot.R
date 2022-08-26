library(grid)
library(VennDiagram)
zeroone <- read.delim("/home1/yyyu/Embryo/H3K9ac/3.0macs2/com_peak_matrix/all_zeroone.matirx",header =T,sep = "\t",stringsAsFactors = F,check.names = FALSE)
zeroone$flag= rownames(zeroone)
VS1=list('ESC'=zeroone[zeroone$k9ac_ESC==1,c('flag')], "TSC"=zeroone[zeroone$k9ac_TSC==1,c('flag')])
D1<-venn.diagram(x=VS1,filename=NULL,lwd=1,lty=1,col=c('red','blue'), cat.cex=1.5, cat.dist= 0.05,margin=0.1,
                 fill=c('red','blue'),cat.col=c('red','blue'), cex=1.5, main = "Distal regions shared by ESC and TSC")#rotation.degree=180, cat.pos=c(135,225),
pdf("./important_plot/Venn_peaks_shared_by_ESC_TSC.pdf",width = 6,height = 6)
grid.newpage()
grid.draw(D1)
dev.off()

#第一种
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

#第二种
library(eulerr)
library(dplyr)

PN5_integrat_01=read.delim('PN5_integrat_01.bed',header = F) %>% dplyr::mutate(flag=row_number())
colnames(PN5_integrat_01)=c("chr","start","end","k9ac","k4me3","k27ac","flag")

PN5_VS=list("k9ac"=PN5_integrat_01[PN5_integrat_01$k9ac==1,c('flag')], 
           "k4me3"=PN5_integrat_01[PN5_integrat_01$k4me3==1,c('flag')], 
           "k27ac"=PN5_integrat_01[PN5_integrat_01$k27ac==1,c('flag')])

PN5_partion <- get.venn.partitions(PN5_VS)
PN5_combinations=c(
  "k9ac"=PN5_partion[PN5_partion["..set.."]=="(k9ac)∖(k4me3∪k27ac)","..count.."],
  "k4me3"=PN5_partion[PN5_partion["..set.."]=="(k4me3)∖(k9ac∪k27ac)","..count.."],
  "k27ac"=PN5_partion[PN5_partion["..set.."]=="(k27ac)∖(k9ac∪k4me3)","..count.."],
  "k9ac&k4me3"=PN5_partion[PN5_partion["..set.."]=="(k9ac∩k4me3)∖(k27ac)","..count.."],
  "k9ac&k27ac"=PN5_partion[PN5_partion["..set.."]=="(k9ac∩k27ac)∖(k4me3)","..count.."],
  "k4me3&k27ac"=PN5_partion[PN5_partion["..set.."]=="(k4me3∩k27ac)∖(k9ac)","..count.."],
  "k9ac&k4me3&k27ac"=PN5_partion[PN5_partion["..set.."]=="k9ac∩k4me3∩k27ac","..count.."]
)

PN5_euler=euler(PN5_combinations)
p_PN5=plot(PN5_euler, 
          quantities =list(type =c("counts")),
          fills=list(fill=c('#CE0535','#0068B7',"#FCD575"),alpha=0.5),
          edges = list(col="black",lwd=2),
          main = list(label=c("PN5 peaks shared by k9ac_k4me3_k27ac"),cex=1))

p_PN5

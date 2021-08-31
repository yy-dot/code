setwd('/home1/yyyu/Embryo/H3K9ac/k9ac_snp/summary/1.0QC/2.0tree')
rm(list = ls())
options(stringsAsFactors = F)

library(ape)
library(ggtree)
library(dendextend)

tree.data <- read.delim("/home1/yyyu/Embryo/H3K9ac/k9ac_snp/4.0correlation/merge/h3k9ac_snp.mutibw.tab", sep = '\t', check.names = F)
dim(tree.data)
tree.data <- tree.data[,c(4:length(colnames(tree.data)))]
tree.data<- tree.data[which(rowSums(tree.data)>0),]

if(F){
  colnames(tree.data)<-c("e2C.M", "e2C.P", "l2C.M", "l2C.P", "4C.M", "4C.P", "8C.M","8C.P", "GV.M", "ICM.M", "ICM.P", "M2.M", "morula.M", "morula.P", "sperm.P", "TE.M", "TE.P", "zygote.M", "zygote.P")
  colname <- c("sperm.P", "M2.M", "GV.M", "zygote.M", "zygote.P", "e2C.M", "e2C.P", "l2C.M", "l2C.P", "4C.M", "4C.P", "8C.M","8C.P", "morula.M", "morula.P", "ICM.M", "ICM.P", "TE.M", "TE.P")
  colnames(tree.data) <- factor(colnames(tree.data),levels = colname)
}

tree.cor <- cor(tree.data, method = "pearson")
tree.dist = dist(tree.cor)
tree.hc = hclust(tree.dist)
tree = as.phylo(tree.hc)#转化成phylo格式
group <- split(colnames(tree.data),c("M","P","M","P","M","P","M","P","M","M","P","M","M","P","P","M","P","M","P"))#给树枝分类
tree <- groupOTU(tree, group)

p1 <- ggtree(tree,aes(color=group))+scale_color_manual(values=c('red','black'))+
  geom_tiplab(angle=90,offset = -0.15)+
  layout_dendrogram()#+geom_text(aes(label=node))

if(F) {
  p1<-ggtree::rotate(p1,26)
  p1<-ggtree::rotate(p1,27)
  p1<-ggtree::rotate(p1,35)
}

p2 <- ggtree(tree,aes(color=group), layout = "daylight")+scale_color_manual(values=c('red','black'))+
  geom_tiplab(angle=90,offset = 0)+
  layout_dendrogram()#+geom_text(aes(label=node))

pdf("snp_tree.pdf",width=6, height=5)
plot(p1)
plot(p2)
dev.off()

#需要的包
library(pheatmap)
library(factoextra)
library(RColorBrewer)
library(dplyr)

#选取方差较大的行
data$cv=apply(data, 1, function(y){sd(y)/mean(y)})
data=data[data$cv>0.25,] %>% dplyr::select(-cv)

#kmean聚类
if(F){
  data.scale=rescale_data(data,-2.5,2.5,2)
  set.seed(3)
  clus=kmeans(as.matrix(data.scale), centers=6, nstart = 8)
  anno_col_temp=data.frame(group=clus$cluster,gene=rownames(data.scale)) %>% arrange(group)
  anno_col=data.frame(group=anno_col_temp$group,row.names = anno_col_temp$gene)
  anno_col$group=as.factor(anno_col$group)
  table(anno_col$group)
}


#kmeans聚类合理性验证
if(F) {
  k.max=15
  wss = sapply(1:k.max,
              function(k){kmeans(data.scale,k,nstart=10)$tot.withinss})
  plot(1:k.max, wss,
      type='b', pch=19, frame=FALSE,
      xlab = 'Number of cluster K',
      ylab = "Total within-clusters sum of squares")

  fviz_cluster(clus, data=data.scale, geom='point', stand=FALSE, ellipse.type='norm')
}


#聚类颜色设置
group_color=colorRampPalette(rev(brewer.pal(n = 12, name ="Paired")))(20)
names(group_color) <- unique(anno_col$group)

#画热图
data.sort=data.scale[rownames(anno_col), ]#让n的行索引与anno_col一致
p1=pheatmap(data.sort,annotation_row = anno_col, 
            border_color=NA, cluster_row=F, cluster_col=F,show_rownames=F,show_colnames = T,
            color = colorRampPalette(c("#0068B7", "white", "#C9121F"))(20),
            annotation_colors=list(group=group_color))

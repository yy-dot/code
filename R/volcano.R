resdata = read.delim('/home1/yyyu/m6A/ALKBH5_knockdown_RNA/ALKBH5_seq/company_result/3_Differential_result/deglist/siALKBH5vssiCtrl_deg.xls',check.names = F)
resdata = na.omit(resdata)
resdata$threshold = factor(ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange) >= 1, ifelse(resdata$log2FoldChange>= 1,'Up','Down'),'NoSignifi'))

vol=ggplot(resdata,aes(x=log2FoldChange,y=-log10(padj),color=threshold))+
  geom_point()+ 
  scale_color_manual(values=c("#0068B7","#7B7C7D","#CF0141"))+#确定点的颜色#FCD575
  geom_text_repel(
    data = resdata[resdata$padj<0.01 & abs(resdata$log2FoldChange)>2,],
    aes(label = gene_name),
    size = 3,
    segment.color = "black", show.legend = FALSE )+ #添加关注的点的基因名
  theme_bw()+#修改图片背景
  theme(
    legend.title = element_blank(),#不显示图例标题
    plot.title = element_text(hjust = 0.5))+
  labs(title="siALKBH5 vs siCtrl",x='log2 FoldChange(Gene Expression)',y='-log10 p-adj')+
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05

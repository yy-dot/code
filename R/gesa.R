setwd("/home1/jche/gao_lab/10.gsea/gsea_output/Epi_d5.5_invivo_invitro")
rm(list = ls())
options(stringsAsFactors = F)

library(tidyverse)
library(DESeq2)
library(parallel)
library(msigdbr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)


compare.list<-list.files(path = "/home1/jche/gao_lab/5.r/SMART-seq/DEseq2/stringtie/same_celltype_diff_time/rds_output",
                         pattern = "*invivo*",
                         full.names = T)
lapply(compare.list, function(x){
  dds<-readRDS(x)
  return(resultsNames(dds))
})

do_gsea<-function(dds_rds,padj=0.5){
  ###将名称同一为ENTREZID
  dds<-readRDS(dds_rds)
  res<-as.data.frame(results(dds)) %>% arrange(desc(log2FoldChange))
  genes<-res$log2FoldChange
  names(genes)<-rownames(res)
  ids<-bitr(geneID = names(genes),fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Mm.eg.db")
  genes<-genes[ids$SYMBOL]
  names(genes)<-ids$ENTREZID
  

  ### WikiPathways
  wiki<-read.gmt("/home1/jche/gao_lab/10.gsea/mmu_gmt/wikipathways-20210310-gmt-Mus_musculus.gmt")
  wiki <- wiki %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
  wpid2gene <- wiki %>% dplyr::select(wpid, gene) #TERM2GENE
  wpid2name <- wiki %>% dplyr::select(wpid, name) #TERM2NAME
  wiki_gsea <- GSEA(genes, TERM2GENE = wpid2gene, TERM2NAME = wpid2name,pvalueCutoff = padj)
  
  ### MSigDb
  #msigdbr_species()
  msigdb <- msigdbr(species = "Mus musculus") %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame
  #dim(msigdb)
  msigdb_gsea <- GSEA(genes, TERM2GENE = msigdb,pvalueCutoff = padj)
  #nrow(em2)
  
  ### GO
  go_cc_gsea <- gseGO(geneList     = genes,
                OrgDb        = org.Mm.eg.db,
                ont          = "CC",
                pvalueCutoff = padj)
  go_bp_gsea <- gseGO(geneList     = genes,
                OrgDb        = org.Mm.eg.db,
                ont          = "BP",
                pvalueCutoff = padj)
  go_mf_gsea <- gseGO(geneList     = genes,
                OrgDb        = org.Mm.eg.db,
                ont          = "MF",
                pvalueCutoff = padj)
  
  ### KEGG
  #search_kegg_organism('mmu', by='kegg_code')
  kegg_gsea <- gseKEGG(geneList     = genes,
                  organism     = 'mmu',
                  pvalueCutoff = padj)
  mkegg_gsea <- gseMKEGG(geneList = genes,
                    organism = 'mmu',
                    pvalueCutoff = padj)
  
  ### Reactome
  reactome_gsea<-ReactomePA::gsePathway(genes,
                                   organism = "mouse",
                                   pvalueCutoff = padj)
  return(list("wiki_gsea"=wiki_gsea,
              "msigdb_gsea"=msigdb_gsea,
              "go_cc_gsea"=go_cc_gsea,
              "go_bp_gsea"=go_bp_gsea,
              "go_mf_gsea"=go_mf_gsea,
              "kegg_gsea"=kegg_gsea,
              "mkegg_gsea"=mkegg_gsea,
              "reactome_gsea"=reactome_gsea))
}

compare.list

### Epi.5.5_invivo_vs_Epi.5.5_invitro
total_result<-do_gsea(compare.list[3],padj = 1)
names(total_result)
sapply(total_result, nrow)
saveRDS(total_result,"gsea_result.rds")

invivo<-lapply(names(total_result), function(i){
  tmp<-setReadable(total_result[[i]], OrgDb = org.Mm.eg.db, keyType="ENTREZID")
  tmp<-as.data.frame(tmp) %>% filter(NES > 1) %>% arrange(qvalues) %>% top_n(30) %>% mutate(class=i)
  if(i=="msigdb_gsea"){
    tmp$Description<-tolower(tmp$Description)
  }
  return(tmp)
})
names(invivo)<-names(total_result)
sapply(invivo, nrow)
invivo<-do.call(rbind,invivo)
write_excel_csv(invivo,path = "invivo_gesa.csv")

invitro<-lapply(names(total_result), function(i){
  tmp<-setReadable(total_result[[i]], OrgDb = org.Mm.eg.db, keyType="ENTREZID")
  tmp<-as.data.frame(tmp) %>% filter(NES < -1) %>% arrange(qvalues) %>% top_n(30) %>% mutate(class=i)
  if(i=="msigdb_gsea"){
    tmp$Description<-tolower(tmp$Description)
  }
  return(tmp)
})
names(invitro)<-names(total_result)
sapply(invitro, nrow)
invitro<-do.call(rbind,invitro)
write_excel_csv(invitro,path = "invitro_gesa.csv")


if(F){
  ## invivo gsea plot
  names(total_result)
  p1<-gseaplot2(total_result[["wiki_gsea"]],
                geneSetID = c("WP113"),
                #title = invivo[invivo$ID=="WP113","Description"],
                pvalue_table = T)
  #ggsave(filename = "gsea_output/EPC.5.5_invivo_vs_Exeup.5.5_invitro/invivo_gsea-TGF_pathway.pdf",width = 8,height = 6)
  names(total_result)
  p2<-gseaplot2(total_result[["msigdb_gsea"]],
                geneSetID = c("GO_ANTERIOR_POSTERIOR_PATTERN_SPECIFICATION"),
                #title = invivo[invivo$ID=="GO_HISTONE_H3_K4_TRIMETHYLATION","Description"],
                pvalue_table = T)
  #ggsave(filename = "gsea_output/EPC.5.5_invivo_vs_Exeup.5.5_invitro/gsea-TGF_pathway.pdf",width = 8,height = 6)
  p3<-gseaplot2(total_result[["go_cc_gsea"]],
                geneSetID = c("GO:0005911"),
                #title = invivo[invivo$ID=="GO_HISTONE_H3_K4_TRIMETHYLATION","Description"],
                pvalue_table = T)
  p4<-gseaplot2(total_result[["go_bp_gsea"]],
                geneSetID = c("GO:0003002"),
                #title = invivo[invivo$ID=="GO_HISTONE_H3_K4_TRIMETHYLATION","Description"],
                pvalue_table = T)
  p5<-gseaplot2(total_result[["go_bp_gsea"]],
                geneSetID = c("GO:0045217"),
                #title = invivo[invivo$ID=="GO_HISTONE_H3_K4_TRIMETHYLATION","Description"],
                pvalue_table = T)
  p6<-gseaplot2(total_result[["go_bp_gsea"]],
                geneSetID = c("GO:0080009"),
                #title = invivo[invivo$ID=="GO_HISTONE_H3_K4_TRIMETHYLATION","Description"],
                pvalue_table = T)
}
if(F){
  invivo<-genes[genes>0]
  invitro<-genes[genes<0]
  ### Cell Marker
  cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Mouse_cell_markers.txt') %>%
    tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
    dplyr::select(cellMarker, geneID) %>%
    dplyr::mutate(geneID = strsplit(geneID, ', '))
  cell_markers
  invivo_cell_markers <- enricher(names(invivo), TERM2GENE=cell_markers, minGSSize=1)
  as.data.frame(invivo_cell_markers)$Description
  invitro_cell_markers <- enricher(names(invitro), TERM2GENE=cell_markers, minGSSize=1)
  as.data.frame(invitro_cell_markers)$Description
  
}

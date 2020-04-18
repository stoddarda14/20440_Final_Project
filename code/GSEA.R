#GSEA script written by Amy Stoddard 4/15/2020
#last modified 4/18/2020
#REMEMBER TO MODIFY FILE NAMES!

#call libraries or install packages
#####
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Mm.eg.db")
library(Seurat)
library(dplyr)
library(Matrix)
library(Seurat)
library(tidyverse)
library(ggplot2)
library("clusterProfiler")
library(org.Mm.eg.db)
library(ReactomePA)
#####

#Prep data function
#####
prep_data <- function(input_file_name)
{
  #read in file
  dge<-read.delim(input_file, header = TRUE, sep = "")
  
  # filter for pvalue of interest
  dge_screened <- dge %>% filter(p_val_adj < .05)
  #tail(dge_screened)
  
  #get conversion to ENTREZID (necessary for GSEA functions)
  genes<-dge_screened['Gene']
  char_genes<-as.character(genes[,1])
  entrez_id <- bitr(char_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  dim(entrez_id)
  
  #remove genes with no entrez id, add entrez ID to dge_filtered DF
  dge_filtered<-subset(dge,Gene %in% entrez_id[,1])
  dge_filtered$EntrezID <-entrez_id[,2]
  
  return(dge_filtered)
}

#####

#set files names and pwd, call prep data
#####
setwd("C:/Users/Amy/Documents/github/20440_Final_Project")
input_file <- "data/DGE_lists/Liver_Cluster_1_Markers_Res_0.25.txt"
input_file_2 <- "data/DGE_lists/Liver_Cluster_2_Markers_Res_0.25.txt"

dge_filtered_1 <-prep_data(input_file)
dge_filtered_2 <-prep_data(input_file_2)

#####

#Run hypergeometric enrichment analysis and plot / save
#####
GSEA_genes<-dge_filtered_1$EntrezID 
#(GSEA_genes)

x <- enrichPathway(gene=GSEA_genes,organism = "mouse", pvalueCutoff=0.01, qvalueCutoff = 0.10, readable=T)
head(as.data.frame(x))

#Plot enrichment results
pdf(file = "figures/GSEA/15_barplot_SI.pdf", width = 20, height = 15)
barplot(x,showCategory=15)
dev.off()

pdf(file = "figures/GSEA/15_dotplot_SI.pdf", width = 20, height = 15)
dotplot(x,showCategory=15)
dev.off()

pdf(file = "figures/GSEA/15_emapplot_SI.pdf", width = 20, height = 15)
emapplot(x)
dev.off()

print(x@result$ID[1])
print(x@result$Description[1])
print(x@result$geneID[1])

write.csv(x@result, file = 'data/GSEA_output/test.csv', row.names = FALSE, col.names = TRUE)

#####

#run gene ontology and plot / save
#####
#run gene ontology
y <- enrichGO(gene=GSEA_genes, OrgDb = 'org.Mm.eg.db', pvalueCutoff=0.01, qvalueCutoff=0.10, readable=T, ont="BP")
head(as.data.frame(y))

#Plot GO results
pdf(file = "figures/GSEA/15_barplot_GO_Liver_1.pdf", width = 20, height = 15)
barplot(y,showCategory=15)
dev.off()

pdf(file = "figures/GSEA/15_dotplot_GO_Liver_1.pdf", width = 20, height = 15)
dotplot(y,showCategory=15)
dev.off()

pdf(file = "figures/GSEA/15_emapplot_GO_Liver_1.pdf", width = 20, height = 15)
emapplot(y)
dev.off()

write.csv(y@result, file = 'data/GSEA_output/GO.csv', row.names = FALSE, col.names = TRUE)

#####


#run comparative enrichment analysis and plot
#####
cluster_compare <- compareCluster(geneCluster=list("1"=dge_filtered_1$EntrezID, "2"=dge_filtered_2$EntrezID), fun= "enrichPathway", organism = "mouse", pvalueCutoff=0.05, readable=T)

head(cluster_compare)
pdf(file = "figures/GSEA/15_dotplot_liver_1_2_compare.pdf", width = 15, height = 10)
dotplot(cluster_compare,showCategory=15)
dev.off()

#####

#run comparative GO analysis and plot
#####
cluster_compare_GO <- compareCluster(geneCluster=list("1"=dge_filtered_1$EntrezID, "2"=dge_filtered$EntrezID), fun= "enrichGO", OrgDb = 'org.Mm.eg.db', pvalueCutoff=0.05, readable=T)

pdf(file = "figures/GSEA/15_dotplot_liver_1_2_compare_GO.pdf", width = 15, height = 10)
dotplot(cluster_compare_GO,showCategory=15)
dev.off()
#####





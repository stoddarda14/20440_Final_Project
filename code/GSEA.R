#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("clusterProfiler")
library(Seurat)
library(dplyr)
library(Matrix)
library(Seurat)
library(tidyverse)
library(ggplot2)
library("clusterProfiler")
#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
library(ReactomePA)

setwd("C:/Users/Amy/Documents/github/20440_Final_Project")
input_file <- "data/DGE_lists/SI_Cluster_0_Markers.txt"


dge<-read.delim(input_file, header = TRUE, sep = "")

# filter for pvalue of interest
genes<-dge['Gene']
char_genes<-as.character(genes[,1])
print(char_genes)
entrez_id <- bitr(char_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dim(entrez_id)

dge_filtered<-subset(dge,Gene %in% entrez_id[,1])
dge_filtered$EntrezID <-entrez_id[,2]

GSEA_genes<-dge_filtered$EntrezID 
head(GSEA_genes)

#Run hypergepmetric enrichment analysis
x <- enrichPathway(gene=GSEA_genes,organism = "mouse", pvalueCutoff=0.01, qvalueCutoff = 0.10, readable=T)
head(as.data.frame(x))

pdf(file = "figures/GSEA/15_barplot_SI.pdf", width = 20, height = 15)
barplot(x,showCategory=15)
dev.off()

pdf(file = "figures/GSEA/15_dotplot_SI.pdf", width = 20, height = 15)
dotplot(x,showCategory=15)
dev.off()

pdf(file = "figures/GSEA/15_emapplot_SI.pdf", width = 20, height = 15)
emapplot(x)
dev.off()

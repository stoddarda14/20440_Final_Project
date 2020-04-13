# Script Use: 20.440 Team Project 
# Team: Amanda Hornick and Amy Stoddard 
# Author: Amanda Hornick, modified by AES for liver 
# Date Last Edited: 3/29/2020
# Script Description: Analylecs of single cell RNA-seq data for murine liver endothelial cells 


# Load Seurat library for single cell sequencing data analytics
library(Seurat)
library(ggplot2)

#file handling
setwd("C:/Users/Amy/Documents/github/first_repo")
input_file <- "data/data_liver.csv"
output_file <- "figures/liver_clustering.pdf"
output_file_2 <- "figures/cluster_expression.pdf"

# Load data 
lec.matrix <- read.csv(input_file, sep = ",", header = TRUE, row.names = NULL)


# Format data to have appropriate row and column names
# Source: https://github.com/satijalab/seurat/issues/1710
lec.matrix.names <- make.unique(as.character(lec.matrix$Feature))
rownames(lec.matrix) <- lec.matrix.names
lec.matrix <- lec.matrix[,-1] # Eliminates column of row names 


# Create Seurat Object 
lec <- CreateSeuratObject(counts = lec.matrix)
lec


# Perform standard log-normalization
# Source:https://satijalab.org/seurat/v3.0/mca.html
#lec <- NormalizeData(lec, normalization.method = "LogNormalize", scale.factor = 10000)


# Calculate variance and mean for each gene in dataset 
# Selecting highly variable genes based on variance mean ratio is a good strategy 
# Select the top 1000 highly variable genes for downstream analylecs
lec <- FindVariableFeatures(lec)


# Calculate and regress out mitochondrial expreslecon per cell 
lec[["percent.mt"]] <- PercentageFeatureSet(lec, pattern = "^mt-")
lec <- ScaleData(lec, vars.to.regress = "percent.mt")


# Dimenleconal reduction (PCA)
# npcs = total # of PCs to compute and store 
# ndims.print = PCs to print genes for 
# nfeatures.print = # genes to print for each PC
lec <- RunPCA(lec, npcs = 50, ndims.print = 1:15, nfeatures.print = 5)


# Visualize dimenleconality reduction with an elbow plot 
#pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/lec_elbow_plot.pdf", width = 5, height = 4) # Save plot to a PDF file 
ElbowPlot(lec, ndims = 50)
#dev.off() # Saves plot 


# Make a heat map 
# dims = dimenlecons to plot
# cells	= a list of cells to plot; if numeric, just plots the top cells
# balanced = plots an equal number of genes with both + and - scores if TRUE
# Plot 15 dimenlecons, elbow appeared to be at 10 dimenlecons
#pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/lec_heatmap.pdf", width = 10, height = 10)
DimHeatmap(lec, dims = c(1:15), cells = 500, balanced = TRUE) 
#dev.off() # Saves plot to location given above 


# Graph-based clustering 
# nn.eps = error bound when performing nearest neighbor search ulecng RANN (default = 0.0, implies exact nearest neighbor search)
# An approx nearest neighbor search increases speed by increalecng the nn.eps parameter, setting at 0 is an exact search
lec <- FindNeighbors(lec, reduction = "pca", dims = 1:15, nn.eps = 0)
# resolution = use a value above 1.0 if you want to obtain a larger number of communities and below 1.0 if you want to obtain a smaller number of communities
# n.start = number of random starts, lowering this reduces clustering time, default is to perform 100 and select the result with the highest modularity
lec <- FindClusters(lec, resolution = .5, n.start = 100) #adjust from 3 -> .5 to get less clusters


# Visualization of clustering with UMAP and tSNE
# min.dist = controls how tightly embedding is allowed to compress points together, larger values ensures embedded points are more evenly distributed, smaller values allow the algorithm to optimize more accurately with regard to local structure, values between 0.001 and 0.5, determines how clumped embedded points are
lec <- RunUMAP(lec, dims = 1:15, min.dist = 0.5)

pdf(output_file, width = 7.5, height = 6)
DimPlot(lec, reduction = "umap", pt.size = 1) + ggtitle(label = "UMAP")
dev.off() # Saves plot to location above 

# Source: https://satijalab.org/seurat/v3.0/mca.html


# Source: https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html

# Discovery of differentially expressed features (cluster biomarkers)
# Find all markers of clusters
# min.pct = only test genes that are detected in min fraction of min.pct cell, meant to speed up the function, default = 0.1 
cluster_nums <- c(0:12)
for (val in cluster_nums)
{
  cluster_markers <- FindMarkers(lec, ident.1 = val, min.pct = 0.25)
  cat(sprintf("Cluster %i", val))
  cat("\n")
  print(cluster_markers[1:5,])
}
# In the future, could save all markers for each cluster and look for stem cell supporting factors. 


# Visualization of genes in specific clusters 
# Choose a few stem cell supporting factors 
pdf(output_file_2, width = 10, height = 10)
p1 <- FeaturePlot(lec, features = c("Wnt2",  "Rspo3", "Vwf", "Jag1"), reduction = "umap", pt.size = 1, cols = c("darkolivegreen1","darkolivegreen") , combine = FALSE)
CombinePlots(plots = p1)
dev.off() # Saves plot to location above 








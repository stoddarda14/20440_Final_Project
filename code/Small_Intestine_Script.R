# Script Use: 20.440 Team Project 
# Team: Amanda Hornick and Amy Stoddard 
# Author: Amanda Hornick 
# Date Last Edited: 3/28/2020
# Script Description: Analysis of single cell RNA-seq data for murine small intestine endothelial cells 


# Load Seurat library for single cell sequencing data analysis
library(Seurat)


# Load data 
si.matrix <- read.csv("C:/Users/Amanda Hornick/myfolder2/data/Small_Intestine_Data.csv", sep = ",", header = TRUE, row.names = NULL)


# Format data to have appropriate row and column names
# Source: https://github.com/satijalab/seurat/issues/1710
si.matrix.names <- make.unique(as.character(si.matrix$Feature))
rownames(si.matrix) <- si.matrix.names
si.matrix <- si.matrix[,-1] # Eliminates column of row names 


# Create Seurat Object 
si <- CreateSeuratObject(counts = si.matrix)
si


# Perform standard log-normalization
# Source:https://satijalab.org/seurat/v3.0/mca.html
si <- NormalizeData(si, normalization.method = "LogNormalize", scale.factor = 10000)


# Calculate variance and mean for each gene in dataset 
# Selecting highly variable genes based on variance mean ratio is a good strategy 
# Select the top 1000 highly variable genes for downstream analysis
si <- FindVariableFeatures(si)


# Calculate and regress out mitochondrial expression per cell 
si[["percent.mt"]] <- PercentageFeatureSet(si, pattern = "^mt-")
si <- ScaleData(si, vars.to.regress = "percent.mt")


# Dimensional reduction (PCA)
# npcs = total # of PCs to compute and store 
# ndims.print = PCs to print genes for 
# nfeatures.print = # genes to print for each PC
si <- RunPCA(si, npcs = 50, ndims.print = 1:15, nfeatures.print = 5)


# Visualize dimensionality reduction with an elbow plot 
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/si_elbow_plot.pdf", width = 5, height = 4) # Save plot to a PDF file 
ElbowPlot(si, ndims = 50)
dev.off() # Saves plot 


# Make a heat map 
# dims = dimensions to plot
# cells	= a list of cells to plot; if numeric, just plots the top cells
# balanced = plots an equal number of genes with both + and - scores if TRUE
# Plot 15 dimensions, elbow appeared to be at 10 dimensions
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/si_heatmap.pdf", width = 10, height = 10)
DimHeatmap(si, dims = c(1:15), cells = 500, balanced = TRUE) 
dev.off() # Saves plot to location given above 


# Graph-based clustering 
# nn.eps = error bound when performing nearest neighbor search using RANN (default = 0.0, implies exact nearest neighbor search)
# An approx nearest neighbor search increases speed by increasing the nn.eps parameter, setting at 0 is an exact search
si <- FindNeighbors(si, reduction = "pca", dims = 1:15, nn.eps = 0)
# resolution = use a value above 1.0 if you want to obtain a larger number of communities and below 1.0 if you want to obtain a smaller number of communities
# n.start = number of random starts, lowering this reduces clustering time, default is to perform 100 and select the result with the highest modularity
si <- FindClusters(si, resolution = 3, n.start = 100)


# Visualization of clustering with UMAP and tSNE
# min.dist = controls how tightly embedding is allowed to compress points together, larger values ensures embedded points are more evenly distributed, smaller values allow the algorithm to optimize more accurately with regard to local structure, values between 0.001 and 0.5, determines how clumped embedded points are
si <- RunUMAP(si, dims = 1:15, min.dist = 0.5)
library(ggplot2)
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/si_clusters.pdf", width = 5, height = 4)
DimPlot(si, reduction = "umap", pt.size = 1) + ggtitle(label = "UMAP")
dev.off() # Saves plot to location above 

# Source: https://satijalab.org/seurat/v3.0/mca.html


# Source: https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html

# Discovery of differentially expressed features (cluster biomarkers)
# Find all markers of clusters
# min.pct = only test genes that are detected in min fraction of min.pct cell, meant to speed up the function, default = 0.1 
cluster_nums <- c(0:12)
for (val in cluster_nums)
{
  cluster_markers <- FindMarkers(si, ident.1 = val, min.pct = 0.25)
  cat(sprintf("Cluster %i", val))
  cat("\n")
  print(cluster_markers[1:5,])
}
# In the future, could save all markers for each cluster and look for stem cell supporting factors. 


# Visualization of genes in specific clusters 
# Choose a few stem cell supporting factors 
# Features for which cells have same value of 0 as reported when attempting to plot: "Wnt4", "Wnt5a", "Wnt5b", "Wnt7a", "Wnt9b", "Wnt11","Notch3"
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/si_featureplots.pdf", width = 10, height = 10)
p1 <- FeaturePlot(si, features = c("Wnt2",  "Notch1", "Notch2", "Notch4"), reduction = "umap", pt.size = 1, cols = c("darkolivegreen1","darkolivegreen") , combine = FALSE)
CombinePlots(plots = p1)
dev.off() # Saves plot to location above 








# Script Use: 20.440 Team Project 
# Team: Amanda Hornick and Amy Stoddard 
# Author: Amanda Hornick 
# Date Last Edited: 4/16/2020
# Script Description: Analysis of single cell RNA-seq data for murine heart endothelial cells 


# Source: https://satijalab.org/seurat/v3.0/mca.html

# Load Seurat library for single cell sequencing data analysis
library(Seurat)
library(ggplot2)
library(tibble)

##########
#file handling
setwd("C:/Users/Amy/Documents/github/20440_Final_Project")
input_file <- "data/raw_datasets/data_heart.csv"
factor_input_file <- "data/raw_datasets/canonical_factors.csv"

# Load data 
heart.matrix <- read.csv(input_file, sep = ",", header = TRUE, row.names = NULL)
##########

# Format data to have appropriate row and column names
# Source: https://github.com/satijalab/seurat/issues/1710
heart.matrix.names <- make.unique(as.character(heart.matrix$Feature))
rownames(heart.matrix) <- heart.matrix.names
heart.matrix <- heart.matrix[,-1] # Eliminates column of row names 


# Create Seurat Object 
heart <- CreateSeuratObject(counts = heart.matrix)
heart


# Perform standard log-normalization
# Source:https://satijalab.org/seurat/v3.0/mca.html
heart <- NormalizeData(heart, normalization.method = "LogNormalize", scale.factor = 10000)


# Calculate variance and mean for each gene in dataset 
# Selecting highly variable genes based on variance mean ratio is a good strategy 
# Select the top 1000 highly variable genes for downstream analysis
heart <- FindVariableFeatures(heart)


# Calculate and regress out mitochondrial expression per cell 
heart[["percent.mt"]] <- PercentageFeatureSet(heart, pattern = "^mt-")
heart <- ScaleData(heart, vars.to.regress = "percent.mt")


# Dimensional reduction (PCA)
# npcs = total # of PCs to compute and store 
# ndims.print = PCs to print genes for 
# nfeatures.print = # genes to print for each PC
# Checked 50 PCs, reaching asymptote, elbow already observed 
heart <- RunPCA(heart, npcs = 50, ndims.print = 1:20, nfeatures.print = 5)


# Visualize dimensionality reduction with an elbow plot 
# Use the same number of dimensions as for PCA
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/heart_elbow_plot.pdf", width = 5, height = 4) # Save plot to a PDF file 
ElbowPlot(heart, ndims = 50, reduction = "pca")
dev.off() # Saves plot 
# Elbow appears to be at PC 15


# Make a heat map 
# dims = dimensions to plot
# cells	= a list of cells to plot; if numeric, just plots the top cells
# balanced = plots an equal number of genes with both + and - scores if TRUE
# Plot 15 dimensions, elbow appeared to be at 10 dimensions
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/heart_heatmap.pdf", width = 10, height = 10)
DimHeatmap(heart, dims = c(1:21), cells = 500, balanced = TRUE) 
dev.off() # Saves plot to location given above 
# Doesn't seem to have a strong signal after 15 PCs 


# Graph-based clustering 
library(ggplot2)
# nn.eps = error bound when performing nearest neighbor search using RANN (default = 0.0, implies exact nearest neighbor search)
# An approx nearest neighbor search increases speed by increasing the nn.eps parameter, setting at 0 is an exact search
# Use 15 PCs as determined before 
heart <- FindNeighbors(heart, reduction = "pca", dims = 1:15, nn.eps = 0)
# resolution = use a value above 1.0 if you want to obtain a larger number of communities and below 1.0 if you want to obtain a smaller number of communities
# n.start = number of random starts, lowering this reduces clustering time, default is to perform 100 and select the result with the highest modularity

# Evaluate what resolution to use by testing different ones in FindClusters and visualizing using RunUMAP
# Want smallest value sensible (0.001) for min.dist because want high accuracy with regard to local structure to 
# better visualize relationships, don't want points dispersed as much and overlapping to get even spacing

# min.dist = controls how tightly embedding is allowed to compress points together, larger values ensures embedded points 
# are more evenly distributed, smaller values allow the algorithm to optimize more accurately with regard to local structure, values between 0.001 and 0.5, determines how clumped embedded points are

# Colors correlate with very separate clusters, may be underclustered, potentially higher density areas in cluster 0 could be broken up
heart_0.0625 <- FindClusters(heart, resolution = 0.0625, n.start = 100)
heart_0.0625 <- RunUMAP(heart_0.0625, dims = 1:15, min.dist = 0.001) 
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/heart_clusters_res_0.0625.pdf", width = 5, height = 4)
DimPlot(heart_0.0625, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP Res = 0.0625")
dev.off() # Saves plot to location above 

# Pretty much the same as lower resolution 
heart_0.125 <- FindClusters(heart, resolution = 0.125, n.start = 100)
heart_0.125 <- RunUMAP(heart_0.125, dims = 1:15, min.dist = 0.001) 
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/heart_clusters_res_0.125.pdf", width = 5, height = 4)
DimPlot(heart_0.125, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP Res = 0.125")
dev.off() # Saves plot to location above 

# Seems to add clusters arbitrarily 
heart_0.25 <- FindClusters(heart, resolution = 0.25, n.start = 100)
heart_0.25 <- RunUMAP(heart_0.25, dims = 1:15, min.dist = 0.001) 
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/heart_clusters_res_0.25.pdf", width = 5, height = 4)
DimPlot(heart_0.25, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP Res = 0.25")
dev.off() # Saves plot to location above 

# Could be splitting 0 too much 
heart_0.5 <- FindClusters(heart, resolution = 0.5, n.start = 100)
heart_0.5 <- RunUMAP(heart_0.5, dims = 1:15, min.dist = 0.001) 
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/heart_clusters_res_0.5.pdf", width = 5, height = 4)
DimPlot(heart_0.5, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP Res = 0.5")
dev.off() # Saves plot to location above 

# Seems a bit much for original cluster 0, espeically between clusters 3 and 4, some overlap
heart_1.0 <- FindClusters(heart, resolution = 1.0, n.start = 100)
heart_1.0 <- RunUMAP(heart_1.0, dims = 1:15, min.dist = 0.001) 
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/heart_clusters_res_1.0.pdf", width = 5, height = 4)
DimPlot(heart_1.0, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP Res = 1.0")
dev.off() # Saves plot to location above

# Definitely overclustered as much overlapping of different clusters 
heart_2.0 <- FindClusters(heart, resolution = 2.0, n.start = 100)
heart_2.0 <- RunUMAP(heart_2.0, dims = 1:15, min.dist = 0.001) 
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/heart_clusters_res_2.0.pdf", width = 5, height = 4)
DimPlot(heart_2.0, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP Res = 2.0")
dev.off() # Saves plot to location above

# Severely overclustered, can't even tell where different clusters are since all colors overlapping 
heart_4.0 <- FindClusters(heart, resolution = 4.0, n.start = 100)
heart_4.0 <- RunUMAP(heart_4.0, dims = 1:15, min.dist = 0.001)
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/heart_clusters_res_4.0.pdf", width = 5, height = 4)
DimPlot(heart_4.0, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP Res = 4.0")
dev.off() # Saves plot to location above


# Use a resolution of 0.5 (7 clusters)
# Don't want to miss differences, can always combine clusters later 

#########
#read in canonical factors
canonical_factors <- read.csv(factor_input_file, sep = ",", header = FALSE, row.names = 1)
factor_vector = character()
n_rows <- dim(canonical_factors)[1]
n_cols <- dim(canonical_factors[1,])[2]
counter <- 1
for (i in 1:n_rows)
{
  for (j in 1:n_cols)
  {
    if (canonical_factors[i,j] != "")
    {
      factor_vector[counter] <- as.character(canonical_factors[i,j])
      counter <- (counter + 1)
    }
  }
}
#print(factor_vector)
#########


# Source: https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html

# Discovery of differentially expressed features (markers for each cluster)
# Identify markers and export to a text file
factor_intersect <- list()
cluster_nums <- c(0:6) # Numbers for clusters as determined earlier
options("max.print"=1000000) # Allows more output to text file 
for (val in cluster_nums)
{
  # min.pct = only test genes that are detected in min fraction of min.pct cell, meant to speed up the function, default = 0.1 
  cluster_markers <- FindMarkers(heart_0.5, ident.1 = val, min.pct = 0.25)
  
  #filter for adj_pvalue < .05
  cluster_markers_screened <- cluster_markers %>% rownames_to_column('gene') %>% filter(p_val_adj < .05) %>% column_to_rownames('gene')
  #get interscection with canonical factors
  factor_intersect[[val+1]]<-intersect(factor_vector,row.names(cluster_markers_screened))
  # Number each filename by cluster markers are for 
  #filename <- sprintf("C:/Users/Amanda Hornick/myfolder2/data/Heart_Cluster_%i_Markers_Res_0.5.txt", val)
  
  # Export cluster markers to text file 
  #sink(filename)
  #print(cluster_markers)
  #sink()
}
out_file1 <- "data/DGE_lists/heart_canonical.csv"
lapply(factor_intersect, function(x) write.table( data.frame(x), out_file1 , append= T, sep=',', row.names=FALSE ))








# Determine what features differentiate clusters that are close together, check whether overclustered 

# Compare clusters 0 and 1 
cluster_markers <- FindMarkers(heart_0.5, ident.1 = 0, ident.2 = 1, min.pct = 0.25)

# Export cluster markers to text file 
sink("C:/Users/Amanda Hornick/myfolder2/data/Heart_Cluster_0_vs_1_Markers_Res_0.5.txt")
print(cluster_markers)
sink()


# Compare clusters 0 and 3 
cluster_markers <- FindMarkers(heart_0.5, ident.1 = 0, ident.2 = 3, min.pct = 0.25)

# Export cluster markers to text file 
sink("C:/Users/Amanda Hornick/myfolder2/data/Heart_Cluster_0_vs_3_Markers_Res_0.5.txt")
print(cluster_markers)
sink()


# Compare clusters 0 and 5 
cluster_markers <- FindMarkers(heart_0.5, ident.1 = 0, ident.2 = 5, min.pct = 0.25)

# Export cluster markers to text file 
sink("C:/Users/Amanda Hornick/myfolder2/data/Heart_Cluster_0_vs_5_Markers_Res_0.5.txt")
print(cluster_markers)
sink()


# Compare clusters 4 and 5 
cluster_markers <- FindMarkers(heart_0.5, ident.1 = 4, ident.2 = 5, min.pct = 0.25)

# Export cluster markers to text file 
sink("C:/Users/Amanda Hornick/myfolder2/data/Heart_Cluster_4_vs_5_Markers_Res_0.5.txt")
print(cluster_markers)
sink()



# Visualization of genes that can differentiate endothelial cells 
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/heart_Prox1.pdf", width = 10, height = 10)
FeaturePlot(heart_0.5, features = "Prox1", reduction = "umap", pt.size = 1, cols = c("darkolivegreen1","darkolivegreen") , combine = FALSE)
dev.off()


#Lyve1: lymphatic 
#Prox1: lymphatic 
#Vwf: scales with size of blood vessels 
#Efnb2: arterial 
#Ephb4: venous 

# Doesn't tell us much since markers pretty evenly spread throughout clusters 





# Relabel clusters based on researching canonical markers for endothelial cell types and comparing to differentially expressed genes
new.cluster.ids <-c("Capillary 1","Artery","Vein","Arteriole","Venule","Capillary 2","Lymphatic")
names(new.cluster.ids) <- levels(heart_0.5)
heart_0.5 <- RenameIdents(heart_0.5, new.cluster.ids)

pdf(file = "figures/clustering/heart_labeled_clusters.pdf", width = 5, height = 4)
DimPlot(heart_0.5, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP Res = 0.5")
dev.off()

pdf(file = "figures/clustering/heart_labeled_on_clusters.pdf", width = 5, height = 4)
DimPlot(heart_0.5, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE, pt.size = 0.5) + NoLegend() + ggtitle(label = "UMAP Res = 0.5")
dev.off()





# Check which clusters express canonical stem cell factors 

#####
#plot canonical factors that are differentially expressed for each cluster!
canonical_factors_to_plot <- character()
for (i in 1:length(factor_intersect))
{
  canonical_factors_to_plot <- append(canonical_factors_to_plot, factor_intersect[[i]])
}
canonical_factors_to_plot <- unique(canonical_factors_to_plot)
out_file_2 <- "figures/canonical_expression/canonical_violin_heart.pdf"
pdf(file = out_file_2, width = 8, height = 8)
VlnPlot(heart_0.5, canonical_factors_to_plot)
dev.off()

out_file_3 <- "figures/canonical_expression/canonical_heatmap_heart.pdf"
pdf(file = out_file_3, width = 10, height = 4 )
DoHeatmap(object = heart_0.5,features = factor_vector, raster=FALSE) #or canonical_factors_to_plot
dev.off()


#####



# Read in canonical stem cell factors 
SC_factors <- read.csv("C:/Users/Amanda Hornick/myfolder2/data/canonical_factors_heart_edited.csv", sep = ",", header = FALSE, row.names = 1)
SC_factors

library(stringr)

# Wnt ligands
Wnt_lig <- vector(mode = "character", length = 8)

# Make vectors for each type of ligands
for (val in 1:8){
  lig <- SC_factors[1,val]
  print(as.character(lig))
  Wnt_lig[val] <- str_to_sentence(as.character(lig)) # Ligands are uppercase in canonical factors spreadsheet but not in EC dataset 
  print(Wnt_lig)
}

# All cells have same value (0) of Wnt4, Wnt5a, Wnt5b, Wnt7a, Wnt9b, Wnt11, Rspo3
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/heart_Wnt_lig_featureplots.pdf", width = 20, height = 20)
pWnt <- FeaturePlot(heart_0.5, features = Wnt_lig, reduction = "umap", pt.size = 1, cols = c("darkolivegreen1","darkolivegreen") , combine = FALSE)
CombinePlots(plots = pWnt)
dev.off() # Saves plot to location above 
# c("WNT2","WNT4","WNT5A","WNT5B","WNT7A","WNT9B","WNT11","RSPO3")




# Fgf ligands 
# All cells have same value (0) of Fgf1
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/heart_Fgf_lig_featureplots.pdf", width = 20, height = 20)
pFgf <- FeaturePlot(heart_0.5, features = "Fgf1", reduction = "umap", pt.size = 1, cols = c("steelblue1","steelblue4") , combine = FALSE)
CombinePlots(plots = pFgf)
dev.off() # Saves plot to location above 






# Notch ligands
Notch_lig <- vector(mode = "character", length = 4)

# Make vectors for each type of ligands
for (val in 1:4){
  lig <- SC_factors[3,val]
  print(as.character(lig))
  Notch_lig[val] <- str_to_sentence(as.character(lig)) # Ligands are uppercase in canonical factors spreadsheet but not in EC dataset 
  print(Notch_lig)
}


pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/heart_Notch_lig_featureplots.pdf", width = 20, height = 20)
pNotch <- FeaturePlot(heart_0.5, features = Notch_lig, reduction = "umap", pt.size = 1, cols = c("violetred1","violetred4") , combine = FALSE)
CombinePlots(plots = pNotch)
dev.off() # Saves plot to location above 






# BMP ligands
BMP_lig <- vector(mode = "character", length = 5)

# Make vectors for each type of ligands
for (val in 1:5){
  lig <- SC_factors[4,val]
  print(as.character(lig))
  BMP_lig[val] <- str_to_sentence(as.character(lig)) # Ligands are uppercase in canonical factors spreadsheet but not in EC dataset 
  print(BMP_lig)
}

# All cells have same value (0) of Bmp2, Bmp6, and Bmp7
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/heart_BMP_lig_featureplots.pdf", width = 20, height = 20)
pBMP <- FeaturePlot(heart_0.5, features = BMP_lig, reduction = "umap", pt.size = 1, cols = c("wheat1","wheat4") , combine = FALSE)
CombinePlots(plots = pBMP)
dev.off() # Saves plot to location above 






# TGFb ligands
TGFb_lig <- vector(mode = "character", length = 5)

# Make vectors for each type of ligands
for (val in 1:5){
  lig <- SC_factors[5,val]
  print(as.character(lig))
  TGFb_lig[val] <- str_to_sentence(as.character(lig)) # Ligands are uppercase in canonical factors spreadsheet but not in EC dataset 
  print(TGFb_lig)
}

# All cells have same value (0) of Inhbb and Gdf15
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/heart_TGFb_lig_featureplots.pdf", width = 20, height = 20)
pTGFb <- FeaturePlot(heart_0.5, features = TGFb_lig, reduction = "umap", pt.size = 1, cols = c("indianred1","indianred4") , combine = FALSE)
CombinePlots(plots = pTGFb)
dev.off() # Saves plot to location above 




# TGFb antagonism ligands
TGFb_ant_lig <- vector(mode = "character", length = 5)

# Make vectors for each type of ligands
for (val in 1:5){
  lig <- SC_factors[6,val]
  print(as.character(lig))
  TGFb_ant_lig[val] <- str_to_sentence(as.character(lig)) # Ligands are uppercase in canonical factors spreadsheet but not in EC dataset 
  print(TGFb_ant_lig)
}

# All cells have same value (0) of Nog
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/heart_TGFb_ant_lig_featureplots.pdf", width = 20, height = 20)
pTGFb_ant <- FeaturePlot(heart_0.5, features = TGFb_ant_lig, reduction = "umap", pt.size = 1, cols = c("mediumpurple1","mediumpurple4") , combine = FALSE)
CombinePlots(plots = pTGFb_ant)
dev.off() # Saves plot to location above 




# HH ligands 
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/heart_HH_lig_featureplots.pdf", width = 20, height = 20)
pHH <- FeaturePlot(heart_0.5, features = "Dhh", reduction = "umap", pt.size = 1, cols = c("chocolate1","chocolate4") , combine = FALSE)
CombinePlots(plots = pHH)
dev.off() # Saves plot to location above 





# EGFR ligands 
# All cells have same value (0) of Tgfa
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/heart_EGFR_lig_featureplots.pdf", width = 20, height = 20)
pEGFR <- FeaturePlot(heart_0.5, features = c("Tgfa","Hbegf"), reduction = "umap", pt.size = 1, cols = c("turquoise1","turquoise4") , combine = FALSE)
CombinePlots(plots = pEGFR)
dev.off() # Saves plot to location above 


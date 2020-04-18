# Script Use: 20.440 Team Project 
# Team: Amanda Hornick and Amy Stoddard 
# Author: Amanda Hornick 
# Date Last Edited: 4/15/2020
# Script Description: Analysis of single cell RNA-seq data for murine liver endothelial cells 


# Source: https://satijalab.org/seurat/v3.0/mca.html

# Load Seurat library for single cell sequencing data analysis
library(Seurat)
library(ggplot2)

##########
#file handling
setwd("C:/Users/Amy/Documents/github/20440_Final_Project")
input_file <- "data/raw_datasets/data_liver.csv"
factor_input_file <- "data/raw_datasets/canonical_factors.csv"

# Load data 
#liver.matrix <- read.csv("C:/Users/Amanda Hornick/myfolder2/data/Liver_Data.csv", sep = ",", header = TRUE, row.names = NULL)
liver.matrix <- read.csv(input_file, sep = ",", header = TRUE, row.names = NULL)
##########

# Format data to have appropriate row and column names
# Source: https://github.com/satijalab/seurat/issues/1710
liver.matrix.names <- make.unique(as.character(liver.matrix$Feature))
rownames(liver.matrix) <- liver.matrix.names
liver.matrix <- liver.matrix[,-1] # Eliminates column of row names 


# Create Seurat Object 
liver <- CreateSeuratObject(counts = liver.matrix)
liver


# Perform standard log-normalization
# Source:https://satijalab.org/seurat/v3.0/mca.html
liver <- NormalizeData(liver, normalization.method = "LogNormalize", scale.factor = 10000)


# Calculate variance and mean for each gene in dataset 
# Selecting highly variable genes based on variance mean ratio is a good strategy 
# Select the top 1000 highly variable genes for downstream analysis
liver <- FindVariableFeatures(liver)


# Calculate and regress out mitochondrial expression per cell 
liver[["percent.mt"]] <- PercentageFeatureSet(liver, pattern = "^mt-")
liver <- ScaleData(liver, vars.to.regress = "percent.mt")


# Dimensional reduction (PCA)
# npcs = total # of PCs to compute and store 
# ndims.print = PCs to print genes for 
# nfeatures.print = # genes to print for each PC
# Checked 50 PCs, reaching asymptote, elbow already observed 
# Checked 100 PCs, still reaching asymptote, elbow already observed
# Checked 5 and 10 PCs, elbow not as easily observable, seems to keep going
liver <- RunPCA(liver, npcs = 50, ndims.print = 1:20, nfeatures.print = 5)


# Visualize dimensionality reduction with an elbow plot 
# Use the same number of dimensions as for PCA
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/liver_elbow_plot.pdf", width = 5, height = 4) # Save plot to a PDF file 
ElbowPlot(liver, ndims = 50, reduction = "pca")
dev.off() # Saves plot 
# Elbow appears to be at PC 17 


# Make a heat map 
# dims = dimensions to plot
# cells	= a list of cells to plot; if numeric, just plots the top cells
# balanced = plots an equal number of genes with both + and - scores if TRUE
# Plot 15 dimensions, elbow appeared to be at 10 dimensions
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/liver_heatmap.pdf", width = 10, height = 10)
DimHeatmap(liver, dims = c(1:21), cells = 500, balanced = TRUE) 
dev.off() # Saves plot to location given above 
# Honestly these heatmaps are kind of useless and don't agree with the elbow method really


# Graph-based clustering 
library(ggplot2)
# nn.eps = error bound when performing nearest neighbor search using RANN (default = 0.0, implies exact nearest neighbor search)
# An approx nearest neighbor search increases speed by increasing the nn.eps parameter, setting at 0 is an exact search
# Use 17 PCs as determined before 
liver <- FindNeighbors(liver, reduction = "pca", dims = 1:17, nn.eps = 0)
# resolution = use a value above 1.0 if you want to obtain a larger number of communities and below 1.0 if you want to obtain a smaller number of communities
# n.start = number of random starts, lowering this reduces clustering time, default is to perform 100 and select the result with the highest modularity

# Evaluate what resolution to use by testing different ones in FindClusters and visualizing using RunUMAP
# Want smallest value sensible (0.001) for min.dist because want high accuracy with regard to local structure to 
# better visualize relationships, don't want points dispersed as much and overlapping to get even spacing

# min.dist = controls how tightly embedding is allowed to compress points together, larger values ensures embedded points 
# are more evenly distributed, smaller values allow the algorithm to optimize more accurately with regard to local structure, values between 0.001 and 0.5, determines how clumped embedded points are

# Definitely seems underclustered on the left side of the plot
liver_0.0625 <- FindClusters(liver, resolution = 0.0625, n.start = 100)
liver_0.0625 <- RunUMAP(liver_0.0625, dims = 1:17, min.dist = 0.001) 
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/liver_clusters_res_0.0625.pdf", width = 5, height = 4)
DimPlot(liver_0.0625, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP Res = 0.0625")
dev.off() # Saves plot to location above 

# Seems to eliminate an important difference in the middle cluster
liver_0.125 <- FindClusters(liver, resolution = 0.125, n.start = 100)
liver_0.125 <- RunUMAP(liver_0.125, dims = 1:17, min.dist = 0.001) 
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/liver_clusters_res_0.125.pdf", width = 5, height = 4)
DimPlot(liver_0.125, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP Res = 0.125")
dev.off() # Saves plot to location above 

# May be missing some variation in the largest cluster, some areas are denser than others, may be clusters
liver_0.25 <- FindClusters(liver, resolution = 0.25, n.start = 100)
liver_0.25 <- RunUMAP(liver_0.25, dims = 1:17, min.dist = 0.001) 
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/liver_clusters_res_0.25.pdf", width = 5, height = 4)
DimPlot(liver_0.25, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP Res = 0.25")
dev.off() # Saves plot to location above 

# May be losing some variation within some of the larger clusters 
liver_0.5 <- FindClusters(liver, resolution = 0.5, n.start = 100)
liver_0.5 <- RunUMAP(liver_0.5, dims = 1:17, min.dist = 0.001) 
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/liver_clusters_res_0.5.pdf", width = 5, height = 4)
DimPlot(liver_0.5, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP Res = 0.5")
dev.off() # Saves plot to location above 

# May be identifying important clusters, unclear if overclustering
liver_1.0 <- FindClusters(liver, resolution = 1.0, n.start = 100)
liver_1.0 <- RunUMAP(liver_1.0, dims = 1:17, min.dist = 0.001) 
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/liver_clusters_res_1.0.pdf", width = 5, height = 4)
DimPlot(liver_1.0, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP Res = 1.0")
dev.off() # Saves plot to location above

# Very overclustered
liver_2.0 <- FindClusters(liver, resolution = 2.0, n.start = 100)
liver_2.0 <- RunUMAP(liver_2.0, dims = 1:17, min.dist = 0.001) 
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/liver_clusters_res_2.0.pdf", width = 5, height = 4)
DimPlot(liver_2.0, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP Res = 2.0")
dev.off() # Saves plot to location above

# Severely overclustered
liver_4.0 <- FindClusters(liver, resolution = 4.0, n.start = 100)
liver_4.0 <- RunUMAP(liver_4.0, dims = 1:17, min.dist = 0.001)
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/liver_clusters_res_4.0.pdf", width = 5, height = 4)
DimPlot(liver_4.0, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP Res = 4.0")
dev.off() # Saves plot to location above







# Source: https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html

# Discovery of differentially expressed features (markers for each cluster)
# Identify markers and export to a text file 
cluster_nums <- c(0:10) # Numbers for clusters as determined earlier
options("max.print"=1000000) # Allows more output to text file 
for (val in cluster_nums)
{
  # min.pct = only test genes that are detected in min fraction of min.pct cell, meant to speed up the function, default = 0.1 
  cluster_markers <- FindMarkers(liver_1.0, ident.1 = val, min.pct = 0.25)
  
  # Number each filename by cluster markers are for 
  filename <- sprintf("C:/Users/Amanda Hornick/myfolder2/data/Liver_Cluster_%i_Markers_Res_1.0.txt", val)
  
  # Export cluster markers to text file 
  sink(filename)
  print(cluster_markers)
  sink()
}





# Redo since seemed overclustered 
# Discovery of differentially expressed features (markers for each cluster)
# Identify markers and export to a text file 
cluster_nums <- c(0:7) # Numbers for clusters as determined earlier
options("max.print"=1000000) # Allows more output to text file 
for (val in cluster_nums)
{
  # min.pct = only test genes that are detected in min fraction of min.pct cell, meant to speed up the function, default = 0.1 
  cluster_markers <- FindMarkers(liver_0.5, ident.1 = val, min.pct = 0.25)
  
  # Number each filename by cluster markers are for 
  filename <- sprintf("C:/Users/Amanda Hornick/myfolder2/data/Liver_Cluster_%i_Markers_Res_0.5.txt", val)
  
  # Export cluster markers to text file 
  sink(filename)
  print(cluster_markers)
  sink()
}
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
print(factor_vector)
#########

library(tibble)

# Redo AGAIN since seemed overclustered 
# Discovery of differentially expressed features (markers for each cluster)
# Identify markers and export to a text file 
factor_intersect <- list()
cluster_nums <- c(0:4) # Numbers for clusters as determined earlier
options("max.print"=1000000) # Allows more output to text file 
for (val in cluster_nums)
{
  # min.pct = only test genes that are detected in min fraction of min.pct cell, meant to speed up the function, default = 0.1 
  cluster_markers <- FindMarkers(liver_0.25, ident.1 = val, min.pct = 0.25)
  
  #filter for adj_pvalue < .05
  cluster_markers_screened <- cluster_markers %>% rownames_to_column('gene') %>% filter(p_val_adj < .05) %>% column_to_rownames('gene')
  #get interscection with canonical factors
  factor_intersect[[val+1]]<-intersect(factor_vector,row.names(cluster_markers_screened))
  # Number each filename by cluster markers are for 
  #filename <- sprintf("C:/Users/Amanda Hornick/myfolder2/data/Liver_Cluster_%i_Markers_Res_0.25.txt", val)
  # Export cluster markers to text file 
  #sink(filename)
  #head(cluster_markers)
  #sink()
}
#####
#save differential / canonical intersection
out_file1 <- "data/DGE_lists/liver_canonical.csv"
lapply(factor_intersect, function(x) write.table( data.frame(x), out_file1 , append= T, sep=',', row.names=FALSE ))
head(cluster_markers_screened['Top2a',])

#this doesn't work
#new <- cluster_markers_screened %>% filter(row.names(cluster_markers_screened) %in% factor_vector)
#new <- cluster_markers_screened[cluster_markers_screened$Name %in% factor_vector ]
#head(new)
#####






# Determine what features differentiate clusters that are close together, check whether overclustered 

# Compare clusters 0 and 1 
cluster_markers <- FindMarkers(liver_0.25, ident.1 = 0, ident.2 = 1, min.pct = 0.25)

# Export cluster markers to text file 
sink("C:/Users/Amanda Hornick/myfolder2/data/Liver_Cluster_0_vs_1_Markers_Res_0.25.txt")
print(cluster_markers)
sink()




# Compare clusters 0 and 2 
cluster_markers <- FindMarkers(liver_0.25, ident.1 = 0, ident.2 = 2, min.pct = 0.25)

# Export cluster markers to text file 
sink("C:/Users/Amanda Hornick/myfolder2/data/Liver_Cluster_0_vs_2_Markers_Res_0.25.txt")
print(cluster_markers)
sink()




# Compare clusters 0 and 3 
cluster_markers <- FindMarkers(liver_0.25, ident.1 = 0, ident.2 = 3, min.pct = 0.25)

# Export cluster markers to text file 
sink("C:/Users/Amanda Hornick/myfolder2/data/Liver_Cluster_0_vs_3_Markers_Res_0.25.txt")
print(cluster_markers)
sink()







# Relabel clusters based on researching canonical markers for endothelial cell types and comparing to differentially expressed genes
new.cluster.ids <-c("Artery 1","Artery 2","Venule","Sinusoidal","Proliferating")
names(new.cluster.ids) <- levels(liver_0.25)
liver_0.25 <- RenameIdents(liver_0.25, new.cluster.ids)

pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/liver_labeled_clusters_res_0.25.pdf", width = 5, height = 4)
DimPlot(liver_0.25, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP Res = 0.25")
dev.off()

pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/liver_labeled_on_clusters_res_0.25.pdf", width = 5, height = 4)
DimPlot(liver_0.25, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE, pt.size = 0.5) + NoLegend() + ggtitle(label = "UMAP Res = 0.25")
dev.off()


#####
#plot canonical factors that are differentially expressed for each cluster!
canonical_factors_to_plot <- character()
for (i in 1:length(factor_intersect))
{
  canonical_factors_to_plot <- append(canonical_factors_to_plot, factor_intersect[[i]])
}
canonical_factors_to_plot <- unique(canonical_factors_to_plot)
out_file_2 <- "figures/canonical_expression/canonical_violin_liver.pdf"
pdf(file = out_file_2, width = 10, height = 20)
VlnPlot(liver_0.25, canonical_factors_to_plot)
dev.off()

out_file_3 <- "figures/canonical_expression/canonical_heatmap_liver.pdf"
pdf(file = out_file_3, width = 10, height = 4 )
DoHeatmap(object = liver_0.25,features = factor_vector, raster=FALSE) #or canonical_factors_to_plot
dev.off()


#####



# Check which clusters express canonical stem cell factors 

# Read in canonical stem cell factors 
SC_factors <- read.csv("C:/Users/Amanda Hornick/myfolder2/data/canonical_factors_edited.csv", sep = ",", header = FALSE, row.names = 1)
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

# All cells have same value (0) of Wnt4, Wnt5a, Wnt5b, Wnt7a, Wnt11
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/liver_Wnt_lig_featureplots.pdf", width = 20, height = 20)
pWnt <- FeaturePlot(liver_0.25, features = Wnt_lig, reduction = "umap", pt.size = 1, cols = c("darkolivegreen1","darkolivegreen") , combine = FALSE)
CombinePlots(plots = pWnt)
dev.off() # Saves plot to location above 





# Fgf ligands 
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/liver_Fgf_lig_featureplots.pdf", width = 20, height = 20)
pFgf <- FeaturePlot(liver_0.25, features = "Fgf1", reduction = "umap", pt.size = 1, cols = c("steelblue1","steelblue4") , combine = FALSE)
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

# All cells have the same value (0) of Jag2
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/liver_Notch_lig_featureplots.pdf", width = 20, height = 20)
pNotch <- FeaturePlot(liver_0.25, features = Notch_lig, reduction = "umap", pt.size = 1, cols = c("violetred1","violetred4") , combine = FALSE)
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

# All cells have same value (0) of Bmp1, Bmp6, and Bmp7
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/liver_BMP_lig_featureplots.pdf", width = 20, height = 20)
pBMP <- FeaturePlot(liver_0.25, features = BMP_lig, reduction = "umap", pt.size = 1, cols = c("wheat1","wheat4") , combine = FALSE)
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

# All cells have same value (0) of Tgfb3 and Gdf15
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/liver_TGFb_lig_featureplots.pdf", width = 20, height = 20)
pTGFb <- FeaturePlot(liver_0.25, features = TGFb_lig, reduction = "umap", pt.size = 1, cols = c("indianred1","indianred4") , combine = FALSE)
CombinePlots(plots = pTGFb)
dev.off() # Saves plot to location above 




# TGFb antagonism ligands
TGFb_ant_lig <- vector(mode = "character", length = 4)

# Make vectors for each type of ligands, take out Nog because not in liver dataset 
for (val in 1:4){
  lig <- SC_factors[6,val+1]
  print(as.character(lig))
  TGFb_ant_lig[val] <- str_to_sentence(as.character(lig)) # Ligands are uppercase in canonical factors spreadsheet but not in EC dataset 
  print(TGFb_ant_lig)
}

# All cells have same value (0) of Nbl1
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/liver_TGFb_ant_lig_featureplots.pdf", width = 20, height = 20)
pTGFb_ant <- FeaturePlot(liver_0.25, features = TGFb_ant_lig, reduction = "umap", pt.size = 1, cols = c("mediumpurple1","mediumpurple4") , combine = FALSE)
CombinePlots(plots = pTGFb_ant)
dev.off() # Saves plot to location above 




# HH ligands 
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/liver_HH_lig_featureplots.pdf", width = 20, height = 20)
pHH <- FeaturePlot(liver_0.25, features = "Dhh", reduction = "umap", pt.size = 1, cols = c("chocolate1","chocolate4") , combine = FALSE)
CombinePlots(plots = pHH)
dev.off() # Saves plot to location above 





# EGFR ligands 
# All cells have same value (0) of Tgfa
pdf(file = "C:/Users/Amanda Hornick/myfolder2/plots/liver_EGFR_lig_featureplots.pdf", width = 20, height = 20)
pEGFR <- FeaturePlot(liver_0.25, features = c("Tgfa","Hbegf"), reduction = "umap", pt.size = 1, cols = c("turquoise1","turquoise4") , combine = FALSE)
CombinePlots(plots = pEGFR)
dev.off() # Saves plot to location above 


#find intersection between 









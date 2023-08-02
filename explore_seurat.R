library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(hdf5r)

# data source: https://www.10xgenomics.com/resources/datasets/5k-human-a0201-b0702-pbmcs-beam-t-2-standard
# tutorial references: https://www.youtube.com/watch?v=5HBzgsz8qyk https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

# Load the PBMC dataset
pbmc.data <- Read10X_h5('pbmc_scrnaseq/5k_BEAM-T_Human_A0201_B0702_PBMC_5pv2_Multiplex_count_raw_feature_bc_matrix.h5')
# multiple modalities present. take a look at them
str(pbmc.data)
# only care about the gene expression data
cts <- pbmc.data$`Gene Expression`
# take a look at the matrix
cts
# rows are the genes. columns are the cell names. "." represent the null value. (it is a sparse matrix)
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = cts, project = "pbmc5k", min.cells = 3, min.features = 200)
pbmc
# An object of class Seurat 
# 17735 features across 3873 samples within 1 assay 
# Active assay: RNA (17735 features, 0 variable features)
# 2 layers present: counts, data

## Standard pre-processing workflow
# 1 QC 
# metrics include: 
# 1) the number of unique genes detected in each cell (Low-quality cells or empty droplets will often have very few genes; Cell doublets or multiplets may exhibit an aberrantly high gene count)
# 2) the total number of molecules detected within a cell (correlates strongly with unique genes)
# 3) The percentage of reads that map to the mitochondrial genome. (Low-quality / dying cells often exhibit extensive mitochondrial contamination; 
# We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features; 
# We use the set of all genes starting with MT- as a set of mitochondrial genes)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Where are QC metrics stored in Seurat?
# The number of unique genes and total molecules are automatically calculated during CreateSeuratObject()
# You can find them stored in the object meta data

# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

# In the example below, we visualize QC metrics, and use these to filter cells.
# We filter cells that have unique feature counts over 2,500 or less than 200
# We filter cells that have >5% mitochondrial counts

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 2 Normalizing the data
# After removing unwanted cells from the dataset, the next step is to normalize the data. 
# By default, we employ a global-scaling normalization method “LogNormalize” that 
# normalizes the feature expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
# Normalized values are stored in pbmc[["RNA"]]@data.
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
# For clarity, in this previous line of code (and in future commands), we provide 
# the default values for certain parameters in the function call. However, this isn’t
# required and the same behavior can be achieved with:
# pbmc <- NormalizeData(pbmc)

# 3 Identification of highly variable features (feature selection)
# The highgle variable genes are useful for the downstream analysis. 
# The procedure in Seurat is described here https://doi.org/10.1016/j.cell.2019.05.031. 
# By default 2000 features will be selected. 
pbmc <- FindVariableFeatures(pbmc)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# 4 Scale the data
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# The results of this are stored in pbmc[["RNA"]]@scale.data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# 5 PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
# print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
# PC_ 1 
# Positive:  PRKCA, NEAT1, ATXN1, MALAT1, ARHGAP15 
# Negative:  EEF1A1, ACTB, S100A4, IL32, PFN1 
# PC_ 2 
# Positive:  TRBV9, TRAV38-2DV8, RAMP1, FILIP1L, NELL2 
# Negative:  TRBV6-5, TRAV24, TRGV10, HLA-DQA2, GZMH 
# PC_ 3 
# Positive:  GNLY, KLF2, GZMB, TRAV24, TRBV6-5 
# Negative:  TRBV13, TRAV41, TRAV8-6, TRGV4, TRGV2 
# PC_ 4 
# Positive:  GPI, ENO1, TPI1, PGK1, GAPDH 
# Negative:  LTB, FLNA, NFKBIA, GZMH, BIRC3 
# PC_ 5 
# Positive:  RRM2, UBE2C, GAPDH, IFITM2, MKI67 
# Negative:  TRAV41, TRBV13, TRAV8-6, TRGV5P, TRGV2 
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# 6 Determine the dimentionality of the dataset
ElbowPlot(pbmc)
# the elbow happens when the PC is at 10-11
# Alternative methods: 
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20)
# There are some fluctuations in the P-values, which makes the choice of PC values uncertain. 
# But there seems PC = 11 can be considered a cutoff value as the p-value significantly drops when PC=12. However PC=15, 6, 9 are also choices worth try. 

# 7 Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
# 
# Number of nodes: 917
# Number of edges: 31223
# 
# Running Louvain algorithm...
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#      Maximum modularity in 10 random starts: 0.7880
#    Number of communities: 7
#    Elapsed time: 0 seconds
head(Idents(pbmc), 7)

# 8 Run non-linear dimensional reduction (UMAP/tSNE)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE)
## save the output
saveRDS(pbmc, file = "pbmc_scrnaseq/pbmc5k_example.rds")

# 9 Finding differentially expressed features (cluster biomarkers)
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene       
# <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>      
#   1 1.28e- 32      0.739 1     0.952 2.27e- 28 0       S100A4     
# 2 2.27e- 13      0.644 0.82  0.56  4.02e-  9 0       TRBV6-5    
# 3 1.75e- 32      1.46  0.816 0.585 3.11e- 28 1       BCL2       
# 4 2.13e- 29      1.43  0.735 0.42  3.78e- 25 1       PRKCA      
# 5 6.33e- 18      0.759 0.831 0.472 1.12e- 13 2       TRAV24     
# 6 1.92e- 11      0.590 0.864 0.569 3.41e-  7 2       TRBV6-5    
# 7 3.33e- 67      1.94  1     0.847 5.91e- 63 3       TPI1       
# 8 8.89e- 56      1.86  0.967 0.708 1.58e- 51 3       PGK1       
# 9 1.35e-148      4.64  1     0.042 2.40e-144 4       TRBV13     
# 10 3.55e-156      3.32  0.956 0.024 6.30e-152 4       TRAV41     
# 11 6.23e-147      3.46  0.76  0.001 1.10e-142 5       TRBV9      
# 12 1.92e-115      3.39  0.64  0.003 3.41e-111 5       TRAV38-2DV8
# 13 4.71e-  4      1.24  0.167 0.759 1   e+  0 6       HSP90B1    
# 14 3.12e-  5      1.11  0.889 0.951 5.53e-  1 6       ARHGDIB    
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
VlnPlot(pbmc, features = c("TRAV38-2DV8", "TRBV9", "TRBV13", "PGK1"))
# Question: why the top ranked markers do not have a significant different expression among clusters??
VlnPlot(pbmc, features = c("TRAV41"))

FeaturePlot(pbmc, features = c("TRBV6-5", "BCL2", "TRAV24", "TPI1", "TRBV13", "TRBV9", "HSP90B1"))
# 10 Assign cell type to clusters
# How to find out the markers for each cluster?
# How to assign the cell types based on the markers? 

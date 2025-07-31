# Script to analyze scRNA-seq data from NSCLC using the Seurat workflow
# Dataset: 20k NSCLC DTCs, 3' v3.1 (from 10X Genomics)

# Install packages (run once)
install.packages("tidyverse")
install.packages("hdf5r")

# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)

# Load raw data from 10X HDF5 file
nsclc.sparse.m <- Read10X_h5(filename = 'C:/Users/mvbis/iCloudDrive/Msc Lab/Analysis/R/Practice/20k_NSCLC_DTC_3p_nextgem_intron_Multiplex_count_raw_feature_bc_matrix.h5')

# Check data structure (contains modalities: gene expression, antibody, multiplex)
str(nsclc.sparse.m)

# Extract gene expression counts
cts <- nsclc.sparse.m$`Gene Expression`

# Preview first 10 genes and cells
cts[1:10, 1:10]

# Create Seurat object (only cells with ≥200 genes and genes in ≥3 cells)
nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)

# Inspect object structure
str(nsclc.seurat.obj)
nsclc.seurat.obj

# -------- 1. Quality Control (QC) ------------------------

# Typical QC metrics:
  # nFeature_RNA: number of genes per cell
  # nCount_RNA: number of UMIs per cell
  # percent.mt: % mitochondrial gene content (cell stress/death)

# Calculate mitochondrial % using genes starting with "MT-"
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Plot: UMI counts vs number of genes (QC scatter)
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +geom_smooth(method = "lm")

# -------- 2. Filter low-quality cells --------------------

# Keep cells with:
# genes between 200–2500
# <5% mitochondrial content
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 &percent.mt < 5)

# Optional: check stats
summary(nsclc.seurat.obj$nFeature_RNA)
summary(nsclc.seurat.obj$percent.mt)
max(nsclc.seurat.obj$nFeature_RNA)
max(nsclc.seurat.obj$percent.mt)

# -------- 3. Normalize data ------------------------------

# Normalize gene expression per cell (default: LogNormalize)
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)

# -------- 4. Identify variable genes ---------------------

# Find 500 most variable genes using "vst" method
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 500)

# Top 10 most variable genes
top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)

# Plot variable genes
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# -------- 5. Scale data ----------------------------------

# Standardize expression (mean=0, variance=1)
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)

# -------- 6. PCA (Linear Dimensionality Reduction) ------

# Run PCA using variable genes
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(nsclc.seurat.obj))

# visualize PCA results
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(nsclc.seurat.obj, dims = 2, cells = 500, balanced = TRUE)
VizDimLoadings(nsclc.seurat.obj, dims = 1:2, reduction = "pca")
DimPlot(nsclc.seurat.obj, reduction = "pca")

# Choose number of PCs to use
ElbowPlot(nsclc.seurat.obj)

# -------- 7. Clustering ----------------------------------

# Build KNN graph using first 15 PCs
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)

# Find clusters at different resolutions
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1))

# Visualize clusters (resolution 0.5)
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.1", label = TRUE)

# Set cluster identity (e.g., use resolution 0.1)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"

# -------- 8. UMAP (Non-linear Dimensionality Reduction) --

# Run UMAP using first 15 PCs
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)

# Plot UMAP clusters
DimPlot(nsclc.seurat.obj, reduction = "umap", group.by = "seurat_clusters")

# Downsample and replot
Idents(nsclc.seurat.obj) <- "seurat_clusters"
subset <- subset(nsclc.seurat.obj, downsample = 1000)
DimPlot(subset, reduction = "umap")

# -------- 9. Marker Gene Detection -----------------------

# Find markers for all clusters
markers <- FindAllMarkers(subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers)

# Top 10 markers per cluster with avg_log2FC > 1
top10 <- markers %>%
  dplyr::filter(avg_log2FC > 1) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 10) %>%
  dplyr::ungroup()

# Heatmap of top markers
DoHeatmap(subset, features = top10$gene) + NoLegend()

# Find markers for a specific cluster (e.g., cluster 2)
cluster13.markers <- FindMarkers(subset, ident.1 = 13)
head(cluster13.markers, n = 5)
FeaturePlot(nsclc.seurat.obj, features = c("CD3D", "CD4", "CD8A", "GZMB", "IFNG", "NKG7", "KLRD1", "GNLY", "PRF1"))

nsclc.seurat.obj <- RenameIdents(nsclc.seurat.obj, '13' = "T CD8+ activated")

DimPlot(nsclc.seurat.obj, label = TRUE)


# Plot specific gene expression (e.g., CCL5)
VlnPlot(subset, features = c("CCL5"))



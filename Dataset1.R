# Install packages (run once)
install.packages("hdf5r")

# Load necessary libraries
library(Seurat)
library(hdf5r)
library(ggplot2)
library(patchwork)  # Needed for plot_annotation()

# Set working directory (where the data files are located
setwd("/home/mica/Desktop/Dataset_1/GSE278456_RAW_LGG")

# Define folder path containing .h5 files
folder_path <- "/home/mica/Desktop/Dataset_1/GSE278456_RAW_LGG"

#  List all .h5 files in the folder with full paths
h5_files <- list.files(path = folder_path, pattern = "\\.h5$", full.names = TRUE)

# Read each .h5 file and create a Seurat object for each dataset
seurat_list <- lapply(h5_files, function(file) {data <- Read10X_h5(filename = file)
  CreateSeuratObject(counts = data, project = basename(file))})

# Name each Seurat object based on file names (without .h5 extension)
names(seurat_list) <- gsub("\\.h5$", "", basename(h5_files))

#Quality Control (QC) for each Seurat object
seurat_list_qc <- lapply(seurat_list, function(seurat_obj) {
  # Calculate percentage of mitochondrial genes
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Filter cells based on QC thresholds:
  # Keep cells with > 200 and < 2500 detected features, and < 5% mitochondrial genes
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  return(seurat_obj)
})


# Create PDF with QC violin plots for all samples
pdf("QC_violin_plots.pdf")

for (i in seq_along(seurat_list_qc)) {
  obj <- seurat_list_qc[[i]]
  name <- names(seurat_list_qc)[i]
  
  # Violin plots of nFeature_RNA, nCount_RNA, percent.mt
  pnl <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 ncol = 3, pt.size = 0.1, combine = TRUE) +
    patchwork::plot_annotation(title = name)
  print(pnl)
}

dev.off()

#Create PDF with FeatureScatter plots for QC metrics
pdf("QC_FeatureScatter_plots.pdf")

for (i in seq_along(seurat_list_qc)) {
  obj <- seurat_list_qc[[i]]
  name <- names(seurat_list_qc)[i]
  
  # Scatter plot: total counts vs % mitochondrial genes
  p1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    ggtitle(paste(name, "- nCount vs percent.mt"))
  
  # Scatter plot: total counts vs number of features
  p2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    ggtitle(paste(name, "- nCount vs nFeature"))
  
  print(p1)
  print(p2)
}

dev.off()

# Normalize data for each Seurat object
seurat_list_qc <- lapply(seurat_list_qc, function(seurat_obj) {
  seurat_obj <- NormalizeData(seurat_obj)
  return(seurat_obj)
})

# Identify variable features (genes) for each object
seurat_list_qc <- lapply(seurat_list_qc, function(obj) {
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
  return(obj)
})

# Show top 20 variable genes from the first sample
head(VariableFeatures(seurat_list_qc[[1]]), 20)  # Muestra los primeros 20 genes variables
VariableFeaturePlot(seurat_list_qc[[1]])


# Select integration features across all samples
features <- SelectIntegrationFeatures(object.list = seurat_list_qc, nfeatures = 3000)

# Scale data and run PCA on each object using the selected features
seurat_list_qc <- lapply(seurat_list_qc, function(obj) {
  obj <- ScaleData(obj, features = features, verbose = FALSE)
  obj <- RunPCA(obj, features = features, verbose = FALSE)
  return(obj)
})

#PROBLEM!!!!: check cells and genes for each seurat obj
for (sample_name in names(seurat_list_qc)) {
  obj <- seurat_list_qc[[sample_name]]
  cat("Sample:", sample_name, "\n")
  cat("# genes:", nrow(obj), "\n")
  cat("# cells:", ncol(obj), "\n\n")
}


# Find integration anchors (to align datasets)
anchors <- FindIntegrationAnchors(object.list = seurat_list_qc, anchor.features = features)

# Integrate data into a single Seurat object
seurat_integrated <- IntegrateData(anchorset = anchors)


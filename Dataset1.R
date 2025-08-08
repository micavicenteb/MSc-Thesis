# Install packages (run once)
install.packages("hdf5r")
install.packages("sctransform") #sctransform improves speed and memory consumption
    #this package improves the speed of the learning procedure.
    install.packages("BiocManager")
    BiocManager::install("glmGamPoi")

# Load necessary libraries
library(Seurat)
library(patchwork) # Needed for plot_annotation()
library(dplyr)
library(ggplot2)
library(hdf5r)
library(glmGamPoi)
library(sctransform)
    object <- SCTransform(object, vst.flavor = "v2") # invoke sctransform - requires Seurat>=4.1
    

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

# Remove datasets that are too small BEFORE processing
seurat_list_qc <- seurat_list_qc[sapply(seurat_list_qc, ncol) > 50 & sapply(seurat_list_qc, nrow) > 1000]


# Procesamiento completo por objeto Seurat
seurat_list_qc <- lapply(seurat_list_qc, function(obj) {
  obj <- SCTransform(obj, vst.flavor = "v2", verbose = FALSE) %>%
    RunPCA(npcs = 20, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:20, n.neighbors = 30, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
    FindClusters(resolution = 0.5, verbose = FALSE)
  return(obj)
})

#crear carpeta
dir.create("Processed_Seurat_Objects", showWarnings = FALSE)
# Guardar los objetos procesados con SCTransform
for (name in names(seurat_list_qc)) {
  saveRDS(seurat_list_qc[[name]], file = file.path("Processed_Seurat_Objects", paste0(name, ".rds")))
}

#objeto para integracion
seurat_list_sct <- lapply(seurat_list_qc, function(obj) {SCTransform(obj, vst.flavor = "v2", verbose = FALSE)})

#INTEGRATION
#Select integration features using only SCT-normalized objects (seurat_list_sct)
features <- SelectIntegrationFeatures(object.list = seurat_list_sct, nfeatures = 3000)

# Prepare SCT objects for integration
seurat_list_sct <- PrepSCTIntegration(object.list = seurat_list_sct, anchor.features = features)

#Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat_list_sct, normalization.method = "SCT", anchor.features = features)

#Integrate data
seurat_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# Save the Seurat object as an RDS file
if (!dir.exists("output")) {dir.create("output")}
saveRDS(seurat_integrated, file = "output/seurat_integrated.rds")

#----------------------CODES I USED TO FIND THE PROBLEM-----------------
# Confirm all objects have the same anchor features
all(sapply(seurat_list_sct, function(obj) all(features %in% rownames(obj))))

# Check cell counts for each object
sapply(seurat_list_sct, ncol)

# Check feature counts for each object
sapply(seurat_list_sct, nrow)

# Print cell and gene counts per sample
for (sample_name in names(seurat_list_sct)) {
  obj <- seurat_list_sct[[sample_name]]
  cat("Sample:", sample_name, "\n")
  cat("# genes:", nrow(obj), "\n")
  cat("# cells:", ncol(obj), "\n\n")}

#use dims ≤ number of PCs possible for thesmallest dataset
max_dims <- min(30, min(sapply(seurat_list_sct, function(x) min(50, ncol(x) - 1))))

#########
# PANDO #
#########

# Load Packages
library(Pando)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(qs)

# Define dataset  
dataset <- "timecourse" # timecourse | combined

# Set data directory
data.dir <- paste0("/g/scb/zaugg/deuner/GRaNIE/tmp/", dataset, ".pp.seuratObject.qs")

# Load seurat object
seurat_object <- qread

# Get motif data
data(motifs)

# Select variable features
seurat_object <- Seurat::FindVariableFeatures(seurat_object, assay='RNA')

# Initiate GRN object and select candidate regions
seurat_object <- initiate_grn(seurat_object)

# Scan candidate regions for TF binding motifs
seurat_object <- find_motifs(
  seurat_object,
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# Infer gene regulatory network
seurat_object <- infer_grn(seurat_object)

# Print inferred coefficients
coef(seurat_object)

# Find gene and regulatory modules 
test_srt <- find_modules(test_srt)

# Print modules
NetworkModules(test_srt)

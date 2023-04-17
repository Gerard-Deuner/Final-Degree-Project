#########
# PANDO #
#########

# https://quadbiolab.github.io/Pando/

# Load Packages
library(Pando)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(qs)
library(Signac)

# Define dataset  
dataset <- "timecourse" # timecourse | combined

# Set data directory
data.dir <- paste0("/g/scb/zaugg/deuner/GRaNIE/tmp/", dataset, ".pp.seuratObject.qs")

# Load seurat object
seurat_object <- qread(data.dir)

# Get motif data
data(motifs)

# Add genome annotation for ChromatinAssay
DefaultAssay(seurat_object) <- "ATAC"
genome(seurat_object) <- "hg38"

# Add annotation
## get peak ranges (timecourse)
peakRanges <- qread("/g/scb/zaugg/marttine/dariaMultiome/peakRanges.timecourse.qs")
Annotation(seurat_object[["ATAC"]]) <- peakRanges

# Select variable features
seurat_object <- Seurat::FindVariableFeatures(seurat_object, assay='RNA')

# Initiate GRN object and select candidate regions
seurat_object <- initiate_grn(seurat_object, 
                              peak_assay = "ATAC",
                              rna_assay = "RNA",
                              exclude_exons = T)

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
modules <- NetworkModules(test_srt)

# Save Network output
write.csv(modules@meta, paste0("/g/scb/zaugg/deuner/Pando/outputdata/", dataset, ".GRN.tsv"))

##############################################
# Seurat Object conversion to Anndata Object #
##############################################

# reference: https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html

# load libraries
library(qs)
library(Seurat)
library(SeuratData)
library(SeuratDisk)

# define dataset I want to convert
dataset <- "timecourse"

# set data directory
data.dir <- "/home/deuner/deuner/GRaNIE/tmp/"

# load the corresponding seurat object
s.obj <- qread(paste0(data.dir, dataset, ".pp.seuratObject.qs"))

# conversion
SaveH5Seurat(s.obj, filename = paste0(data.dir, dataset, ".h5Seurat"))
Convert(paste0(data.dir, dataset, ".h5Seurat"), dest = "h5ad")

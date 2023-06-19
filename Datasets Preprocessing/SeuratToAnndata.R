##############################################
# Seurat Object conversion to Anndata Object #
##############################################

# reference: https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html

# load libraries
library(qs)
library(Seurat)
library(SeuratData)
library(SeuratDisk)

# read dataset argument from command line
args <- commandArgs(trailingOnly = TRUE)

# choose dataset
dataset <- args[1]

# set data directory
data.dir <- "/home/deuner/deuner/GRaNIE/tmp/"

# load the corresponding seurat object
s.obj <- qread(paste0(data.dir, dataset, ".pp.nomicro.seuratObject.qs"))
DefaultAssay(s.obj) <- "RNA"

# convert factor columns to character
i <- sapply(s.obj@meta.data, is.factor)
s.obj@meta.data[i] <- lapply(s.obj@meta.data[i], as.character)

# conversion
SaveH5Seurat(s.obj, filename = paste0(data.dir, dataset, ".nomicro.h5Seurat"), overwrite = T)
Convert(paste0(data.dir, dataset, ".nomicro.h5Seurat"), dest = "h5ad", overwrite = T)

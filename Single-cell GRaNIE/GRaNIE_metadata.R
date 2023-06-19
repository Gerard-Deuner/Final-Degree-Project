############################################################
# Run GRaNIE for a specific metadata column. e.g. celltype #
############################################################

# Load libraries
library(GRaNIE)
library(Seurat)
library(qs)
library(dplyr)
library(Signac)

# read arguments from command line
args <- commandArgs(trailingOnly = TRUE)

# choose dataset
dataset <- args[1] # timecourse | combined

# choose correlation method
correlation.method <- args[2] # pearson | spearman

# choose metadata column  
meta <- args[3] # celltype, celltype_wnn, etc

# Set up source of helper functions 
source("/g/scb/zaugg/deuner/GRaNIE/code/GRaNIE_helper_functions.R")

# Set genome assembly version
genomeAssembly = "hg38"

# Set up the main directory
path = "/g/scb/zaugg/deuner/GRaNIE"

# Use Zaugg internal TFBS folder
TFBS_folder = NULL

# Load feature file that gives Ensembl IDs and gene names to translate names to Ensembl IDs. 
file_RNA_features = paste0("/g/zaugg/carnold/Projects/GRN_pipeline/misc/singleCell/sharedMetadata/features_RNA_", genomeAssembly, ".tsv.gz")

# Load the seurat object
s.obj <- qread(paste0(path,"/tmp/", dataset, ".pp.seuratObject.qs"))

# Path to the output directory
seurat_outputFolder = paste0(path,"/outputdata/", dataset, "_", correlation.method,"_", meta)

# Prepare data
s.obj = prepareSeuratData_GRaNIE(s.obj, outputDir = seurat_outputFolder, pseudobulk_source = meta,
                                 file_RNA_features = file_RNA_features,
                                 clusterResolutions = res,
                                 prepareData = TRUE,
                                 saveSeuratObject = TRUE)

# runGRaNIE for that resolution
GRN = runGRaNIE(
  datasetName = "Timecourse_Dataset",
  dir_output = paste0(seurat_outputFolder,"/output_pseudobulk_", meta, "_RNA_limma_quantile_ATAC_DESeq2_sizeFactors"), 
  file_peaks = paste0(seurat_outputFolder,"/atac.pseudobulkFromClusters_", meta, ".tsv.gz"), 
  file_rna = paste0(seurat_outputFolder,"/rna.pseudobulkFromClusters_", meta, ".tsv.gz"), 
  file_metadata = paste0(seurat_outputFolder,"/metadata_", meta, ".tsv.gz"),
  genomeAssembly = "hg38", 
  nCores = 8, 
  runNetworkAnalyses = FALSE,
  correlation.method = correlation.method) 


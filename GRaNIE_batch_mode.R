#############################################################################################
# GRaNIE Batch Mode using SPEARMAN or PEARSON Correlation on TIMECOURSE or COMBINED dataset #
#############################################################################################

# Load libraries
library(GRaNIE)
library(Seurat)
library(qs)
library(dplyr)
library(Signac)

# read arguments from command line
args <- commandArgs(trailingOnly = TRUE)

# choose dataset
dataset <- args[1]

# choose correlation method
correlation.method <- args[2]

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
seurat_outputFolder = paste0(path,"/outputdata/batch_mode/", dataset,"_batch_mode_", correlation.method)

# Prepare data
s.obj = prepareSeuratData_GRaNIE(s.obj, outputDir = seurat_outputFolder, pseudobulk_source = "cluster", 
                                      assayName_RNA = "RNA", assayName_ATAC= "ATAC",
                                      file_RNA_features = file_RNA_features,
                                      prepareData = TRUE,
                                      saveSeuratObject = TRUE)

# runGRaNIE in batch mode
runGRaNIE_batchMode(datasetName = paste(dataset, "Dataset", sep = "_"),
                    inputDir = seurat_outputFolder, 
                    outputDir = paste0(seurat_outputFolder, "/Batch_Mode_Outputs"),
                    genomeAssembly = genomeAssembly,
                    TFBS_folder = TFBS_folder,
                    nCores = 8,
                    TF_peak.fdr.threshold = 0.2,
                    peak_gene.fdr.threshold = 0.1,
                    runNetworkAnalyses = FALSE,
                    forceRerun = TRUE,
                    correlation.method = correlation.method) 

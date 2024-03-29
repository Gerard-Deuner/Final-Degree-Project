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
dataset <- args[1] # timecourse | combined

# choose correlation method
corr.method <- args[2] # pearson | spearman

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
s.obj <- qread(paste0(path,"/tmp/", dataset, ".pp.nomicro.seuratObject.qs"))

# Path to the output directory
seurat_outputFolder = paste0(path,"/outputdata/batch_mode/", dataset,"_batch_mode_", corr.method, "_nomicro")
dir.create(paste0(path,"/outputdata/batch_mode/", dataset,"_batch_mode_", corr.method, "_nomicro"))

# Prepare data
s.obj = prepareSeuratData_GRaNIE(s.obj, outputDir = seurat_outputFolder, pseudobulk_source = "cluster",
                                 file_RNA_features = file_RNA_features,
                                 prepareData = TRUE,
                                 saveSeuratObject = TRUE)

# create outputs directory
dir.create(paste0(seurat_outputFolder, "/Batch_Mode_Outputs"))
 
# runGRaNIE in batch mode
runGRaNIE_batchMode(datasetName = paste(dataset, "Dataset", sep = "_"),
                    inputDir = seurat_outputFolder,
                    outputDir = paste0(seurat_outputFolder, "/Batch_Mode_Outputs"),
                    clusterResolutions = c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2)),
                    genomeAssembly = genomeAssembly,
                    TFBS_folder = TFBS_folder,
                    nCores = 8,
                    TF_peak.fdr.threshold = 0.2,
                    peak_gene.fdr.threshold = 0.1,
                    runNetworkAnalyses = TRUE,
                    forceRerun = TRUE,
                    correlation.method = corr.method)


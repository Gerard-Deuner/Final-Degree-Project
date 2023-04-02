################################
# scGRaNIE ON COMBINED DATASET #
################################

######################
# SET UP ENVIRONMENT #
######################

# Load libraries
library(GRaNIE)
library(Seurat)
library(qs)
library(dplyr)


########################################################
# STANDARD RUN - CELL TYPE PSEUDOBULKING - 14 CLUSTERS #
########################################################

# Set up source of helper functions
source("/g/zaugg/carnold/Projects/GRN_pipeline/misc/singleCell/functions_singleCell_GRaNIE.R")

# Set genome assembly version
genomeAssembly = "hg38"

# Set up the main directory
path = "/g/scb/zaugg/deuner/GRaNIE"

# Use Zaugg internal TFBS folder
TFBS_folder = NULL

# Load feature file that gives Ensembl IDs and gene names to translate names to Ensembl IDs. 
file_RNA_features = paste0("/g/zaugg/carnold/Projects/GRN_pipeline/misc/singleCell/sharedMetadata/features_RNA_", genomeAssembly, ".tsv.gz")

# Load the seurat object
seuratFile = "timecourse.pp.seuratObject.qs"
timecourse.s <- qread(paste0(path,"/tmp/",seuratFile))

# Path to output directory
seurat_outputFolder = paste0(path,"/outputdata/timecourse_celltype")

# Prepare data
timecourse = prepareSeuratData_GRaNIE(timecourse.s, outputDir = seurat_outputFolder, pseudobulk_source = "celltype", 
                                      assayName_RNA = "SCT", assayName_ATAC = "ATAC",
                                      file_RNA_features = file_RNA_features,
                                      saveSeuratObject = TRUE,
                                      prepareData = FALSE)

# Run GRaNIE
GRN = runGRaNIE(
  datasetName = "Timecourse_Dataset",
  dir_output = paste0(seurat_outputFolder,"/output_pseudobulk_celltype_RNA_limma_quantile_ATAC_DESeq2_sizeFactors"), 
  file_peaks = paste0(seurat_outputFolder,"/atac.pseudobulkFromClusters_celltype.tsv.gz"), 
  file_rna = paste0(seurat_outputFolder,"/rna.pseudobulkFromClusters_celltype.tsv.gz"), 
  file_metadata = paste0(seurat_outputFolder,"/metadata_celltype.tsv.gz"),
  genomeAssembly = "hg38", 
  nCores = 8, 
  runNetworkAnalyses = FALSE) 


#########################################################
# BATCH MODE - RUN IT FOR DIFFERENT CLUSTER RESOLUTIONS # (on SWNN graph)
#########################################################

# Set up source of helper functions (spearman correlation)
source("/g/zaugg/deuner/GRaNIE/code/GRaNIE_spearman_helper_functions")

# Path to the output directory
seurat_outputFolder = paste0(path,"/outputdata/timecourse_batch_mode")

# Prepare data
timecourse = prepareSeuratData_GRaNIE(timecourse.s, outputDir = seurat_outputFolder, pseudobulk_source = "cluster", 
                                      assayName_RNA = "SCT", assayName_ATAC= "ATAC",
                                      file_RNA_features = file_RNA_features,
                                      prepareData = TRUE,
                                      saveSeuratObject = TRUE)

# create outputs directory
dir.create("/g/scb/zaugg/deuner/GRaNIE/outputdata/timecourse_batch_mode/Batch_Mode_Outputs")

# Run GRaNIE
runGRaNIE_batchMode(datasetName = "Timecourse_Dataset",
                    inputDir = seurat_outputFolder, 
                    outputDir = paste0(seurat_outputFolder, "/Batch_Mode_Outputs"),
                    genomeAssembly = genomeAssembly,
                    TFBS_folder = TFBS_folder,
                    nCores = 8,
                    TF_peak.fdr.threshold = 0.2,
                    peak_gene.fdr.threshold = 0.1,
                    runNetworkAnalyses = FALSE,
                    forceRerun = TRUE) 


##############################################################################
# STANDARD RUN - CELL TYPE PSEUDOBULKING - REMOVE TINY CLUSTER - 13 CLUSTERS #
##############################################################################

# Path to the output directory
seurat_outputFolder = paste0(path,"/outputdata/timecourse_celltype_notinycluster")

# Load the seurat object
timecourse.s.2 <- qread(paste0(path,"/tmp/",seuratFile))

# Filter out cells belonging to the tiny cluster
timecourse.s.2 <- subset(timecourse.s.2, subset = celltype != "diff - tiny")
levels(timecourse.s.2@meta.data$celltype)
WhichCells(timecourse.s.2, idents = "diff - tiny")
Idents(timecourse.s.2)

# Prepare data
timecourse = prepareSeuratData_GRaNIE(timecourse.s.2, outputDir = seurat_outputFolder, pseudobulk_source = "celltype", 
                                      assayName_RNA = "SCT", assayName_ATAC= "ATAC",
                                      file_RNA_features = file_RNA_features,
                                      saveSeuratObject = TRUE,
                                      prepareData = FALSE)

# Run GRaNIE
GRN = runGRaNIE(
  datasetName = "Timecourse_Dataset",
  dir_output = paste0(seurat_outputFolder,"/output_pseudobulk_celltype_RNA_limma_quantile_ATAC_DESeq2_sizeFactors"), 
  file_peaks = paste0(seurat_outputFolder,"/atac.pseudobulkFromClusters_celltype.tsv.gz"), 
  file_rna = paste0(seurat_outputFolder,"/rna.pseudobulkFromClusters_celltype.tsv.gz"), 
  file_metadata = paste0(seurat_outputFolder,"/metadata_celltype.tsv.gz"),
  genomeAssembly = "hg38", 
  nCores = 8, 
  runNetworkAnalyses = FALSE) 


#################################################################################
# STANDARD RUN - CELL TYPE PSEUDOBULKING - REMOVE 2 WEIRD CLUSTER - 12 CLUSTERS #
#################################################################################

# Path to the output directory
seurat_outputFolder = paste0(path,"/outputdata/timecourse_celltype_clear_clusters")

# Load the seurat object
timecourse.s.3 <- qread(paste0(path,"/tmp/",seuratFile))

# Filter out cells belonging to the tiny cluster
timecourse.s.3 <- subset(timecourse.s.2, subset = celltype %in% c("diff - tiny", "diff - RUNX1"))
timecourse.s.3 <- subset(timecourse.s.2, subset = celltype != "diff - tiny")
timecourse.s.3 <- subset(timecourse.s.2, subset = celltype != "diff - RUNX1")
levels(timecourse.s.3@meta.data$celltype)
WhichCells(timecourse.s.2, idents = "diff - tiny")
WhichCells(timecourse.s.2, idents = "diff - RUNX1")
Idents(timecourse.s.2)

# Prepare data
timecourse = prepareSeuratData_GRaNIE(timecourse.s.3, outputDir = seurat_outputFolder, pseudobulk_source = "celltype", 
                                      assayName_RNA = "SCT", assayName_ATAC= "ATAC",
                                      file_RNA_features = file_RNA_features,
                                      saveSeuratObject = TRUE,
                                      prepareData = FALSE)

# Run GRaNIE
GRN = runGRaNIE(
  datasetName = "Timecourse_Dataset",
  dir_output = paste0(seurat_outputFolder,"/output_pseudobulk_celltype_RNA_limma_quantile_ATAC_DESeq2_sizeFactors"), 
  file_peaks = paste0(seurat_outputFolder,"/atac.pseudobulkFromClusters_celltype.tsv.gz"), 
  file_rna = paste0(seurat_outputFolder,"/rna.pseudobulkFromClusters_celltype.tsv.gz"), 
  file_metadata = paste0(seurat_outputFolder,"/metadata_celltype.tsv.gz"),
  genomeAssembly = "hg38", 
  nCores = 8, 
  runNetworkAnalyses = FALSE) 


#####################################################
# CELLTYPE CLUSTERING ON WNN CLUSTERS - 16 CLUSTERS #
#####################################################

# Set up source of helper functions
source("/g/zaugg/carnold/Projects/GRN_pipeline/misc/singleCell/functions_singleCell_GRaNIE.R")

# Set genome assembly version
genomeAssembly = "hg38"

# Set up the main directory
path = "/g/scb/zaugg/deuner/GRaNIE"

# Use Zaugg internal TFBS folder
TFBS_folder = NULL

# Load feature file that gives Ensembl IDs and gene names to translate names to Ensembl IDs. 
file_RNA_features = paste0("/g/zaugg/carnold/Projects/GRN_pipeline/misc/singleCell/sharedMetadata/features_RNA_", genomeAssembly, ".tsv.gz")

# Load the seurat object
seuratFile = "timecourse.pp.seuratObject.qs"
timecourse.s <- qread(paste0(path,"/tmp/",seuratFile))

# Path to output directory
seurat_outputFolder = paste0(path,"/outputdata/timecourse_celltype_wnn")

# Prepare data
timecourse = prepareSeuratData_GRaNIE(timecourse.s, outputDir = seurat_outputFolder, pseudobulk_source = "celltype_wnn", 
                                      assayName_RNA = "SCT", assayName_ATAC = "ATAC",
                                      file_RNA_features = file_RNA_features,
                                      saveSeuratObject = TRUE,
                                      prepareData = FALSE)

# Run GRaNIE
GRN = runGRaNIE(
  datasetName = "Timecourse_Dataset",
  dir_output = paste0(seurat_outputFolder,"/output_pseudobulk_celltype_wnn_RNA_limma_quantile_ATAC_DESeq2_sizeFactors"), 
  file_peaks = paste0(seurat_outputFolder,"/atac.pseudobulkFromClusters_celltype_wnn.tsv.gz"), 
  file_rna = paste0(seurat_outputFolder,"/rna.pseudobulkFromClusters_celltype_wnn.tsv.gz"), 
  file_metadata = paste0(seurat_outputFolder,"/metadata_celltype_wnn.tsv.gz"),
  genomeAssembly = "hg38", 
  nCores = 8, 
  runNetworkAnalyses = FALSE) 


######################################################
# SPEARMAN CORRELATION CELL TYPE PSEUDOBULKING (WNN) #
######################################################

# Set up source of helper functions
source("/g/zaugg/carnold/Projects/GRN_pipeline/misc/singleCell/functions_singleCell_GRaNIE.R")

# Set genome assembly version
genomeAssembly = "hg38"

# Set up the main directory
path = "/g/scb/zaugg/deuner/GRaNIE"

# Use Zaugg internal TFBS folder
TFBS_folder = NULL

# Load feature file that gives Ensembl IDs and gene names to translate names to Ensembl IDs. 
file_RNA_features = paste0("/g/zaugg/carnold/Projects/GRN_pipeline/misc/singleCell/sharedMetadata/features_RNA_", genomeAssembly, ".tsv.gz")

# Load the seurat object
seuratFile = "timecourse.pp.seuratObject.qs"
timecourse.s <- qread(paste0(path,"/tmp/",seuratFile))

# Path to output directory
seurat_outputFolder = paste0(path,"/outputdata/timecourse_celltype_spearman_wnn")

# Prepare data
timecourse = prepareSeuratData_GRaNIE(timecourse.s, outputDir = seurat_outputFolder, pseudobulk_source = "celltype_wnn", 
                                      assayName_RNA = "SCT", assayName_ATAC = "ATAC",
                                      file_RNA_features = file_RNA_features,
                                      saveSeuratObject = TRUE,
                                      prepareData = FALSE)

# Run GRaNIE
GRN = runGRaNIE.spearman(
  datasetName = "Timecourse_Dataset",
  dir_output = paste0(seurat_outputFolder,"/output_pseudobulk_celltype_wnn_RNA_limma_quantile_ATAC_DESeq2_sizeFactors"), 
  file_peaks = paste0(seurat_outputFolder,"/atac.pseudobulkFromClusters_celltype_wnn.tsv.gz"), 
  file_rna = paste0(seurat_outputFolder,"/rna.pseudobulkFromClusters_celltype_wnn.tsv.gz"), 
  file_metadata = paste0(seurat_outputFolder,"/metadata_celltype_wnn.tsv.gz"),
  genomeAssembly = "hg38", 
  nCores = 8, 
  runNetworkAnalyses = FALSE) 


# Adapt runGRaNIE function
runGRaNIE.spearman <- function(dir_output = "output_GRaNIE", 
                      datasetName = "undescribed",
                      file_peaks, file_rna, file_metadata,
                      TFBS_folder = NULL,
                      genomeAssembly = "hg38",
                      normalization_peaks = "DESeq2_sizeFactors", 
                      idColumn_peaks = "peakID",
                      normalization_rna = "limma_quantile", 
                      idColumn_RNA =  "ENSEMBL",
                      includeSexChr = FALSE,
                      minCV = 0,
                      minNormalizedMean_peaks = 5,
                      minNormalizedMean_RNA = 1,
                      minSizePeaks = 5,
                      promoterRange = 250000, 
                      useGCCorrection = FALSE,
                      TF_peak.fdr.threshold = 0.2,
                      peak_gene.fdr.threshold = 0.1,
                      runNetworkAnalyses = FALSE, 
                      nCores = 8,
                      forceRerun = TRUE
) {
  
  library(tidyverse)
  library(GRaNIE)
  
  file_GRN = paste0(dir_output, "/GRN.qs")
  
  if (file.exists(file_GRN) & !forceRerun) {
    cat(" Skip, output file already present...")
    next
  }
  
  
  countsATAC   = read_tsv(file_peaks)
  countsRNA    = read_tsv(file_rna)
  metadata.all = read_tsv(file_metadata) 
  
  
  # Arbitrary list with information and metadata that is stored within the GRN object
  metadata.l = list(name = datasetName,
                    file_peaks = file_peaks,
                    file_rna  = file_rna,
                    genomeAsembly = genomeAssembly,
                    file_metadata = file_metadata 
  )
  
  
  GRN = initializeGRN(objectMetadata = metadata.l, 
                      outputFolder = dir_output,
                      genomeAssembly = genomeAssembly)
  
  
  GRN = addData(GRN,
                counts_peaks = countsATAC, normalization_peaks = normalization_peaks, idColumn_peaks = idColumn_peaks,
                counts_rna = countsRNA, normalization_rna = normalization_rna, idColumn_RNA = idColumn_RNA,
                sampleMetadata = metadata.all, allowOverlappingPeaks = TRUE, 
                forceRerun = forceRerun)
  
  GRN = plotPCA_all(GRN, data = c("rna", "peaks"), topn = 500, type = "normalized", removeFiltered = FALSE, forceRerun = forceRerun)
  
  # Should the pipeline be run for only a subset of TFs or all? The special keyword "all" will use all TF that are found in the HOCOMOCO folder; however, if only a subset should be considered, specify the subset here with c() and the TF names, as shown below
  TFs = "all"
  
  if (is.null(TFBS_folder)) {
    # Base directory of the folder with the TFBS predictions. 
    # The TFBS predictions are expected as *.bed files as well as a translation table with the name translationTable.csv
    # We provide all files here: https://www.embl.de/download/zaugg/GRN/hg19_hg38_mm10_PWMScan.zip (7.5 GB)
    # Make sure they are in the same genome assembly as the peaks data
    hocomocoVersion = if_else(genomeAssembly == "hg38", "v11", "v10")
    motifFolder = paste0("/g/zaugg/zaugg_shared/annotations/TFBS/", genomeAssembly, "/PWMScan_HOCOMOCO", hocomocoVersion)
    
  } else {
    motifFolder = TFBS_folder
  }
  
  
  GRN = addTFBS(GRN, motifFolder = motifFolder, TFs = TFs, filesTFBSPattern = "_TFBS", fileEnding = ".bed", forceRerun = forceRerun)
  
  ######################
  # PEAKS TFBS OVERLAP #
  ######################
  
  GRN = overlapPeaksAndTFBS(GRN, nCores = nCores, forceRerun = forceRerun)
  
  qs::qsave(GRN, file_GRN)
  
  
  # Chromosomes to include for peaks and peak-gene associations. This should be a vector of chromosome names
  if (includeSexChr) {
    if (stringr::str_starts(genomeAssembly, "mm")) {
      chrToKeep_peaks = c(paste0("chr", 1:19), "chrX", "chrY")
    } else {
      chrToKeep_peaks = c(paste0("chr", 1:22), "chrX", "chrY")
    }
    
  } else {
    if (stringr::str_starts(genomeAssembly, "mm")) {
      chrToKeep_peaks = c(paste0("chr", 1:19))
    } else {
      chrToKeep_peaks = c(paste0("chr", 1:22))
    }
  }
  
  GRN = filterData(GRN, minNormalizedMean_peaks = minNormalizedMean_peaks, minNormalizedMeanRNA = minNormalizedMean_RNA, 
                   chrToKeep_peaks = chrToKeep_peaks, minSize_peaks = minSizePeaks,
                   minCV_peaks = minCV, minCV_genes = minCV, forceRerun = forceRerun)
  
  GRN = addConnections_TF_peak(GRN, connectionTypes = c("expression"), plotDiagnosticPlots = FALSE, plotDetails = FALSE, 
                               corMethod = "spearman", maxFDRToStore = 0.3,
                               useGCCorrection = useGCCorrection, percBackground_size = 75, percBackground_resample = TRUE,
                               forceRerun = forceRerun)
  
  
  file_input_TADs = ""
  overlapTypeGene = "TSS"
  
  GRN = addConnections_peak_gene(GRN,
                                 overlapTypeGene = overlapTypeGene,
                                 corMethod = "spearman", shuffleRNACounts = TRUE,
                                 promoterRange = promoterRange, TADs = NULL,
                                 nCores = nCores, plotDiagnosticPlots = TRUE,
                                 forceRerun = forceRerun)
  
  
  GRN = filterGRNAndConnectGenes(GRN, TF_peak.fdr.threshold = TF_peak.fdr.threshold, 
                                 TF_peak.connectionTypes = "expression" ,
                                 peak_gene.fdr.threshold = peak_gene.fdr.threshold,
                                 gene.types = c("protein_coding"),
                                 allowMissingTFs = FALSE, allowMissingGenes = FALSE,
                                 peak_gene.r_range = c(0,1))
  
  file_connections = paste0(dir_output, "/connections_TFPeak", TF_peak.fdr.threshold, "_peakGene", peak_gene.fdr.threshold, ".tsv.gz")
  
  GRN = add_TF_gene_correlation(GRN, nCores = nCores)
  connections.df = getGRNConnections(GRN, 
                                     include_TF_gene_correlations = TRUE, 
                                     include_peakMetadata = TRUE, 
                                     include_TFMetadata = TRUE, 
                                     include_geneMetadata = TRUE)
  
  write_tsv(connections.df, file_connections)
  
  if (runNetworkAnalyses) {
    GRN = performAllNetworkAnalyses(GRN)
  }
  
  
  # GRN = visualizeGRN(GRN, plotAsPDF = FALSE)
  
  qs::qsave(GRN, file_GRN)
  
  GRN
  
}


#########################################################
# RUN GRaNIE FOR MISSING CLUSTER RESOLUTIONS (SPEARMAN) #
#########################################################

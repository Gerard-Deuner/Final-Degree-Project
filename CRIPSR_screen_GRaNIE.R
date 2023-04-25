###################################
# scGRaNIE on NPC and Neuron data # 
###################################

#Load libraries
library(qs)
library(Seurat)
library(Matrix)

#Set directory and cell type
dir <- "/g/scb/zaugg/claringb/scz_CRISPR_screen/NGN2_NPC/scifi-RNA/output/processing_by_Mikael_Umut/"

## Combine seurat objects
#Read Seurat objects
npc <- qread(paste0(dir, "/6.Aligned/NPC/NPCseuratObject.filtered.SCTransformed.clustered.co-embedding.qs"))
neur <- qread(paste0(dir, "/6.Aligned/Neuron/neuronseuratObject.filtered.SCTransformed.clustered.co-embedding.qs"))

#Clear metadata for a cleaner output
columns.to.remove <- c("leidenSub0.4", "leidenSub0.6", "leidenSub0.8", "leidenSub1.0",
                       "leidenSub1.4", "leidenSub1.6", "umap1", "umap2")

for(i in columns.to.remove) {
  npc[[i]] <- NULL
  neur[[i]] <- NULL
}

#Combine into one seurat object
seur.rna <- merge(npc, y = neur, add.cell.ids = c("NPC", "Neuron"))

#Add ATAC data
peaks <- readMM("/g/scb/zaugg/marttine/RNA_ATAC_integration/CRISPR/input/atac.counts.mtx.gz")
rownames(peaks) <- read.csv("/g/scb/zaugg/marttine/RNA_ATAC_integration/CRISPR/input/peaks.atac.csv.gz")$x
colnames(peaks) <- read.csv("/g/scb/zaugg/marttine/RNA_ATAC_integration/CRISPR/input/cell.atac.csv.gz")$id

## create a seurat object
seur.atac <- CreateSeuratObject(counts = peaks, assay = 'ATAC', project = 'CRISPR')

#Combine rna and atac seurat object
seur <- merge(seur.rna, seur.atac)

#Save seurat object
qsave(seur, "/g/scb/zaugg/deuner/GRaNIE/tmp/npc_neuron.seuratObject.qs")

#Read seurat object
seur <- qread("/g/scb/zaugg/deuner/GRaNIE/tmp/npc_neuron.seuratObject.qs")

## Run scGRaNIE
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


# Path to the output directory
seurat_outputFolder = paste0(path,"/outputdata/", "NPC_Neuron_leiden1.2")

# Prepare data
seur = prepareSeuratData_GRaNIE(seur, outputDir = seurat_outputFolder, pseudobulk_source = "leidenSub1.2",
                                 file_RNA_features = file_RNA_features,
                                 prepareData = TRUE,
                                 saveSeuratObject = TRUE)

# runGRaNIE for the specified metadata column
GRN = runGRaNIE(
  datasetName = "CRISPR_NPC_Neuron_Dataset",
  dir_output = paste0(seurat_outputFolder,"/output_pseudobulk_leidenSub1.2_RNA_limma_quantile_ATAC_DESeq2_sizeFactors"), 
  file_peaks = paste0(seurat_outputFolder,"/atac.pseudobulkFromClusters_leidenSub1.2.tsv.gz"), 
  file_rna = paste0(seurat_outputFolder,"/rna.pseudobulkFromClusters_leidenSub1.2.tsv.gz"), 
  file_metadata = paste0(seurat_outputFolder,"/metadata_leidenSub1.2.tsv.gz"),
  genomeAssembly = "hg38", 
  nCores = 8, 
  runNetworkAnalyses = TRUE,
  correlation.method = "spearman") 

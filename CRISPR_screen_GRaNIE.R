###################################
# scGRaNIE on NPC and Neuron data # 
###################################

#Load libraries
library(qs)
library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(GRaNIE)
library(Signac)
library(data.table)
library(stringr)

# load cluster embedding for each cell and modalities
emb <- fread("/g/scb/zaugg/marttine/RNA_ATAC_integration/CRISPR/iterations/embeddingsAndLeidenClusters.csv")
rownames(emb) <- emb$V1
head(emb)

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

#Get ATAC data
peaks <- readMM("/g/scb/zaugg/marttine/RNA_ATAC_integration/CRISPR/input/atac.counts.mtx.gz")
rownames(peaks) <- read.csv("/g/scb/zaugg/marttine/RNA_ATAC_integration/CRISPR/input/peaks.atac.csv.gz")$x
colnames(peaks) <- read.csv("/g/scb/zaugg/marttine/RNA_ATAC_integration/CRISPR/input/cell.atac.csv.gz")$id

# Create seurat object for ATAC 
seur.atac <- CreateSeuratObject(counts = peaks, assay = 'ATAC', project = 'CRISPR')

# Add metadata to seurat objects

seur.atac <- seur.atac %>%
  AddMetaData(emb %>% dplyr::filter(modality == "ATAC"))

# adapt cell names of metadata so they match the seur.rna cell names (in seur.rna sample names are NPC and Neuron and in metadata they are lowercased)
emb <- emb %>%
  dplyr::filter(modality == "RNA") %>%
  mutate(V1 = rownames(seur.rna@meta.data))
rownames(emb) <- emb$V1

seur.rna <- seur.rna %>%
  AddMetaData(emb %>% dplyr::filter(modality == "RNA"))

( colnames(seur.rna) == emb %>% dplyr::filter(modality == "RNA") %>% pull(V1) ) %>% table()


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
seur = prepareDataSepModalities_GRaNIE(seur.rna, seur.atac, outputDir = seurat_outputFolder, pseudobulk_source = "leidenSub1.2",
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



# Prepare Data function from separate modalities
prepareDataSepModalities_GRaNIE <- function(seur.rna, seur.atac,
                                     outputDir = "pseudobulk", saveSeuratObject = TRUE,
                                     file_RNA_features = "/g/zaugg/carnold/Projects/GRN_pipeline/misc/singleCell/sharedMetadata/features_RNA_hg38.tsv.gz", 
                                     assayName_RNA = "RNA", assayName_ATAC= "ATAC", 
                                     prepareData = TRUE, SCT_nDimensions = 50,  dimensionsToIgnore_LSI_ATAC = 1,
                                     pseudobulk_source = "cluster",
                                     clusteringAlgorithm = 3, clusterResolutions = c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2)), 
                                     doDimPlots = TRUE,
                                     run_MultiK = FALSE,
                                     forceRerun = FALSE
) {
  
  start = Sys.time()
  
  # checkmate::assertClass(seu.s, "Seurat")
  # checkmate::assertSubset(c(assayName_RNA, assayName_ATAC), names(seu.s@assays))
  # checkmate::assertFlag(prepareData)
  # checkmate::assertFileExists(file_RNA_features, access = "r")
  # checkmate::assertSubset(pseudobulk_source, c("cluster", colnames(seu.s@meta.data)))
  # 
  # checkmate::assertFlag(saveSeuratObject)
  # checkmate::assertFlag(run_MultiK)
  # checkmate::assertSubset(clusteringAlgorithm, c(1:4), empty.ok = FALSE)
  # checkmate::assertIntegerish(clusteringAlgorithm, len = 1)
  # checkmate::assertFlag(doDimPlots)
  # checkmate::assertIntegerish(SCT_nDimensions, lower = 5, upper = 500)
  # checkmate::assertIntegerish(dimensionsToIgnore_LSI_ATAC, lower = 0, upper = SCT_nDimensions - 1)
  # checkmate::assertFlag(forceRerun)
  
  if (!dir.exists(outputDir)) {
    dir.create(outputDir)
  }
  
  GRaNIE:::.startLogger(paste0(outputDir, "/prepareData_Seurat_GRaNIE.R.log") , "INFO",  removeOldLog = FALSE)
  
  if (saveSeuratObject & !GRaNIE:::is.installed("qs")) {
    stop("Package qs is not installed but needed due to saveSeuratObject = TRUE")
  }
  
  if (saveSeuratObject) {
    # TODO check whether qs is installed
  }
  
  
  
  # Currently hard-coded options
  # See SCT transform
  returnOnlyVarGenes = FALSE
  
  # Needed to retrieve ENSEMBL IDs for gene names
  if (!is.null(file_RNA_features)) {
    futile.logger::flog.info(paste0("Reading and checking the RNA features file."))
    
    features = readr::read_tsv(file_RNA_features, col_names = FALSE, col_types = "ccffii")
    colnames(features) = c("ens_id", "gene", "view", "chr", "start", "end")
    
  }
  
  
  
  if (prepareData) {
    
    futile.logger::flog.info(paste0("Preparing RNA and ATAC data (RNA: SCTransform, PCA, UMPA; ATAC: TFID, FindTopFeatures, SVD, UMAP; combined: FindMultiModalNeighbors, UMAP). This may take a while."))
    
    Seurat::DefaultAssay(seu.rna) <- assayName_RNA
    
    dimToUseRNA = seq_len(SCT_nDimensions)
    
    # Please note that this matrix is non-sparse, and can therefore take up a lot of memory if stored for all genes. 
    # To save memory, we store these values only for variable genes, by setting the return.only.var.genes = TRUE by default in the SCTransform() function call.
    seu.s <- Seurat::SCTransform(seu.rna, verbose = FALSE, return.only.var.genes = returnOnlyVarGenes) %>% 
      Seurat::RunPCA() %>% 
      Seurat::RunUMAP(dims = dimToUseRNA, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
    
    
    Seurat::DefaultAssay(seu.atac) <- assayName_ATAC
    seu.s <- Signac::RunTFIDF(seu.atac)
    seu.s <- Signac::FindTopFeatures(seu.atac, min.cutoff = 'q0')
    seu.s <- Signac::RunSVD(seu.atac)
    
    
    dimToUse_ATAC = setdiff(1:SCT_nDimensions, dimensionsToIgnore_LSI_ATAC)
    
    seu.s <- Seurat::RunUMAP(seu.atac, reduction = 'lsi', dims = dimToUse_ATAC, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
    
    # # We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering
    # seu.s <- Seurat::FindMultiModalNeighbors(seu.s, reduction.list = list("pca", "lsi"), dims.list = list(dimToUseRNA, dimToUse_ATAC),
    #                                          weighted.nn.name = "weighted.nn",
    #                                          knn.graph.name = "wknn",
    #                                          snn.graph.name = "wsnn")
    # 
    # seu.s <- Seurat::RunUMAP(seu.s, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
    
    
  } # end if prepareData
  
  # Currently the only opton that works. Count aggregation
  sumCounts  = TRUE
  aggregationType = dplyr::if_else(sumCounts == TRUE, "", "_mean")
  
 
    futile.logger::flog.info(paste0("Creating pseudobulk based on the following pre-existing cell identity column: ", pseudobulk_source))
    
    file_in_rna  = paste0(outputDir, "/rna.pseudobulkFromClusters_" , pseudobulk_source, aggregationType, ".tsv.gz")
    file_in_atac = paste0(outputDir, "/atac.pseudobulkFromClusters_", pseudobulk_source, aggregationType, ".tsv.gz")
    
    if (file.exists(file_in_rna) & file.exists(file_in_atac)) {
      seu.s
    }
    
    # Normalize levels and make the names compatible
    seu.rna@meta.data[[pseudobulk_source]] = gsub(" ", "_", seu.rna@meta.data[[pseudobulk_source]])
    seu.rna@meta.data[[pseudobulk_source]] = gsub("-", "_", seu.rna@meta.data[[pseudobulk_source]])
    
    seu.atac@meta.data[[pseudobulk_source]] = gsub(" ", "_", seu.atac@meta.data[[pseudobulk_source]])
    seu.atac@meta.data[[pseudobulk_source]] = gsub("-", "_", seu.atac@meta.data[[pseudobulk_source]])
    
    Seurat::Idents(seu.rna) =  pseudobulk_source
    Seurat::Idents(seu.atac) =  pseudobulk_source
    
    # if (doDimPlots) {
    #   pdfFile = paste0(outputDir, "/dimPlot_combined_", pseudobulk_source, ".pdf")
    #   .doDimPlots(seu.s, groupBy = pseudobulk_source, pdfFile = pdfFile)
    # }
    
    file_metadata = paste0(outputDir, "/metadata_", pseudobulk_source, ".tsv.gz")
    .writeMetadata(seu.rna, file_metadata)
    
    
    # RNA #
    futile.logger::flog.info(paste0(" Aggregate and prepare RNA counts for each cluster"))
    rna.pseudobulk.clean = .aggregateCounts(seu.rna, assayName_RNA, groupBy = pseudobulk_source, ID_column = "gene", sumCounts = sumCounts)
    
    # Merge with the actual features and their correct mappings.
    rna.pseudobulk.clean2 = .addEnsemblIDs(rna.pseudobulk.clean, mapping = features)
    
    futile.logger::flog.info(paste0(" Writing RNA counts to file ", file_in_rna))
    readr::write_tsv(rna.pseudobulk.clean2, file_in_rna)
    
    .printScarcity(rna.pseudobulk.clean2, type = "RNA")
    
    
    # ATAC #
    futile.logger::flog.info(paste0(" Aggregate and prepare ATAC counts for each cluster"))
    atac.pseudobulk.clean = .aggregateCounts(seu.atac, assayName_ATAC, groupBy = pseudobulk_source, sumCounts = sumCounts, ID_column = "peakID")
    
    # Replace the first hyphen with a colon
    atac.pseudobulk.clean$peakID = sub("-", ":", atac.pseudobulk.clean$peakID)
    
    futile.logger::flog.info(paste0(" Writing ATAC counts to file ", file_in_atac))
    readr::write_tsv(atac.pseudobulk.clean, file_in_atac)
    
    .printScarcity(atac.pseudobulk.clean, type = "ATAC")
    
  
  
  
  if (saveSeuratObject) {
    
    file_seurat_rna = paste0(outputDir, "/seuratObject.rna.qs")
    qs::qsave(seu.rna,  file_seurat)
    file_seurat_atac = paste0(outputDir, "/seuratObject.atac.qs")
    qs::qsave(seu.atac,  file_seurat)
    
  } 
  
  GRaNIE:::.printExecutionTime(start, prefix = "") 
  futile.logger::flog.info(paste0("Finished successfully, all output files have been saved in ", outputDir, ". Returning the Seurat object. Happy GRaNIE'ing!"))
  

}

#' @import ggplot2
.doDimPlots <- function(seu.s, groupBy = "seurat_clusters", pdfFile) {
  
  
  if (!is.null(pdfFile)) {
    futile.logger::flog.info(paste0(" Writing DimPlot to file ", pdfFile))
    grDevices::pdf(pdfFile, width = 15, height = 10)
  }
  
  p1 <- Seurat::DimPlot(seu.s, reduction = "umap.rna", group.by = groupBy, label = TRUE, label.size = 2.5, repel = TRUE)  + ggtitle("RNA")
  p2 <- Seurat::DimPlot(seu.s, reduction = "umap.atac", group.by = groupBy, label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
  p3 <- Seurat::DimPlot(seu.s, reduction = "wnn.umap", group.by = groupBy, label = TRUE, label.size = 2.5, repel = TRUE)  + ggtitle("WNN")
  plot(p1 + p2 + p3 & Seurat::NoLegend() & theme(plot.title = element_text(hjust = 0.5)))
  
  if (!is.null(pdfFile)) {
    grDevices::dev.off()
  }
  
  
}

.writeMetadata <- function(seu.s, file_metadata) {
  
  nClusters = levels(seu.s)
  futile.logger::flog.info(paste0(" Number of clusters found: ",  length(nClusters)))
  
  futile.logger::flog.info(paste0(" Cells per cluster: "), table(Seurat::Idents(seu.s)), capture = TRUE)
  
  metadata = tibble::tibble(sampleName = levels(seu.s), 
                            nCells = table(Seurat::Idents(seu.s)) %>% as.vector())
  
  
  futile.logger::flog.info(paste0(" Writing metadata to file ", file_metadata))
  readr::write_tsv(metadata, file_metadata)
}

.aggregateCounts <- function(seu.s, assayName, groupBy = "ident", sumCounts = TRUE, slotName = "counts", ID_column = "peakID") {
  
  checkmate::assertSubset(slotName, c("scale.data", "counts", "data"))
  
  if (sumCounts) {
    
    pseudobulk = Seurat::AggregateExpression(
      seu.s,
      assays = assayName,
      features = NULL,
      return.seurat = FALSE,
      group.by = groupBy,
      add.ident = NULL,
      slot =  slotName, # If slot is set to 'data', this function assumes that the data has been log normalized and therefore feature values are exponentiated prior to aggregating so that sum is done in non-log space. Otherwise, if slot is set to either 'counts' or 'scale.data', no exponentiation is performed prior to aggregating
    ) 
    
  } else {
    
    # seu.s[["SCT"]]@scale.data contains the residuals (normalized values), and is used directly as input to PCA. 
    # The ‘corrected’ UMI counts are stored in seu.s[["SCT"]]@counts
    # We store log-normalized versions of these corrected counts in seu.s[["SCT"]]@data
    
    # assayRNA = "SCT"
    
    # The residuals for this model are normalized values, and can be positive or negative. 
    # Positive residuals for a given gene in a given cell indicate that we observed more UMIs than expected given the gene’s average expression 
    # in the population and cellular sequencing depth, while negative residuals indicate the converse.
    
    # Take the mean of the normalized values per cluster
    pseudobulk = Seurat::AverageExpression(
      seu.s,
      assays = assayName,
      features = NULL,
      return.seurat = FALSE,
      group.by = groupBy,
      add.ident = NULL,
      slot = "counts", # If slot is set to 'data', this function assumes that the data has been log normalized and therefore feature values are exponentiated prior to aggregating so that sum is done in non-log space. Otherwise, if slot is set to either 'counts' or 'scale.data', no exponentiation is performed prior to aggregating
    ) 
    
  }
  
  # atac.before = GetAssayData(object = seu.s, assay = "ATAC", slot = "counts")
  # atac.before2 = GetAssayData(object = seu.s, assay = "ATAC", slot = "data")
  # atac.before3 = GetAssayData(object = seu.s, assay = "ATAC", slot = "scale.data")
  # 
  # # Sanity check
  # for (i in 1:10) {
  #   stopifnot(identical(sum(atac.before[i,]), sum(atac.pseudobulk$ATAC[i,]))) 
  # }
  
  # Prepare data for GRN format
  pseudobulk[[assayName]] %>%
    as.data.frame() %>%
    tibble::rownames_to_column(ID_column) %>%
    tibble::as_tibble() 
  
  
}


.addEnsemblIDs <- function(rna.pseudobulk.clean, mapping) {
  
  futile.logger::flog.info(paste0(" Mapping gene names to Ensembl IDs using the provided features file..."))
  rna.pseudobulk.clean = dplyr::left_join(rna.pseudobulk.clean, mapping, by = "gene")
  
  rna.pseudobulk.clean.filt =  dplyr::filter(rna.pseudobulk.clean, !is.na(.data$ens_id))
  nRows_missing = nrow(rna.pseudobulk.clean) - nrow(rna.pseudobulk.clean.filt)
  genesNA = rna.pseudobulk.clean$gene[which(!rna.pseudobulk.clean$gene %in% rna.pseudobulk.clean.filt$gene)]
  futile.logger::flog.info(paste0(" Deleted the following ", nRows_missing, " rows because no Ensembl ID could be found: ", paste0(genesNA, collapse = ", ")))
  
  # Check ambiguities for gene names, some gene names may be ambiguous
  geneNameFreq = table(rna.pseudobulk.clean.filt$gene)
  genesDelete = names(which(geneNameFreq > 1))
  
  futile.logger::flog.info(paste0(" Deleted the following ", length(genesDelete), " genes because the gene name was not unique and appeared multiple times: ", paste0(genesDelete, collapse = ",")))
  
  rna.pseudobulk.clean = dplyr::filter(rna.pseudobulk.clean.filt, !.data$gene %in% genesDelete) %>%
    dplyr::select("ens_id", tidyselect::everything(), -"view", -"chr", -"start", -"end", -"gene") %>%
    dplyr::rename(ENSEMBL = "ens_id")
  
  rna.pseudobulk.clean
  
}

.printScarcity <- function(pseudobulk, type = "RNA") {
  
  checkmate::assertSubset(type, c("RNA", "ATAC"))
  colname = dplyr::if_else(type == "RNA", "ENSEMBL", "peakID")
  
  pseudobulk.m = pseudobulk %>%
    tibble::column_to_rownames(colname) %>%
    as.matrix()
  
  spasity_fraction = (length(pseudobulk.m) - Matrix::nnzero(pseudobulk.m)) / length(pseudobulk.m)
  
  futile.logger::flog.info(paste0(" Sparsity ", type, ": ",  spasity_fraction))
}

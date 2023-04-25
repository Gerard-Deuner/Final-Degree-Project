############################################
# HELPER FUNCTIONS FOR SINGLE-CELL GRaNIE  #
############################################

# regress out mt.percent and ribo.perncent of SCT
# choose between pearson or spearman correlation

# prepapre the date (in normalisation step I regress out both mitochondrial and ribosomal genes)
prepareSeuratData_GRaNIE <- function(seu.s, outputDir = "pseudobulk", saveSeuratObject = TRUE,
                                     file_RNA_features = "/g/zaugg/carnold/Projects/GRN_pipeline/misc/singleCell/sharedMetadata/features_RNA_hg38.tsv.gz", 
                                     assayName_RNA = "RNA", assayName_ATAC= "ATAC", 
                                     prepareData = TRUE, SCT_nDimensions = 50,  dimensionsToIgnore_LSI_ATAC = 1,
                                     pseudobulk_source = "cluster",
                                     clusteringAlgorithm = 3, clusterResolutions = c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2)), 
                                     doDimPlots = TRUE,
                                     run_MultiK = FALSE
) {
  
  start = Sys.time()
  
  library(tidyverse)
  library(Seurat)
  library(checkmate)
  library(futile.logger)
  
  checkmate::assertClass(seu.s, "Seurat")
  checkmate::assertSubset(c(assayName_RNA, assayName_ATAC), names(seu.s@assays))
  checkmate::assertFlag(prepareData)
  checkmate::assertFileExists(file_RNA_features, access = "r")
  checkmate::assertSubset(pseudobulk_source, c("cluster", colnames(seu.s@meta.data)))
  
  checkmate::assertFlag(saveSeuratObject)
  checkmate::assertFlag(run_MultiK)
  checkmate::assertSubset(clusteringAlgorithm, c(1:4), empty.ok = FALSE)
  checkmate::assertIntegerish(clusteringAlgorithm, len = 1)
  checkmate::assertFlag(doDimPlots)
  checkmate::assertIntegerish(SCT_nDimensions, lower = 5, upper = 500)
  checkmate::assertIntegerish(dimensionsToIgnore_LSI_ATAC, lower = 0, upper = SCT_nDimensions - 1)
  
  if (!dir.exists(outputDir)) {
    dir.create(outputDir)
  }
  
  .startLogger(paste0(outputDir, "/prepareData_Seurat_GRaNIE.R.log") , "INFO",  removeOldLog = FALSE)
  
  if (saveSeuratObject & !is.installed("qs")) {
    stop("Package qs is not installed but needed due to saveSeuratObject = TRUE")
  }
  
  if (saveSeuratObject) {
    library(qs)
  }
  
  
  
  # Currently hard-coded options
  # See SCT transform
  returnOnlyVarGenes = FALSE
  
  # Needed to retrieve ENSEMBL IDs for gene names
  if (!is.null(file_RNA_features)) {
    futile.logger::flog.info(paste0("Reading and checking the RNA features file."))
    
    features = read_tsv(file_RNA_features, col_names = FALSE, col_types = "ccffii")
    colnames(features) = c("ens_id", "gene", "view", "chr", "start", "end")
    
  }
  
  
  
  if (prepareData) {
    
    futile.logger::flog.info(paste0("Preparing RNA and ATAC data (RNA: SCTransform, PCA, UMPA; ATAC: TFID, FindTopFeatures, SVD, UMAP; combined: FindMultiModalNeighbors, UMAP). This may take a while."))
    
    DefaultAssay(seu.s) <- assayName_RNA
    
    dimToUseRNA = seq_len(SCT_nDimensions)
    
    # Please note that this matrix is non-sparse, and can therefore take up a lot of memory if stored for all genes. 
    # To save memory, we store these values only for variable genes, by setting the return.only.var.genes = TRUE by default in the SCTransform() function call.
    seu.s <- SCTransform(seu.s, verbose = FALSE, return.only.var.genes = returnOnlyVarGenes, vars.to.regress = c("percent.mt", "percent.ribo")) %>% 
      RunPCA() %>% 
      RunUMAP(dims = dimToUseRNA, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
    
    
    DefaultAssay(seu.s) <- assayName_ATAC
    seu.s <- RunTFIDF(seu.s)
    seu.s <- FindTopFeatures(seu.s, min.cutoff = 'q0')
    seu.s <- RunSVD(seu.s)
    
    
    dimToUse_ATAC = setdiff(1:SCT_nDimensions, dimensionsToIgnore_LSI_ATAC)
    
    seu.s <- RunUMAP(seu.s, reduction = 'lsi', dims = dimToUse_ATAC, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
    
    # We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering
    seu.s <- FindMultiModalNeighbors(seu.s, reduction.list = list("pca", "lsi"), dims.list = list(dimToUseRNA, dimToUse_ATAC),
                                     weighted.nn.name = "weighted.nn",
                                     knn.graph.name = "wknn",
                                     snn.graph.name = "wsnn")
    
    seu.s <- RunUMAP(seu.s, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
    
    
  } # end if prepareData
  
  # Currently the only opton that works. Count aggregation
  sumCounts  = TRUE
  aggregationType = if_else(sumCounts == TRUE, "", "_mean")
  
  if (pseudobulk_source == "cluster") {
    
    futile.logger::flog.info(paste0("Running clustering for the following resolutions: ", paste0(clusterResolutions, collapse = ", "), ". This may take a while."))
    
    for (clusterResolution in unique(clusterResolutions)) {
      
      file_in_rna  = paste0(outputDir, "/rna.pseudobulkFromClusters_res", clusterResolution, aggregationType, ".tsv.gz")
      file_in_atac = paste0(outputDir, "/atac.pseudobulkFromClusters_res", clusterResolution, aggregationType, ".tsv.gz")
      
      if (file.exists(file_in_rna) & file.exists(file_in_atac)) {
        next
      }
      futile.logger::flog.info(paste0(" Resolution ", clusterResolution))
      
      seu.s <- FindClusters(seu.s, graph.name = "wsnn", 
                            resolution = clusterResolution,
                            algorithm = clusteringAlgorithm, verbose = FALSE)
      
      clusterColName = paste0("wsnn_res.", clusterResolution)
      Idents(seu.s) = seu.s[[clusterColName]] # before: clusterColName
      seu.s@meta.data[[clusterColName]] = paste0("cluster", seu.s@meta.data[[clusterColName]])
      seu.s@meta.data$seurat_clusters = paste0("cluster", seu.s@meta.data$seurat_clusters)
      
      if (doDimPlots) {
        pdfFile = paste0(outputDir, "/dimPlot_combined_clusterResolution", clusterResolution, ".pdf")
        .doDimPlots(seu.s, groupBy = "seurat_clusters", pdfFile = pdfFile)
      }
      
      file_metadata = paste0(outputDir, "/metadata_res", clusterResolution, ".tsv.gz")
      .writeMetadata(seu.s, file_metadata)
      
      
      ### Count aggregation
      futile.logger::flog.info(paste0(" Aggregate and prepare RNA counts for each cluster"))
      rna.pseudobulk.clean = .aggregateCounts(seu.s, assayName_RNA, groupBy = clusterColName, sumCounts, ID_column = "gene") #previous: groupBy = "ident"
      
      # rna.before = GetAssayData(object = seu.s, assay = "RNA", slot = "counts")
      # rna.before2 = GetAssayData(object = seu.s, assay = "RNA", slot = "data")
      # rna.before3 = GetAssayData(object = seu.s, assay = "RNA", slot = "scale.data")
      # 
      # # Sanity check
      # for (i in 1:10) {
      #   stopifnot(identical(sum(rna.before[i,]), sum(rna.pseudobulk$RNA[i,]))) 
      # }
      
      # Merge with the actual features and their correct mappings.
      rna.pseudobulk.clean = .addEnsemblIDs(rna.pseudobulk.clean, mapping = features)
      
      futile.logger::flog.info(paste0(" Writing RNA counts to file ", file_in_rna))
      write_tsv(rna.pseudobulk.clean, file_in_rna)
      
      .printScarcity(rna.pseudobulk.clean, type = "RNA")
      
      ########
      # ATAC #
      ########
      
      futile.logger::flog.info(paste0(" Aggregate and prepare ATAC counts for each cluster"))
      atac.pseudobulk.clean = .aggregateCounts(seu.s, assayName_ATAC, groupBy = clusterColName, sumCounts, ID_column= "peakID") #previous: groupBy = "ident"
      
      
      # Replace the first hyphen with a colon
      atac.pseudobulk.clean$peakID = sub("-", ":", atac.pseudobulk.clean$peakID)
      
      futile.logger::flog.info(paste0(" Writing ATAC counts to file ", file_in_atac))
      write_tsv(atac.pseudobulk.clean, file_in_atac)
      
      .printScarcity(atac.pseudobulk.clean, type = "ATAC")
      
      
    } # end for (clusterResolution in clusterResolutions)
    
    if (run_MultiK) {
      
      futile.logger::flog.info(paste0(" Run MultiK"))
      
      library(MultiK)
      multik <- MultiK(seu.s, reps=10)
      DiagMultiKPlot(multik$k, multik$consensus)
      bestk = 3
      # CLustering
      clusters <- getClusters(seu.s, bestk)
      # Run SigClust at optimal K level:
      pval <- CalcSigClust(seu.s, clusters$clusters)
      # Make diagnostic plots (this includes a dendrogram of cluster centroids with the pairwise SigClust p values mapped on the nodes, 
      # and a heatmap of the pairwise SigClust p values)
      PlotSigClust(seu.s, clusters$clusters, pval)
    }
    
    
    
  }  # end if (pseudobulk_source == "cluster")
  else {
    
    futile.logger::flog.info(paste0("Creating pseudobulk based on the following pre-existing cell identity column: ", pseudobulk_source))
    
    file_in_rna  = paste0(outputDir, "/rna.pseudobulkFromClusters_" , pseudobulk_source, aggregationType, ".tsv.gz")
    file_in_atac = paste0(outputDir, "/atac.pseudobulkFromClusters_", pseudobulk_source, aggregationType, ".tsv.gz")
    
    # CAREFUL WITH THIS!
    if (file.exists(file_in_rna) & file.exists(file_in_atac)) {
      next
    }
    
    # Normalize levels and make the names compatible
    seu.s@meta.data[[pseudobulk_source]] = gsub(" ", "_", seu.s@meta.data[[pseudobulk_source]])
    seu.s@meta.data[[pseudobulk_source]] = gsub("-", "_", seu.s@meta.data[[pseudobulk_source]])
    
    
    Idents(seu.s) = seu.s[[pseudobulk_source]]
    
    if (doDimPlots) {
      pdfFile = paste0(outputDir, "/dimPlot_combined_", pseudobulk_source, ".pdf")
      .doDimPlots(seu.s, groupBy = pseudobulk_source, pdfFile = pdfFile)
    }
    
    file_metadata = paste0(outputDir, "/metadata_", pseudobulk_source, ".tsv.gz")
    .writeMetadata(seu.s, file_metadata)
    
    
    # RNA #
    futile.logger::flog.info(paste0(" Aggregate and prepare RNA counts for each cluster"))
    rna.pseudobulk.clean = .aggregateCounts(seu.s, assayName_RNA, groupBy = pseudobulk_source, ID_column = "gene", sumCounts) #previously: groupBy = "ident"
    
    # Merge with the actual features and their correct mappings.
    rna.pseudobulk.clean2 = .addEnsemblIDs(rna.pseudobulk.clean, mapping = features)
    
    futile.logger::flog.info(paste0(" Writing RNA counts to file ", file_in_rna))
    write_tsv(rna.pseudobulk.clean2, file_in_rna)
    
    .printScarcity(rna.pseudobulk.clean2, type = "RNA")
    
    
    # ATAC #
    futile.logger::flog.info(paste0(" Aggregate and prepare ATAC counts for each cluster"))
    atac.pseudobulk.clean = .aggregateCounts(seu.s, assayName_ATAC, groupBy = pseudobulk_source, sumCounts, ID_column= "peakID") #previously: groupBy = "ident"
    
    # Replace the first hyphen with a colon
    atac.pseudobulk.clean$peakID = sub("-", ":", atac.pseudobulk.clean$peakID)
    
    futile.logger::flog.info(paste0(" Writing ATAC counts to file ", file_in_atac))
    # added by me so rna and atac sample names match 
    #names(atac.pseudobulk.clean)[2:ncol(atac.pseudobulk.clean)] <- levels(seu.s)
    write_tsv(atac.pseudobulk.clean, file_in_atac)
    
    .printScarcity(atac.pseudobulk.clean, type = "ATAC")
    
  }
  
  
  if (saveSeuratObject) {
    
    file_seurat = paste0(outputDir, "/seuratObject.qs")
    qsave(seu.s,  file_seurat)
    
  } 
  
  .printExecutionTime(start, prefix = "") 
  futile.logger::flog.info(paste0("Finished successfully, all output files have been saved in ", outputDir, ". Returning the Seurat object. Happy GRaNIE'ing!"))
  
  seu.s
}

.doDimPlots <- function(seu.s, groupBy = "seurat_clusters", pdfFile) {
  
  
  if (!is.null(pdfFile)) {
    futile.logger::flog.info(paste0(" Writing DimPlot to file ", pdfFile))
    pdf(pdfFile, width = 15, height = 10)
  }
  
  p1 <- DimPlot(seu.s, reduction = "umap.rna", group.by = groupBy, label = TRUE, label.size = 2.5, repel = TRUE)  + ggtitle("RNA")
  p2 <- DimPlot(seu.s, reduction = "umap.atac", group.by = groupBy, label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
  p3 <- DimPlot(seu.s, reduction = "wnn.umap", group.by = groupBy, label = TRUE, label.size = 2.5, repel = TRUE)  + ggtitle("WNN")
  plot(p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5)))
  
  if (!is.null(pdfFile)) {
    dev.off()
  }
  
  
}

.writeMetadata <- function(seu.s, file_metadata) {
  
  nClusters = levels(seu.s)
  futile.logger::flog.info(paste0(" Number of clusters found: ",  length(nClusters)))
  
  futile.logger::flog.info(paste0(" Cells per cluster: "), table(Idents(seu.s)), capture = TRUE)
  
  metadata = tibble(sampleName = levels(seu.s), 
                    nCells = table(Idents(seu.s)) %>% as.vector())
  
  
  futile.logger::flog.info(paste0(" Writing metadata to file ", file_metadata))
  write_tsv(metadata, file_metadata)
}
                                                          
.aggregateCounts <- function (seu.s, assayName, groupBy = "ident", sumCounts, slotName = "counts", ID_column = "peakID") {
  
  checkmate::assertSubset(slotName, c("scale.data", "counts", "data"))
  
  if (sumCounts) {
    
    pseudobulk = AggregateExpression(
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
    pseudobulk = AverageExpression(
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
    rownames_to_column(ID_column) %>%
    as_tibble() 
  
  
}


.addEnsemblIDs <- function(rna.pseudobulk.clean, mapping) {
  
  futile.logger::flog.info(paste0(" Mapping gene names to Ensembl IDs using the provided features file..."))
  rna.pseudobulk.clean = left_join(rna.pseudobulk.clean, mapping, by = "gene")
  
  rna.pseudobulk.clean.filt =  dplyr::filter(rna.pseudobulk.clean, !is.na(ens_id))
  nRows_missing = nrow(rna.pseudobulk.clean) - nrow(rna.pseudobulk.clean.filt)
  genesNA = rna.pseudobulk.clean$gene[which(! rna.pseudobulk.clean$gene %in% rna.pseudobulk.clean.filt$gene)]
  futile.logger::flog.info(paste0(" Deleted the following ", nRows_missing, " rows because no Ensembl ID could be found: ", paste0(genesNA, collapse = ", ")))
  
  # Check ambiguities for gene names, some gene names may be ambiguous
  geneNameFreq = table(rna.pseudobulk.clean.filt$gene)
  genesDelete = names(which(geneNameFreq > 1))
  
  futile.logger::flog.info(paste0(" Deleted the following ", length(genesDelete), " genes because the gene name was not unique and appeared multiple times: ", paste0(genesDelete, collapse= ",")))
  
  rna.pseudobulk.clean = dplyr::filter(rna.pseudobulk.clean.filt, !gene %in% genesDelete) %>%
    dplyr::select(ens_id, everything(), -view, -chr, -start, -end, -gene) %>%
    dplyr::rename(ENSEMBL = ens_id)
  
  rna.pseudobulk.clean
  
}

.printScarcity <- function (pseudobulk, type = "RNA") {
  
  checkmate::assertSubset(type, c("RNA", "ATAC"))
  colname = dplyr::if_else(type == "RNA", "ENSEMBL", "peakID")
  
  pseudobulk.m = pseudobulk %>%
    column_to_rownames(colname) %>%
    as.matrix()
  
  spasity_fraction = (length(pseudobulk.m) - Matrix::nnzero(pseudobulk.m)) / length(pseudobulk.m)
  
  flog.info(paste0(" Sparsity ", type, ": ",  spasity_fraction))
}



.printExecutionTime <- function(startTime, prefix = " ") {
  
  endTime  <-  Sys.time()
  futile.logger::flog.info(paste0(prefix, "Finished successfully. Execution time: ", round(endTime - startTime, 1), " ", units(endTime - startTime)))
}


#' @import utils
is.installed <- function(mypkg){
  is.element(mypkg, installed.packages()[,1])
} 


.startLogger <- function(logfile, level = "INFO", removeOldLog = TRUE, appenderName = "consoleAndFile", verbose = TRUE) {
  
  checkmate::assertChoice(level, c("TRACE", "DEBUG", "INFO", "WARN", "ERROR", "FATAL"))
  checkmate::assertFlag(removeOldLog)
  checkmate::assertChoice(appenderName, c("console", "file", "consoleAndFile"))
  checkmate::assertFlag(verbose)
  
  if (appenderName != "console") {
    checkmate::assertDirectory(dirname(logfile), access = "w")
    if (file.exists(logfile)) {
      file.remove(logfile)
    }
  }
  
  # LEVELS: TRACE, DEBUG, INFO, WARN, ERROR, FATAL
  invisible(futile.logger::flog.threshold(level))
  
  
  if (appenderName == "console") {
    invisible(futile.logger::flog.appender(futile.logger::appender.console()))
  } else if (appenderName == "file")  {
    invisible(futile.logger::flog.appender(futile.logger::appender.file(file = logfile)))
  } else {
    invisible(futile.logger::flog.appender(futile.logger::appender.tee(file = logfile)))
  }
  
  
}

runGRaNIE_batchMode <- function (datasetName, 
                                 inputDir, outputDir,
                                 clusterResolutions = c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2)),
                                 idColumn_peaks = "peakID",  idColumn_RNA = "ENSEMBL",
                                 genomeAssembly = "hg38",
                                 TFBS_folder = NULL,
                                 nCores = 8, 
                                 normRNA_all =  c("limma_quantile"),   normATAC_all =  c("DESeq2_sizeFactors"),
                                 includeSexChr = FALSE,
                                 minCV = 0,
                                 minNormalizedMean_peaks = 5,
                                 minNormalizedMean_RNA = 1,
                                 minSizePeaks = 5,
                                 promoterRange = 250000, useGCCorrection = FALSE,
                                 TF_peak.fdr.threshold = 0.2,
                                 peak_gene.fdr.threshold = 0.1,
                                 runNetworkAnalyses = FALSE, 
                                 forceRerun = TRUE,
                                 correlation.method) {
  
  
  library(tidyverse)
  library(qs)
  library(GRaNIE)
  library(futile.logger)
  
  checkmate::assertFlag(includeSexChr)
  checkmate::assertFlag(useGCCorrection)
  checkmate::assertFlag(runNetworkAnalyses)
  checkmate::assertFlag(forceRerun)
  
  
  # normRNA_all =  c("limma_quantile", "DESeq2_sizeFactors")
  # normATAC_all =  c("limma_cyclicloess", "EDASeq_GC_peaks",  "DESeq2_sizeFactors")
  
  
  for (subsetName in c("pseudobulk")) {
    
    for (clusterResolution in clusterResolutions) { 
      
      
      for (normRNA in normRNA_all) {
        
        for (normATAC in normATAC_all ) {
          
          cat("Running GRaNIE for cluster resolution ", clusterResolution, " (normalization RNA: ", normRNA, ", normalization ATAC: ", normATAC)
          
          dir_output = paste0(outputDir, "/output_", subsetName, "_clusterRes", clusterResolution, "_RNA_", normRNA, "_ATAC_", normATAC)
          
          file_in_rna  = paste0(inputDir, "/rna.pseudobulkFromClusters_res", clusterResolution, ".tsv.gz")
          file_in_atac = paste0(inputDir, "/atac.pseudobulkFromClusters_res", clusterResolution, ".tsv.gz")
          file_metadata = paste0(inputDir, "/metadata_res", clusterResolution, ".tsv.gz")
          
          GRN = runGRaNIE(  dir_output, 
                            datasetName = paste0(datasetName, "_pseudobulk_res", clusterResolution),
                            file_peaks = file_in_atac, 
                            file_rna = file_in_rna, 
                            file_metadata = file_metadata,
                            genomeAssembly = genomeAssembly,
                            normalization_peaks = normATAC, 
                            idColumn_peaks = idColumn_peaks,
                            normalization_rna = normRNA, 
                            idColumn_RNA = idColumn_RNA,
                            includeSexChr = includeSexChr,
                            TFBS_folder = TFBS_folder,
                            nCores = nCores, 
                            minCV = minCV,
                            minNormalizedMean_peaks = minNormalizedMean_peaks,
                            minNormalizedMean_RNA = minNormalizedMean_RNA,
                            minSizePeaks = minSizePeaks,
                            promoterRange = promoterRange, 
                            useGCCorrection = useGCCorrection,
                            TF_peak.fdr.threshold = TF_peak.fdr.threshold,
                            peak_gene.fdr.threshold = peak_gene.fdr.threshold,
                            runNetworkAnalyses = runNetworkAnalyses, 
                            forceRerun = forceRerun,
                            correlation.method = corr.method
          )
          
        } # end for ATAC normalization
      } # end for RNA normalization
      
      
    } # end for cluster resolution
    
    
  }
  
}



runGRaNIE <- function(dir_output = "output_GRaNIE", 
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
                      forceRerun = TRUE,
                      correlation.method
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
  
  # added by me so rna and atac sample names match 
  # names(countsATAC)[2:ncol(countsATAC)] <- names(countsRNA)[2:ncol(countsRNA)]
  
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
                               corMethod = correlation.method, maxFDRToStore = 0.3,
                               useGCCorrection = useGCCorrection, percBackground_size = 75, percBackground_resample = TRUE,
                               forceRerun = forceRerun)
  
  
  file_input_TADs = ""
  overlapTypeGene = "TSS"
  
  GRN = addConnections_peak_gene(GRN,
                                 overlapTypeGene = overlapTypeGene,
                                 corMethod = correlation.method, shuffleRNACounts = TRUE,
                                 promoterRange = promoterRange, TADs = NULL,
                                 nCores = nCores, plotDiagnosticPlots = TRUE,
                                 forceRerun = forceRerun)
  
  
  GRN = filterGRNAndConnectGenes(GRN, TF_peak.fdr.threshold = TF_peak.fdr.threshold, 
                                 TF_peak.connectionTypes = "expression" ,
                                 peak_gene.fdr.threshold = peak_gene.fdr.threshold,
                                 gene.types = c("protein_coding"),
                                 allowMissingTFs = FALSE, allowMissingGenes = FALSE,
                                 peak_gene.r_range = c(0,1))
  # run analysis on FDRs
  GRN = generateStatsSummary(GRN, TF_peak.fdr = c(0.05, 0.1, 0.2, 0.25, 0.3), TF_peak.connectionTypes = "all",
                             peak_gene.fdr = c(0.1, 0.2, 0.3), peak_gene.r_range = c(0, 1), allowMissingGenes = c(FALSE,TRUE), 
                             allowMissingTFs = c(FALSE), gene.types = c("protein_coding", "lincRNA"),
                             forceRerun = TRUE)
  
  GRN = plot_stats_connectionSummary(GRN, type = "heatmap", plotAsPDF = TRUE, pages = 3)
  GRN = plot_stats_connectionSummary(GRN, type = "boxplot", plotAsPDF = TRUE, pages = 1)
  
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

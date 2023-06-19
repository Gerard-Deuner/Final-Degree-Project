##############################
# STEP By STEP eGRN ANALYSIS #
##############################

# load libraries
library(Seurat)
library(qs)
library(dplyr)
library(GRaNIE)

# load eGRN
GRN <- qread("/g/scb/zaugg/deuner/GRaNIE/outputdata/combined_pearson_Res_16/output_pseudobulk_wsnn_res.16_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/GRN.qs")

# build eGRN graph
GRN = build_eGRN_graph(GRN, forceRerun = TRUE)

# visualize the eGRN
GRN = visualizeGRN(GRN, plotAsPDF = FALSE, maxEdgesToPlot = 20000)

# general network statistics
GRN = plotGeneralGraphStats(GRN, plotAsPDF = FALSE)

# calculate general enrichment
GRN = calculateGeneralEnrichment(GRN)

# general network enrichment
GRN = plotGeneralEnrichment(GRN, plotAsPDF = FALSE, pages = 1)

# generate graph comunities and their summarizing statistics
GRN = calculateCommunitiesStats(GRN, clustering = "louvain", forceRerun = FALSE)

# plot general structure & connectivity statistics for each community in a filtered GRN
GRN = plotCommunitiesStats(GRN, plotAsPDF = FALSE)

# run an enrichment analysis for the genes in each community in the filtered GRN object
GRN = calculateCommunitiesEnrichment(GRN)

# plot community-based enrichment results for a filtered GRN object
GRN = plotCommunitiesEnrichment(GRN, plotAsPDF = FALSE)

# run an enrichment analysis for the set of genes connected to a particular TF or sets of TFs in the filtered GRN object
GRN = calculateTFEnrichment(GRN)

# plot TF-based GO enrichment results
GRN = plotTFEnrichment(GRN, plotAsPDF = FALSE)


###########################
# AUTOMATIC eGRN ANALYSIS #
###########################

file <- "combined_batch_mode_spearman_nomicro"
for (f in list.files(paste0("/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/", file , "/Batch_Mode_Outputs/"))){
   GRN = qread(paste0("/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/", file , "/Batch_Mode_Outputs/", f, "/GRN.qs"))
   GRN@connections$all.filtered$`0` %>% nrow()
   GRN@connections$all.filtered$`0` <- GRN@connections$all.filtered$`0` %>% filter(TF_peak.fdr < 0.2, peak_gene.p_adj < 0.1) 
   GRN@connections$all.filtered$`0` %>% nrow()
   GRN@connections$all.filtered$`0` %>% View
   GRN = performAllNetworkAnalyses(GRN, ontology = c("GO_BP"), forceRerun = TRUE)
}
  

GRN = qread("/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/combined_batch_mode_spearman_nomicro/Batch_Mode_Outputs/output_pseudobulk_clusterRes8_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/GRN.qs")
GRN = performAllNetworkAnalyses(GRN, ontology = c("GO_BP"), forceRerun = TRUE)
qsave(GRN, file = "/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/combined_batch_mode_spearman_nomicro/Batch_Mode_Outputs/output_pseudobulk_clusterRes8_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/GRN2.qs")

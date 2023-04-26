##############################
# STEP By STEP eGRN ANALYSIS #
##############################

# load libraries
library(Seurat)
library(qs)
library(dplyr)
library(GRaNIE)

# # load eGRN
# GRN <- qread("/g/scb/zaugg/deuner/GRaNIE/outputdata/combined_pearson_Res_16/output_pseudobulk_wsnn_res.16_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/GRN.qs")
# 
# # build eGRN graph
# GRN = build_eGRN_graph(GRN, forceRerun = TRUE)
# 
# # visualize the eGRN
# GRN = visualizeGRN(GRN, plotAsPDF = FALSE, maxEdgesToPlot = 20000)
# 
# # general network statistics
# GRN = plotGeneralGraphStats(GRN, plotAsPDF = FALSE)
# 
# # calculate general enrichment
# GRN = calculateGeneralEnrichment(GRN)
# 
# # general network enrichment
# GRN = plotGeneralEnrichment(GRN, plotAsPDF = FALSE, pages = 1)
# 
# # generate graph comunities and their summarizing statistics
# GRN = calculateCommunitiesStats(GRN, clustering = "louvain", forceRerun = FALSE)
# 
# # plot general structure & connectivity statistics for each community in a filtered GRN
# GRN = plotCommunitiesStats(GRN, plotAsPDF = FALSE)
# 
# # run an enrichment analysis for the genes in each community in the filtered GRN object
# GRN = calculateCommunitiesEnrichment(GRN)
# 
# # plot community-based enrichment results for a filtered GRN object
# GRN = plotCommunitiesEnrichment(GRN, plotAsPDF = FALSE)
# 
# # run an enrichment analysis for the set of genes connected to a particular TF or sets of TFs in the filtered GRN object
# GRN = calculateTFEnrichment(GRN)
# 
# # plot TF-based GO enrichment results
# GRN = plotTFEnrichment(GRN, plotAsPDF = FALSE)

###########################
# AUTOMATIC eGRN ANALYSIS #
###########################

for (file1 in list.files("/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/")) {
  for (file2 in list.files(paste0("/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/", file1 , "/Batch_Mode_Outputs/"))){
    GRN = qread(paste0("/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/", file1 , "/Batch_Mode_Outputs/", file2, "/GRN.qs"))
    GRN = performAllNetworkAnalyses(GRN, ontology = c("GO_BP"), forceRerun = TRUE)
  }
}
  

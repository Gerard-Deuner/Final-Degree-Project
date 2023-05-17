##############################
# CHOSEN RESOLUTION ANALYSIS #
##############################

# load libraries
library(qs)
library(GRaNIE)
library(dplyr)

# set up GRN directory
dir <- "/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/combined_batch_mode_spearman_nomicro/Batch_Mode_Outputs/output_pseudobulk_clusterRes8_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/"

# load GRN object
GRN = qread(paste0(dir, "GRN.qs"))

# construct the eGRN graph
GRN = build_eGRN_graph(GRN, forceRerun = TRUE)

# visualize the filtered eGRN
GRN = visualizeGRN(GRN, plotAsPDF = FALSE, maxEdgesToPlot = 1200)
GRN = visualizeGRN(GRN, plotAsPDF = FALSE, maxEdgesToPlot = 1200, layout = "kk", colorby = "community", nCommunitiesMax = 15)
GRN = visualizeGRN(GRN, plotAsPDF = FALSE, maxEdgesToPlot = 1200, layout = "kk", colorby = "community", nCommunitiesMax = 60)
GRN = visualizeGRN(GRN, plotAsPDF = FALSE, maxEdgesToPlot = 1200, layout = "kk")
GRN = visualizeGRN(GRN, plotAsPDF = FALSE, maxEdgesToPlot = 1600, graph = "TF-peak-gene")
GRN = visualizeGRN(GRN, plotAsPDF = FALSE, maxEdgesToPlot = 1600, graph = "TF-peak-gene", layout = "kk")

# TF enrichment
GRN = plotTFEnrichment(GRN, plotAsPDF = FALSE, n = 3, pages = c(5), maxWidth_nchar_plot = 150)



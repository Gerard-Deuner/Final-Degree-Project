########################
# GRaNIE GRNs ANALYSIS #
########################

library(qs)
library(GRaNIE)
library(dplyr)

#############################################################################################################################################
# Load GRN object from GRaNIE ran on timecourse dataset using the redefined clusters as pseudobulk source and standard GRaNIE normalization #
#############################################################################################################################################

# load eGRN
GRN1 <- qread("/g/scb/zaugg/deuner/GRaNIE/outputdata/timecourse_celltype/output_pseudobulk_celltype_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/GRN.qs")
GRN1

# build eGRN
GRN1 = build_eGRN_graph(GRN1, forceRerun = TRUE)
GRN1 = visualizeGRN(GRN1, plotAsPDF = FALSE)

# get nr. of unique TFs, nr. of unique peaks, nr. of unique genes, nr of TF-peak links, nr. of peak-gene links, nr. of TF-peak-gene links

GRN1.uniqueTFs <- GRN1@connections$all.filtered$`0` %>% pull("TF.name") %>% unique() %>% as.data.frame()
GRN1.uniquePEAKs <- GRN1@connections$all.filtered$`0` %>% pull("peak.ID") %>% unique() %>% as.data.frame()
GRN1.uniqueGENEs <- GRN1@connections$all.filtered$`0` %>% pull("gene.name") %>% unique() %>% as.data.frame()
GRN1.TF_peak.links <- GRN1@connections$all.filtered$`0` %>% as.data.frame() %>% dplyr::select(c("TF.name", "peak.ID")) %>% distinct()
GRN1.peak_gene.links <- GRN1@connections$all.filtered$`0` %>% as.data.frame() %>% dplyr::select(c("gene.name", "peak.ID")) %>% distinct()
GRN1.TF_peak_gene.links <- GRN1@connections$all.filtered$`0` %>% as.data.frame() %>% dplyr::select(c("TF.name", "peak.ID", "gene.name")) %>% distinct()

#rna.data1 <- read_tsv("/g/scb/zaugg/deuner/GRaNIE/outputdata/Timecourse_redefined_celltype_diff_res.0.5/rna.pseudobulkFromClusters_celltype.tsv.gz")


#################################################################################################################################
# Load GRN object from GRaNIE ran on timecourse dataset using the redefined clusters as pseudobulk source and preprocessed data #
#################################################################################################################################

# load eGRN
GRN2 <- qread("/g/scb/zaugg/deuner/GRaNIE/outputdata/Timecourse_redefined_celltype_preprocessedData/output_pseudobulk_redefined_celltype_pp/GRN.qs")
GRN2

# build eGRN
GRN2 = build_eGRN_graph(GRN2, forceRerun = TRUE)
GRN2 = visualizeGRN(GRN2, plotAsPDF = FALSE)

GRN2@data$peaks$counts

# get nr. of unique TFs, nr. of unique peaks, nr. of unique genes, nr of TF-peak links, nr. of peak-gene links, nr. of TF-peak-gene links

GRN2.uniqueTFs <- GRN2@connections$all.filtered$`0` %>% pull("TF.name") %>% unique() %>% as.data.frame()
GRN2.uniquePEAKs <- GRN2@connections$all.filtered$`0` %>% pull("peak.ID") %>% unique() %>% as.data.frame()
GRN2.uniqueGENEs <- GRN2@connections$all.filtered$`0` %>% pull("gene.name") %>% unique() %>% as.data.frame()
GRN2.TF_peak.links <- GRN2@connections$all.filtered$`0` %>% as.data.frame() %>% dplyr::select(c("TF.name", "peak.ID")) %>% distinct()
GRN2.peak_gene.links <- GRN2@connections$all.filtered$`0` %>% as.data.frame() %>% dplyr::select(c("gene.name", "peak.ID")) %>% distinct()
GRN2.TF_peak_gene.links <- GRN2@connections$all.filtered$`0` %>% as.data.frame() %>% dplyr::select(c("TF.name", "peak.ID", "gene.name")) %>% distinct()

# rna.data <- read_tsv("/g/scb/zaugg/deuner/GRaNIE/outputdata/Timecourse_redefined_celltype_preprocessedData/rna.pseudobulkFromClusters_celltype.tsv.gz")


###################
# VISUALIZE STATS #
###################

# dataframe containing eGRNs stats
eGRN.stats <- data.frame("eGRN" = c("GRN1", "GRN2"),
                         "unique TFs" = c(nrow(GRN1.uniqueTFs), nrow(GRN2.uniqueTFs)),
                         "unique PEAKS" = c(nrow(GRN1.uniquePEAKs), nrow(GRN2.uniquePEAKs)),
                         "unique GENEs" = c(nrow(GRN1.uniqueGENEs), nrow(GRN2.uniqueGENEs)),
                         "TF-PEAK links" = c(nrow(GRN1.TF_peak.links), nrow(GRN2.TF_peak.links)),
                         "PEAK-GENE links" = c(nrow(GRN1.peak_gene.links), nrow(GRN2.peak_gene.links)),
                         "TF-PEAK-GENE links" = c(nrow(GRN1.TF_peak_gene.links), nrow(GRN2.TF_peak_gene.links))
                         )
eGRN.stats

# convert it into long format
eGRN.stats.long <- tidyr::gather(
                      eGRN.stats,
                      - eGRN,
                      key =
                        "Stat",
                      value =
                        "Value"
                      )
eGRN.stats.long

# barplot of the metrics
ggplot(eGRN.stats.long, aes(Stat, Value, fill = eGRN)) +
  geom_bar(stat="identity", position = "dodge", color = "black") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10)) + 
  labs(x = "eGRNs Metrics", title = "eGRN Metrics comparison between different normalisations", caption = "GRN1 = standard preprocessing, GRN2 = accurate preprocessing") + 
  scale_fill_manual(values = c("#98B4D4", "#FF847C"))
  

# See which Genes, TFs and Peaks overlap between both eGRNs
intersect(GRN1.uniqueTFs, GRN2.uniqueTFs)

GRN2@data$RNA
GRN2@data$peaks


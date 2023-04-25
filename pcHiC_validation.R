#############################
# pcHi-C Validation script #
#############################

## Description: This script performs promoter capture Hi-C validation on the resulting peak to gene links to any gene regulatory network.
## Input: path of the csv file (or similar) containing peak-gene links information / path of the GRNs to analyse 
## Output: plots showing the analyses

## Data cell types: Neuron, NPC and iPSC

# load libraries
library(dplyr)
library(qs)
library(Seurat)
library(ggplot2)
library(pROC)

# read arguments from command line
args <- commandArgs(trailingOnly = TRUE)

# set up dataset
dataset <- args[1] # timecourse | combined

# set up correlation method
corr.method <- args[2] # pearson | spearman

# set up which cell type links are used for validation
cell.type <- args[3] # neuron | NPC | iPSC | all

# set output dir
out.dir <- paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/results/", cell.type, "_", dataset, "_", corr.method, "/")
dir.create(out.dir) # create dir if it doesn't exist already

# path of GRNs
GRNs.dir <- paste0("/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/", dataset, "_batch_mode_", corr.method, "/Batch_Mode_Outputs/")

# PCHiC links path
links.dir <- paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/links/", dataset, "_", cell.type, "_pchic_links.tsv")
links <- read.csv(links.dir) 

# subset just useful columns 
links <- links %>%
  dplyr::select(gene, peak)
head(links)

# convert gene symbols to ENSEMBL IDs
require('biomaRt')
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'hgnc_symbol',
    'ensembl_gene_id',
    'gene_biotype'),
  uniqueRows = TRUE)
colnames(annotLookup)[1] <- "gene"
links <- links %>% inner_join(annotLookup, by = "gene", multiple = "all") 
colnames(links)[3] <- "gene.ENSEMBL"
head(annotLookup)
head(links)

# set vector with the resolutions
resolutions <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))

# set colors for ROC curves
colors <- c("#FFC107", "#2196F3", "#4CAF50", "#FF5722", "#9C27B0", "#3F51B5", "#FFEB3B", "#8BC34A", "#E91E63", 
            "#673AB7", "#00BCD4", "#FF9800", "#CDDC39", "#795548", "#607D8B", "#9E9E9E", "#F44336", "#00ACC1")
                     
# set vecor of auc positioning in y axis
aucy <- seq(10, 110, by = 5)

# iterate over GRNs
j = 0

# boolean for first plot
first = TRUE

# vectors to store TPs and FPs
TP.vec <- c()
FP.vec <- c()

# dataframe where all peak-gene links will be stored (independent of the resolution)
TP.all <- data.frame()

# list where all the GRN links dfs will be stored no matter if they are validated or not
GRN.links.all.list <- list()

for (res in resolutions){
  # set index
  j <- j + 1
  
  # set up GRN directory
  GRN.dir <- paste0(GRNs.dir, "output_pseudobulk_clusterRes", res, "_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/")
  # read GRN object
  GRN <- qread(paste0(GRN.dir, "GRN.qs"))
  
  # get gene-peak connections from GRN
  #GRN_links <- GRN@connections$all.filtered$`0` %>% as.data.frame() %>% dplyr::select(peak.ID, gene.name, peak_gene.r, peak_gene.p_raw)
  GRN.links <- GRN@connections$peak_genes$`0` %>% dplyr::filter(peak_gene.p_raw < 0.05) %>% as.data.frame() %>% dplyr::select(peak.ID, gene.ENSEMBL, peak_gene.r, peak_gene.p_raw)
  
  # adapt peaks format to links format
  form.peak <- rep("", nrow(GRN.links))
  for (i in 1:nrow(GRN.links)){
    peak <- as.character(GRN.links$peak.ID[i])
    split.peak <- strsplit(peak, split = ":")
    new.peak <- paste(split.peak[[1]][1], split.peak[[1]][2], sep = "-")
    form.peak[i] <- new.peak
  }
  GRN.links$id <- seq(1:nrow(GRN.links))
  GRN.links$peak.ID <- form.peak 
  names(GRN.links)[1:2] <- c("peak", "gene.ENSEMBL")
  
  # Binary classification
  # get rows that match in GRN links and HiC links by gene name and peak id
  TP <- inner_join(GRN.links, links, by = c("peak", "gene.ENSEMBL"), multiple = "all")
  TP$res <- rep(resolutions[j], nrow(TP))
  TP.all <- rbind(TP.all, TP)
  # set 1 if they are in HiC data and 0 if not
  bc <- ifelse(GRN.links$id %in% TP$id, yes = 1, no = 0)
  # add column to GRN_links with that information
  GRN.links$bc <- bc
  # create an confusion matrix 
  conf.m <- table(bc)
  
  # store GRN links data in GRNs list
  GRN.links.all.list <- append(GRN.links.all.list, assign(paste("GRN_links", res, sep = "."), GRN.links))
  
  # Compute ROC curve
  par(pty = "s")
  
  # print metrics
  print("")
  print(paste0("Resolution: ", res))
  print(paste0("GRN links: ", nrow(GRN.links)))
  print(paste0("TP: ", nrow(TP)))
  FP <- nrow(GRN.links) - nrow(TP)
  print(paste0("FP: ", FP))
  print("")
  
  TP.vec <- c(TP.vec, nrow(TP))
  FP.vec <- c(FP.vec, FP)
  
  if ((first == TRUE) & (sum(bc) > 0)){
    ROC <- roc(GRN.links$bc, GRN.links$peak_gene.r, plot = TRUE, legacy.axes = TRUE, percent = TRUE,
               xlab="False Positive Percentage | 1 - Specificity", ylab="True Positive Percentage | Sensitivity",
               col = colors[j], lwd=2, print.auc = TRUE, main = "ROC", add = FALSE, print.auc.x = 115, print.auc.y = aucy[j])
    legend("bottomright",  legend = resolutions, col = colors, ncol = 2)
    first = FALSE
    print(ci.auc(ROC))
  }
  if ((first == FALSE) & (sum(bc) > 0)){ 
    ROC <- roc(GRN.links$bc, GRN.links$peak_gene.r, plot = TRUE, legacy.axes = TRUE, percent = TRUE,
               xlab="False Positive Percentage | 1 - Specificity", ylab="True Positive Percentage | Sensitivity",
               col = colors[j], lwd=2, print.auc = TRUE, main = "ROC", add = TRUE, print.auc.x = 115, print.auc.y = aucy[j])
    legend("bottomright", legend = resolutions, col = colors, lwd = 4, ncol = 2)
    print(ci.auc(ROC))
  }
  
}

# Save all the GRNs links found in the pcHiC data
dir.create(paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/validated_links/"))
GRN.links.all <- do.call("rbind", GRN.links.all.list)
head(GRN.links.all)
write.csv(GRN.links.all, paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/validated_links/", dataset, "_", corr.method, ".csv"))


# barplots of TP vs FP
TPFP.df <- data.frame(TPFP = c(TP.vec, FP.vec), cat = c(rep("TP", 19), rep("FP", 19)), res = c(resolutions <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2)), resolutions <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))))
p1 <- ggplot(TPFP.df, aes(x = as.factor(res), y = TPFP, fill = cat)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  theme_classic() + 
  scale_fill_manual(values = c("#E3B448", "#3A6B35"))
ggsave(paste0(out.dir, "TPvsFP.pdf"), p1, device = "pdf")

# compute ratio of TPs vs FPs
TPFP.vec <- TP.vec/FP.vec
ratio.df <- data.frame(res = resolutions, ratio = TPFP.vec)
p2 <- ggplot(ratio.df, aes(x = as.factor(res), y = ratio)) + 
  geom_bar(stat = "identity", fill = "#993300", color = "black") + 
  theme_classic() + 
  labs(x = "Resolutions", y = "TP / FP", title = "TP and FP ratio per resolution") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
ggsave(paste0(out.dir, "TPvsFP_ratio.pdf"), p2, device = "pdf")

# Check if the TPs are shared 
TP.peak.gene <- TP.all %>% dplyr::select(c("peak", "gene.ENSEMBL", "res"))

# count occurrence of each link
TP.peak.gene.counts <- TP.peak.gene %>% dplyr::select(c("peak", "gene.ENSEMBL", "res")) %>% group_by(peak, gene.ENSEMBL) %>% count 

# create a column for the pasted peak-gene info
TP.peak.gene$link <- paste(TP.peak.gene$peak, TP.peak.gene$gene, sep = "_")

# join counts info to the dataframe
TP.peak.gene <- inner_join(TP.peak.gene, TP.peak.gene.counts, by = c("peak", "gene.ENSEMBL"))

# check whether the link is shared or not
TP.peak.gene$shared <- ifelse(TP.peak.gene$n > 1, TRUE, FALSE)

# plot count of shared links vs non-shared links
p3 <- ggplot(TP.peak.gene %>% distinct(peak, gene.ENSEMBL, .keep_all = TRUE), aes(shared, fill = shared)) + 
  geom_bar(stat = "count") + 
  theme_classic() + 
  scale_fill_manual(values = c("#E3B448", "#3A6B35")) + 
  labs(title = "Are the links shared? (unique links)") + 
  geom_text(stat='count', aes(label=..count..), vjust=-1) + 
  theme(legend.position = "none")

p4 <- ggplot(TP.peak.gene, aes(shared, fill = shared)) + 
  geom_bar(stat = "count") + 
  theme_classic() + 
  scale_fill_manual(values = c("#E3B448", "#3A6B35")) + 
  labs(title = "Are the links shared? (all links)") + 
  geom_text(stat='count', aes(label=..count..), vjust=-1)

ggsave(paste0(out.dir, "shared_links_unique.pdf"), p3, device = "pdf")
ggsave(paste0(out.dir, "shared_links_all.pdf"), p4, device = "pdf")

# plot duplicated vs unique links per cluster resolution
p5 <- ggplot(TP.peak.gene, aes(as.factor(res), fill = shared)) + 
  geom_bar(stat = "count", position = "dodge") + 
  theme_classic() +
  scale_fill_manual(values = c("#E3B448", "#3A6B35")) + 
  labs(title = "Are the links shared?", x = "Cluster Resolutions")
ggsave(paste0(out.dir, "duplicate_vs_unique_links.pdf"), p5, device = "pdf")


# plot most frequent links across resolutions
TP.peak.gene.counts$link <- paste(TP.peak.gene.counts$peak, TP.peak.gene.counts$gene.ENSEMBL, sep = "_")
TP.peak.gene.counts$shared <- ifelse(TP.peak.gene.counts$n > 1, TRUE, FALSE)

p6 <- ggplot(TP.peak.gene.counts %>% arrange(desc(n)) %>% head(10), aes(x = link, y = n, fill = link)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") + 
  scale_fill_viridis_d() + 
  labs(y = "count", title = "Most frequent links across resolutions", caption = "Top 10 Most Frequent Links") 
ggsave(paste0(out.dir, "most_freq_links.pdf"), p6, device = "pdf")

# plot only duplicated links and their corresponding cluster resolutions
library(RColorBrewer)
palette3_info <- brewer.pal.info[brewer.pal.info$category == "qual", ]  # Extract color info
palette3_all <- unlist(mapply(brewer.pal,                     # Create vector with all colors
                              palette3_info$maxcolors,
                              rownames(palette3_info)))
palette3_all 
palette3 <- sample(palette3_all, 19)                    # Sample colors
palette3

p7 <- ggplot(TP.peak.gene %>% dplyr::filter(shared == TRUE) , aes(x = link, fill = as.factor(res))) + 
  geom_bar(stat = "count") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y = "count", title = "To which eGRNs do the shared links belong?", fill = "Resolutions", caption = "35 Different Links") + 
  scale_fill_manual(values = palette3)
ggsave(paste0(out.dir, "links_GRN_mapping.pdf"), p7, device = "pdf")

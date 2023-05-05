#####################################################
# OVERLAP ANALYSIS ON DIFFERENT CLUSTER RESOLUTIONS #
#####################################################

# load libraries
library(GRaNIE)
library(dplyr)
library(qs)
library(ggplot2)
library(ggpubr)
library(pROC)
library(EnsDb.Hsapiens.v79)
library(UpSetR)
library(data.table)
library(ggrepel)
library(stringr)

# Set up path of the cluster resolutions eGRNs
path <- "/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/combined_batch_mode_spearman_nomicro/Batch_Mode_Outputs/"

# set vector with the resolutions
resolutions <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))

# create data.frame to store each GRN metrics
GRN.metrics <- data.frame(matrix(ncol = 7, nrow = 0))
metrics.names <- c("res", "uniqueTFs", "uniquePEAKs", "uniqueGENEs", "TF.PEAK.links", "PEAK.GENE.links", "TF.PEAK.GENE.links")
colnames(GRN.metrics) <- metrics.names

# iterate over GRNs
for (res in resolutions){
  
  # set up GRN directory
  GRN_dir <- paste0(path, "output_pseudobulk_clusterRes", res, "_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/")
  # read GRN object
  GRN <- qread(paste0(GRN_dir, "GRN.qs"))
  
  # retrieve different links from the GRN
  # TF-Peak
  GRN.TF_peak.links <- GRN@connections$all.filtered$`0` %>% as.data.frame() %>% dplyr::select(c("TF.name", "peak.ID")) %>% distinct()
  # Peak-gene
  GRN.peak_gene.links <- GRN@connections$all.filtered$`0` %>% as.data.frame() %>% dplyr::select(c("gene.name", "peak.ID")) %>% distinct()
  # TF-peak-gene 
  GRN.TF_peak_gene.links <- GRN@connections$all.filtered$`0` %>% as.data.frame() %>% dplyr::select(c("TF.name", "peak.ID", "gene.name")) %>% distinct()

  # store the links in independent dataframes (probably there's a better way to do this without using so much memory)
  assign(paste("TF.PEAK.GENE", res, sep = "."), GRN.TF_peak_gene.links)
  assign(paste("PEAK.GENE", res, sep = "."), GRN.peak_gene.links)
  assign(paste("TF.PEAK", res, sep = "."), GRN.TF_peak.links)
  
}


###################
# Peak-Gene links #
###################

# create binary df with all peak-genes links and resolutions as columns
df_list <- list(PEAK.GENE.0.1, PEAK.GENE.0.25, PEAK.GENE.0.5, PEAK.GENE.0.75, PEAK.GENE.1, PEAK.GENE.2, PEAK.GENE.3, PEAK.GENE.4,
                 PEAK.GENE.5, PEAK.GENE.6, PEAK.GENE.7, PEAK.GENE.8, PEAK.GENE.9, PEAK.GENE.10, PEAK.GENE.12, PEAK.GENE.14,
                 PEAK.GENE.16, PEAK.GENE.18, PEAK.GENE.20)

# binary df indicating if link (row) is present in GRN (col)
# a bunary df is required to do an upset plot
peak.gene.links <- data.frame(matrix(ncol = 19, nrow = 0))
colnames(peak.gene.links) <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))

# define resolutions to analyze
resolutions <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))

# vector to store link names (p.ex: MAP2_chr1:198197-198653)
names <- c()

# iterators
i = 1
z = 1

# iterate over GRN specific Peak-gene dfs
for (df in df_list){
  res <- resolutions[i]
  # iterate over the links
  for (j in 1:nrow(df)){
    # extract the link name
    name <- paste(df[j,1], df[j,2], sep = "_")
    # if that link has already appeared in another GRN then add 1 to the corresponding cell 
    if (name %in% rownames(peak.gene.links)) {
      peak.gene.links[name, i] = 1
    # otherwise, just create a new row for that link and assign a 1 to the correpponding GRN column
    } else {
      binary <- rep(0, 19)
      binary[i] <- 1
      names(binary) <- resolutions
      peak.gene.links <- rbind(peak.gene.links, binary)
      rownames(peak.gene.links)[z] <- name
      z = z + 1
    }
    binary <- c()
    name <- ""
  }
  i = i+1
}

# assign again colnames (somehow they are removed)
colnames(peak.gene.links) <- resolutions
head(peak.gene.links)

# upset plot
upset(peak.gene.links, sets = as.character(resolutions[1:length(resolutions)]), keep.order = TRUE, nintersects = NA)

# upset plot for links shared in at least 15 resolutions 
peak.gene.links.min15 <- peak.gene.links %>% 
  dplyr::filter(rowSums(.) >= 15)
colnames(peak.gene.links.min15) <- as.character(resolutions)

upset(peak.gene.links.min15, sets = as.character(resolutions[1:length(resolutions)]), keep.order = TRUE, nintersects = NA)

# inspect genes and peaks
peak.gene.df <- str_split_fixed(rownames(peak.gene.links.min15), "_", 2)
peak.gene.df[,1]
table(peak.gene.df[,1])
table(peak.gene.df[,2])

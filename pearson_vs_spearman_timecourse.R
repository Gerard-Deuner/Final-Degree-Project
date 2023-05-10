#####################################
# TIMECOURSE-SPEARMAN eGRNs ANALYSES #
#####################################

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

#############################################
# ANALYSIS ON DIFFERENT CLUSTER RESOLUTIONS #
#############################################

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
  
  # vector to store the GRN metrics
  metrics <- c(res)
  
  # retrieve metrics from the GRN
  GRN.uniqueTFs <- GRN@connections$all.filtered$`0` %>% pull("TF.name") %>% unique() %>% as.data.frame()
  GRN.uniqueTFs.len <- GRN.uniqueTFs %>% nrow()
  metrics <- c(metrics, GRN.uniqueTFs.len)
  GRN.uniquePEAKs <- GRN@connections$all.filtered$`0` %>% pull("peak.ID") %>% unique() %>% as.data.frame()
  GRN.uniquePEAKs.len <- GRN.uniquePEAKs %>% nrow()
  metrics <- c(metrics, GRN.uniquePEAKs.len)
  GRN.uniqueGENEs <- GRN@connections$all.filtered$`0` %>% pull("gene.name") %>% unique() %>% as.data.frame()
  GRN.uniqueGENEs.len <- GRN.uniqueGENEs %>% nrow()
  metrics <- c(metrics, GRN.uniqueGENEs.len)
  GRN.TF_peak.links <- GRN@connections$all.filtered$`0` %>% as.data.frame() %>% dplyr::select(c("TF.name", "peak.ID")) %>% distinct()
  GRN.TF_peak.links.len <- GRN.TF_peak.links %>% nrow()
  metrics <- c(metrics, GRN.TF_peak.links.len)
  GRN.peak_gene.links <- GRN@connections$all.filtered$`0` %>% as.data.frame() %>% dplyr::select(c("gene.name", "peak.ID")) %>% distinct()
  GRN.peak_gene.links.len <- GRN.peak_gene.links%>% nrow()
  metrics <- c(metrics, GRN.peak_gene.links.len)
  GRN.TF_peak_gene.links <- GRN@connections$all.filtered$`0` %>% as.data.frame() %>% dplyr::select(c("TF.name", "peak.ID", "gene.name")) %>% distinct()
  GRN.TF_peak_gene.links.len <- GRN.TF_peak_gene.links %>% nrow()
  metrics <- c(metrics, GRN.TF_peak_gene.links.len)
  
  # store the GRN metrics in the dataframe containing all GRNs metrics
  names(metrics) <- metrics.names
  GRN.metrics <- GRN.metrics %>% rbind(metrics)
  
  # store the TF-PEAK-GENE links in independent dataframes
  assign(paste("TF.PEAK.GENE", res, sep = "."), GRN.TF_peak_gene.links)
  assign(paste("PEAK.GENE", res, sep = "."), GRN.peak_gene.links)
  assign(paste("TF.PEAK", res, sep = "."), GRN.TF_peak.links)
  assign(paste("TFs", res, sep = "."), GRN.uniqueTFs)
  assign(paste("PEAKs", res, sep = "."), GRN.uniquePEAKs)
  assign(paste("GENEs", res, sep = "."), GRN.uniqueGENEs)
  
  # remove GRN object from memory
  rm(GRN)
  # reset metrics vector
  metrics <- c()
}

# rename the columns because somehow they do not show up
colnames(GRN.metrics) <- metrics.names
GRN.metrics$res <- as.factor(GRN.metrics$res)

colnames(GRN.metrics)

######################################################
# PLOT SOME METRICS OF DIFFERENT CLUSTER RESOLUTIONS #
######################################################

p1 <- ggplot(GRN.metrics, aes(res, uniqueTFs, fill = res)) + 
  geom_bar(stat = "identity") +
  labs(title = "Number of unique TFs per eGRN", x = "Cluster Resolution") +
  scale_fill_viridis_d(option = "inferno") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

p2 <- ggplot(GRN.metrics, aes(res, uniquePEAKs, fill = res)) + 
  geom_bar(stat = "identity") +
  labs(title = "Number of unique PEAKs per eGRN", x = "Cluster Resolution") +
  scale_fill_viridis_d(option = "inferno") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

p3 <- ggplot(GRN.metrics, aes(res, uniqueGENEs, fill = res)) + 
  geom_bar(stat = "identity") +
  labs(title = "Number of unique GENEs per eGRN", x = "Cluster Resolution") +
  scale_fill_viridis_d(option = "inferno") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

p4 <- ggplot(GRN.metrics, aes(res, TF.PEAK.links, fill = res)) + 
  geom_bar(stat = "identity") +
  labs(title = "Number of TF-PEAK links per eGRN", x = "Cluster Resolution") +
  scale_fill_viridis_d(option = "inferno") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

p5 <- ggplot(GRN.metrics, aes(res, PEAK.GENE.links, fill = res)) + 
  geom_bar(stat = "identity") +
  labs(title = "Number of PEAK-GENE links per eGRN", x = "Cluster Resolution") +
  scale_fill_viridis_d(option = "inferno") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

p6 <- ggplot(GRN.metrics, aes(res, TF.PEAK.GENE.links, fill = res)) + 
  geom_bar(stat = "identity") +
  labs(title = "Number of TF-PEAK-GENE links per eGRN", x = "Cluster Resolution") +
  scale_fill_viridis_d(option = "inferno") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

ggarrange(p1, p2, p3, p4, p5, p6 + rremove("x.text"), 
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 3)

###########################
# OVERLAPPING CONNECTIONS #
###########################

######################
# TF-peak-gene links #  
######################

# create binary df with all links and resolutions as columns
df_list <- list(TF.PEAK.GENE.0.1, TF.PEAK.GENE.0.25, TF.PEAK.GENE.0.5, TF.PEAK.GENE.0.75, TF.PEAK.GENE.1, TF.PEAK.GENE.2, TF.PEAK.GENE.3, TF.PEAK.GENE.4,
                TF.PEAK.GENE.5, TF.PEAK.GENE.6, TF.PEAK.GENE.7, TF.PEAK.GENE.8, TF.PEAK.GENE.9, TF.PEAK.GENE.10, TF.PEAK.GENE.12, TF.PEAK.GENE.14,
                TF.PEAK.GENE.16, TF.PEAK.GENE.18, TF.PEAK.GENE.20)


all.links <- data.frame(matrix(ncol = 19, nrow = 0))
colnames(all.links) <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))

resolutions <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))
names <- c()

i = 1
z = 1
for (df in df_list){
  res <- resolutions[i]
  for (j in 1:nrow(df)){
    name <- paste(df[j,1], df[j,2], sep = "_")
    name <- paste(name, df[j,3], sep = "_")
    if (name %in% rownames(all.links)) {
      all.links[name, i] = 1
    } else {
      binary <- rep(0, 19)
      binary[i] <- 1
      names(binary) <- resolutions
      all.links <- rbind(all.links, binary)
      rownames(all.links)[z] <- name
      z = z + 1
    }
    binary <- c()
    name <- ""
  }
  i = i+1
}

colnames(all.links) <- resolutions
all.links

# upset plot
all.links.min3 <- all.links %>% dplyr::filter(rowSums(.) > 15)

all.links.min3 %>% colnames() %>% class() == resolutions %>%  class()

colnames(all.links.min3) <- as.character(resolutions)

upset(all.links, sets = as.character(resolutions[1:length(resolutions)]), keep.order = TRUE, nintersects = NA)

upset(all.links.min3, sets = as.character(resolutions[1:length(resolutions)]), keep.order = TRUE, nintersects = NA)

# Take a look at links shared min among 15 resolutions
all.links.min3


# overlap analysis on optimal resolutions
upset(all.links[,5:14], sets = as.character(resolutions[5:14]), keep.order = TRUE, nintersects = NA)

# links present in all of them
upset(all.links[,5:14] %>% dplyr::filter(rowSums(.) > 7), sets = as.character(resolutions[5:14]), keep.order = TRUE, nintersects = NA)

shared.links <- all.links[,5:14] %>% dplyr::filter(rowSums(.) == 10)

shared.links.df <- str_split_fixed(rownames(shared.links), "_", 3)

# check TFs
table(shared.links.df[,1])

# check genes
shared.links.df[,3]


# check peaks
shared.links.df[,2]





# add resolution id column in all tf-peak-gene dfs
df.lens <- c()
for (df in df_list){
  df.lens <- c(df.lens, nrow(df)) 
}

# Merge all TF-peak-gene interactions
merged.all.links <- do.call("rbind", df_list)
merged.all.links <- merged.all.links %>% cbind(res = rep(resolutions, times = df.lens))
merged.all.links

# Reduce TF.name to just the name
new.TFs <- c()
for (TF in merged.all.links$TF.name){
  new.TFs <- c(new.TFs, strsplit(TF, "[.]")[[1]][1])
}
merged.all.links$TF <- new.TFs

ggplot(merged.all.links %>% dplyr::filter(res %in% resolutions[10:19]), aes(TF, fill = as.factor(res))) + 
  geom_bar(stat = "count", position = "stack") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



###################
# Peak-Gene links #
###################

# create binary df with all peak-genes links and resolutions as columns
df_list2 <- list(PEAK.GENE.0.1, PEAK.GENE.0.25, PEAK.GENE.0.5, PEAK.GENE.0.75, PEAK.GENE.1, PEAK.GENE.2, PEAK.GENE.3, PEAK.GENE.4,
                 PEAK.GENE.5, PEAK.GENE.6, PEAK.GENE.7, PEAK.GENE.8, PEAK.GENE.9, PEAK.GENE.10, PEAK.GENE.12, PEAK.GENE.14,
                 PEAK.GENE.16, PEAK.GENE.18, PEAK.GENE.20)

peak.gene.links <- data.frame(matrix(ncol = 19, nrow = 0))
colnames(peak.gene.links) <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))
peak.gene.links
resolutions <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))
names <- c()

i = 1
z = 1
for (df in df_list2){
  res <- resolutions[i]
  for (j in 1:nrow(df)){
    name <- paste(df[j,1], df[j,2], sep = "_")
    if (name %in% rownames(peak.gene.links)) {
      peak.gene.links[name, i] = 1
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

colnames(peak.gene.links) <- resolutions
peak.gene.links

# upset plot
upset(peak.gene.links, sets = as.character(resolutions), keep.order = TRUE, nintersects = NA)


peak.gene.links.min3 <- peak.gene.links %>% dplyr::filter(rowSums(.) >= 15)

peak.gene.links.min3 %>% colnames() %>% class() == resolutions %>%  class()

colnames(peak.gene.links.min3) <- as.character(resolutions)

upset(peak.gene.links.links, sets = as.character(resolutions[1:length(resolutions)]), keep.order = TRUE, nintersects = NA)

upset(peak.gene.links.min3, sets = as.character(resolutions[1:length(resolutions)]), keep.order = TRUE, nintersects = NA)
peak.gene.links.min3

library(stringr)
peak.gene.df <- str_split_fixed(rownames(peak.gene.links.min3), "_", 2)
peak.gene.df[,1]
table(peak.gene.df[,1])
table(peak.gene.df[,2])

###################
# TF-Peak links #
###################

# create binary df with all peak-genes links and resolutions as columns
df_list3 <- list(TF.PEAK.0.1, TF.PEAK.0.25, TF.PEAK.0.5, TF.PEAK.0.75, TF.PEAK.1, TF.PEAK.2, TF.PEAK.3, TF.PEAK.4,
                 TF.PEAK.5, TF.PEAK.6, TF.PEAK.7, TF.PEAK.8, TF.PEAK.9, TF.PEAK.10, TF.PEAK.12, TF.PEAK.14,
                 TF.PEAK.16, TF.PEAK.18, TF.PEAK.20)

tf.peak.links <- data.frame(matrix(ncol = 19, nrow = 0))
colnames(tf.peak.links) <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))
tf.peak.links
resolutions <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))
names <- c()

i = 1
z = 1
for (df in df_list3){
  res <- resolutions[i]
  for (j in 1:nrow(df)){
    name <- paste(df[j,1], df[j,2], sep = "_")
    if (name %in% rownames(tf.peak.links)) {
      tf.peak.links[name, i] = 1
    } else {
      binary <- rep(0, 19)
      binary[i] <- 1
      names(binary) <- resolutions
      tf.peak.links <- rbind(tf.peak.links, binary)
      rownames(tf.peak.links)[z] <- name
      z = z + 1
    }
    binary <- c()
    name <- ""
  }
  i = i+1
}

colnames(tf.peak.links) <- resolutions
tf.peak.links

tf.peak.links.min3 <- tf.peak.links %>% dplyr::filter(rowSums(.) >= 15)

tf.peak.links.min3 %>% colnames() %>% class() == resolutions %>%  class()

colnames(tf.peak.links.min3) <- as.character(resolutions)

upset(tf.peak.links, sets = as.character(resolutions[1:length(resolutions)]), keep.order = TRUE, nintersects = NA)

upset(tf.peak.links.min3, sets = as.character(resolutions[1:length(resolutions)]), keep.order = TRUE, nintersects = NA)
tf.peak.links.min3

# upset plot
upset(tf.peak.links, sets = as.character(resolutions), keep.order = TRUE, nintersects = NA)

tf.peak.df <- str_split_fixed(rownames(tf.peak.links.min3), "_", 2)
tf.peak.df[,1]
table(tf.peak.df[,1])
table(tf.peak.df[,2])
tf.peak.df[,2]

##############
# unique TFs #
##############

# create binary df with all peak-genes links and resolutions as columns
df_list4 <- list(TFs.0.1, TFs.0.25, TFs.0.5, TFs.0.75, TFs.1, TFs.2, TFs.3, TFs.4,
                 TFs.5, TFs.6, TFs.7, TFs.8, TFs.9, TFs.10, TFs.12, TFs.14,
                 TFs.16, TFs.18, TFs.20)

tfs <- data.frame(matrix(ncol = 19, nrow = 0))
colnames(tfs) <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))
tfs
resolutions <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))
names <- c()
tfs

i = 1
z = 1
for (df in df_list4){
  res <- resolutions[i]
  for (j in 1:nrow(df)){
    name <- paste0(df[j,1], df[j,2])
    if (name %in% rownames(tfs)) {
      tfs[name, i] = 1
    } else {
      binary <- rep(0, 19)
      binary[i] <- 1
      names(binary) <- resolutions
      tfs <- rbind(tfs, binary)
      rownames(tfs)[z] <- name
      z = z + 1
    }
    binary <- c()
    name <- ""
  }
  i = i+1
}


colnames(tfs) <- resolutions
tfs

# upset plot
upset(tfs, sets = as.character(resolutions), keep.order = TRUE, nintersects = NA)



library(stringr)
peak.gene.df <- str_split_fixed(rownames(peak.gene.links.min3), "_", 2)
peak.gene.df[,1]
table(peak.gene.df[,1])
table(peak.gene.df[,2])

################
# unique PEAKs #
################

# create binary df with all peak-genes links and resolutions as columns
df_list5 <- list(PEAKs.0.1, PEAKs.0.25, PEAKs.0.5, PEAKs.0.75, PEAKs.1, PEAKs.2, PEAKs.3, PEAKs.4,
                 PEAKs.5, PEAKs.6, PEAKs.7, PEAKs.8, PEAKs.9, PEAKs.10, PEAKs.12, PEAKs.14,
                 PEAKs.16, PEAKs.18, PEAKs.20)

peaks <- data.frame(matrix(ncol = 19, nrow = 0))
colnames(peaks) <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))
peaks
resolutions <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))
names <- c()
peaks

i = 1
z = 1
for (df in df_list5){
  res <- resolutions[i]
  for (j in 1:nrow(df)){
    name <- paste0(df[j,1], df[j,2])
    if (name %in% rownames(peaks)) {
      peaks[name, i] = 1
    } else {
      binary <- rep(0, 19)
      binary[i] <- 1
      names(binary) <- resolutions
      peaks <- rbind(peaks, binary)
      rownames(peaks)[z] <- name
      z = z + 1
    }
    binary <- c()
    name <- ""
  }
  i = i+1
}


colnames(peaks) <- resolutions
peaks

# upset plot
upset(peaks, sets = as.character(resolutions), keep.order = TRUE, nintersects = NA)

################
# unique GENEs #
################

# create binary df with all peak-genes links and resolutions as columns
df_list6 <- list(GENEs.0.1, GENEs.0.25, GENEs.0.5, GENEs.0.75, GENEs.1, GENEs.2, GENEs.3, GENEs.4,
                 GENEs.5, GENEs.6, GENEs.7, GENEs.8, GENEs.9, GENEs.10, GENEs.12, GENEs.14,
                 GENEs.16, GENEs.18, GENEs.20)

genes <- data.frame(matrix(ncol = 19, nrow = 0))
colnames(genes) <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))
genes
resolutions <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))
names <- c()
genes

i = 1
z = 1
for (df in df_list6){
  res <- resolutions[i]
  for (j in 1:nrow(df)){
    name <- paste0(df[j,1], df[j,2])
    if (name %in% rownames(genes)) {
      genes[name, i] = 1
    } else {
      binary <- rep(0, 19)
      binary[i] <- 1
      names(binary) <- resolutions
      genes <- rbind(genes, binary)
      rownames(genes)[z] <- name
      z = z + 1
    }
    binary <- c()
    name <- ""
  }
  i = i+1
}


colnames(genes) <- resolutions
genes

# upset plot
upset(genes, sets = as.character(resolutions), keep.order = TRUE, nintersects = NA)


###########################################
# VALIDATION USING HI-C DATA FROM NEURONS #
###########################################

# get pcHiC links
links.dir <- "/g/scb/zaugg/deuner/GRaNIE/validationdata/timecourse_pchic_links.csv"
tc.links <- read.csv(links.dir)
head(tc.links)

# Set up path of the cluster resolutions eGRNs
path <- "/g/scb2/zaugg/deuner/GRaNIE/outputdata/timecourse_batch_mode/Batch_Mode_Outputs/"

# set vector with the resolutions
resolutions <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))

# set colors for ROC curves
colors <- c("darkred","deeppink4" , "deeppink", "pink","violet", 
                     "purple", "darkmagenta",  "darkblue", "blue", "darkcyan", 
                     "lightblue", "turquoise", "lightgreen", "green", "darkolivegreen1",
                     "yellow", "darkgoldenrod2",  "chocolate", "red")
                     
# set vecor of auc positioning in y axis
aucy <- seq(10, 110, by = 5)

# iterate over GRNs
j = 0

# boolean for first plot
first = TRUE

# vectors to store TPs and FPs
TP.vec <- c()
FP.vec <- c()

TP.all <- data.frame()

for (res in resolutions){
  # set index
  j <- j + 1
  
  # set up GRN directory
  GRN_dir <- paste0(path, "output_pseudobulk_clusterRes", res, "_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/")
  # read GRN object
  GRN <- qread(paste0(GRN_dir, "GRN.qs"))
  
  # get gene-peak connections from GRN
  GRN_links <- GRN@connections$all.filtered$`0` %>% as.data.frame() %>% dplyr::select(peak.ID, gene.name, peak_gene.r, peak_gene.p_raw)
  
  # adapt peaks format to tc.links format
  form.peak <- rep("", nrow(GRN_links))
  for (i in 1:nrow(GRN_links)){
    peak <- as.character(GRN_links$peak.ID[i])
    split_peak <- strsplit(peak, split = ":")
    new_peak <- paste(split_peak[[1]][1], split_peak[[1]][2], sep = "-")
    form.peak[i] <- new_peak
  }
  GRN_links$id <- seq(1:nrow(GRN_links))
  GRN_links$peak.ID <- form.peak 
  names(GRN_links)[1:2] <- c("peak", "gene")
  
  # Binary classification
  # get rows that match in GRN links and HiC links by gene name and peak id
  TP <- inner_join(GRN_links, tc.links, by = c("peak", "gene"), multiple = "all")
  TP$res <- rep(resolutions[j], nrow(TP))
  TP.all <- rbind(TP.all, TP)
  # set 1 if they are in HiC data and 0 if not
  bc <- ifelse(GRN_links$id %in% TP$id, yes = 1, no = 0)
  # add column to GRN_links with that information
  GRN_links$bc <- bc
  # create an confusion matrix 
  conf.m <- table(bc)
  
  # store GRN links data in a variable
  assign(paste("GRN_links", res, sep = "."), GRN_links)
  
  # Compute ROC curve
  par(pty = "s")
  
  # print metrics
  print("")
  print(paste0("Resolution: ", res))
  print(paste0("GRN links: ", nrow(GRN_links)))
  print(paste0("TP: ", nrow(TP)))
  FP <- nrow(GRN_links) - nrow(TP)
  print(paste0("FP: ", FP))
  print("")
  
  TP.vec <- c(TP.vec, nrow(TP))
  FP.vec <- c(FP.vec, FP)
  
  if ((first == TRUE) & (sum(bc) > 0)){
    roc(GRN_links$bc, GRN_links$peak_gene.r, plot = TRUE, legacy.axes = TRUE, percent = TRUE,
        xlab="False Positive Percentage | 1 - Specificity", ylab="True Positive Percentage | Sensitivity",
        col = colors[j], lwd=2, print.auc = TRUE, main = "ROC", add = FALSE, print.auc.x = 115, print.auc.y = aucy[j])
    legend("bottomright",  legend = resolutions, col = colors, ncol = 2)
    first = FALSE
  }
  if ((first == FALSE) & (sum(bc) > 0)){ 
    roc(GRN_links$bc, GRN_links$peak_gene.r, plot = TRUE, legacy.axes = TRUE, percent = TRUE,
        xlab="False Positive Percentage | 1 - Specificity", ylab="True Positive Percentage | Sensitivity",
        col = colors[j], lwd=2, print.auc = TRUE, main = "ROC", add = TRUE, print.auc.x = 115, print.auc.y = aucy[j])
    legend("bottomright", legend = resolutions, col = colors, lwd = 4, ncol = 2)
  }
}


# barplots of TP vs FP
TPFP.df <- data.frame(TPFP = c(TP.vec, FP.vec), cat = c(rep("TP", 19), rep("FP", 19)), res = c(resolutions <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2)), resolutions <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))))
ggplot(TPFP.df, aes(x = as.factor(res), y = TPFP, fill = cat)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  theme_classic() + 
  scale_fill_manual(values = c("#E3B448", "#3A6B35"))

TP.all

# Check if the TPs are shared 
TP.peak.gene <- TP.all %>% dplyr::select(c("peak", "gene", "res"))

# count occurrence of each link
TP.peak.gene.counts <- TP.peak.gene %>% dplyr::select(c("peak", "gene", "res")) %>% group_by(peak, gene) %>% count 

# create a column for the pasted peak-gene info
TP.peak.gene$link <- paste(TP.peak.gene$peak, TP.peak.gene$gene, sep = "_")

# join counts info to the dataframe
TP.peak.gene <- inner_join(TP.peak.gene, TP.peak.gene.counts, by = c("peak", "gene"))

# check whether the link is shared or not
TP.peak.gene$shared <- ifelse(TP.peak.gene$n > 1, TRUE, FALSE)

# plot count of shared links vs non-shared links
p1 <- ggplot(TP.peak.gene %>% distinct(peak, gene, .keep_all = TRUE), aes(shared, fill = shared)) + 
  geom_bar(stat = "count") + 
  theme_classic() + 
  scale_fill_manual(values = c("#E3B448", "#3A6B35")) + 
  labs(title = "Are the links shared? (unique links)") + 
  geom_text(stat='count', aes(label=..count..), vjust=-1) + 
  theme(legend.position = "none")

p2 <- ggplot(TP.peak.gene, aes(shared, fill = shared)) + 
  geom_bar(stat = "count") + 
  theme_classic() + 
  scale_fill_manual(values = c("#E3B448", "#3A6B35")) + 
  labs(title = "Are the links shared? (all links)") + 
  geom_text(stat='count', aes(label=..count..), vjust=-1)

p1 + p2

# plot duplicated vs unique links per cluster resolution
ggplot(TP.peak.gene, aes(as.factor(res), fill = shared)) + 
  geom_bar(stat = "count", position = "dodge") + 
  theme_classic() +
  scale_fill_manual(values = c("#E3B448", "#3A6B35")) + 
  labs(title = "Are the links shared?", x = "Cluster Resolutions")


# plot most frequent links across resolutions
TP.peak.gene.counts$link <- paste(TP.peak.gene.counts$peak, TP.peak.gene.counts$gene, sep = "_")
TP.peak.gene.counts$shared <- ifelse(TP.peak.gene.counts$n > 1, TRUE, FALSE)

ggplot(TP.peak.gene.counts %>% arrange(desc(n)) %>% head(10), aes(x = link, y = n, fill = link)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") + 
  scale_fill_viridis_d() + 
  labs(y = "count", title = "Most frequent links across resolutions", caption = "Top 10 Most Frequent Links") 


# plot only duplicated links and their corresponding cluster resolutions
library(RColorBrewer)
palette3_info <- brewer.pal.info[brewer.pal.info$category == "qual", ]  # Extract color info
palette3_all <- unlist(mapply(brewer.pal,                     # Create vector with all colors
                              palette3_info$maxcolors,
                              rownames(palette3_info)))
palette3_all 
palette3 <- sample(palette3_all, 16)                    # Sample colors
palette3

ggplot(TP.peak.gene %>% dplyr::filter(shared == TRUE) , aes(x = link, fill = as.factor(res))) + 
  geom_bar(stat = "count") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y = "count", title = "To which eGRNs do the shared links belong?", fill = "Resolutions", caption = "26 Different Links") + 
  scale_fill_manual(values = palette3)




#####################################################################################
# CHECK IF MOST OVERLAPPING LINKS IN DIFFERENT CLUSTERS ARE PRESENT IN PCHI-C LINKS #
#####################################################################################

for (deg.th in 1:14){
  # get highly represented links across resolutions
  high.rep.pg.links <-  all.links[rowSums(all.links[,1:ncol(all.links)]) > deg.th,]
  
  # overview of dataframes
  ## links from different cluster resolutions
  head(merged.all.links) # all links merged
  head(high.rep.pg.links)  # links - binary table indicating in which cluster resolution are inferred
  ## pcHiC links
  head(TP.all) # links of timecourse dataset found in pcHi-C interactions , chr1, start1, and start2 are the peak coordinates
  
  # create vector for each high degree link to indicate if this link is found in the pchic links or not
  high.rep.links.pcHiC <- high.rep.pg.links
  pcHiC <- rep(0, nrow(high.rep.links.pcHiC))
  high.rep.links.pcHiC$pcHiC <- pcHiC
  
  # format TP.all peak-gene links
  pchic.links <- c()
  for (i in 1:nrow(TP.all)){
    split.link <- strsplit(TP.all[i, "peak"], split = "-")
    pchic.link <- paste(split.link[[1]][1], split.link[[1]][2], sep = ":")
    pchic.link <- paste(pchic.link, split.link[[1]][3], sep = "-")
    pchic.link <- paste(pchic.link, TP.all[i, "gene"], sep = ";")
    pchic.links <- c(pchic.links, pchic.link)
  }
  pchic.links
  
  # find links that are in pchic data
  links.in.pchic <- c()
  for (i in 1:nrow(high.rep.pg.links)){
    split.link <- strsplit(rownames(high.rep.pg.links)[i], split = "_")
    rep.link <- paste(split.link[[1]][2], split.link[[1]][3], sep = ";")
    if (rep.link %in% pchic.links) {
      high.rep.links.pcHiC[i, "pcHiC"] = 1
      links.in.pchic <- c(links.in.pchic, rep.link) 
    }
  }
  
  print(links.in.pchic)
  
  
  
  
  p <- ggplot(high.rep.links.pcHiC, aes(pcHiC, fill = as.factor(pcHiC))) + 
    geom_bar(stat = "count") + 
    labs(caption = paste("Degree Threshold:", deg.th), title = "High degree links mapping with pcHi-C links", fill = "pcHi-C") + 
    theme_bw() + 
    annotate("text", x = 1, y = max(table(high.rep.links.pcHiC$pcHiC))/2, label = as.character(length(links.in.pchic)), size = 7) + 
    annotate("text", x = 0, y = max(table(high.rep.links.pcHiC$pcHiC))/2, label = as.character(nrow(high.rep.links.pcHiC)-length(links.in.pchic)), size = 7) + 
    scale_fill_manual(values = c("#E3B448", "#3A6B35"))
  
  assign(paste0("p", deg.th), p)
}
p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12

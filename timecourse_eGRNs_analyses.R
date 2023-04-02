##############################
# TIMECOURSE GRaNIE ANALYSES #
##############################

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

#%%%%%%%%%%%%%%%%%%%%%#
# PEARSON CORRELATION #
#%%%%%%%%%%%%%%%%%%%%%#

#############################################
# ANALYSIS ON DIFFERENT CLUSTER RESOLUTIONS #
#############################################

# Set up path of the cluster resolutions eGRNs
path <- "/g/scb/zaugg/zaugg_shared/data/10xMultiome/Daria/timecourse/GRaNIE/"

# set vector with the resolutions
resolutions <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))

# creat data.frame to store each GRN metrics
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
all.links.min3 <- all.links %>% dplyr::filter(rowSums(.) > 3)

all.links.min3 %>% colnames() %>% class() == resolutions %>%  class()

colnames(all.links.min3) <- as.character(resolutions)

upset(all.links, sets = as.character(resolutions[1:length(resolutions)]), keep.order = TRUE, nintersects = NA)

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

# upset plot
upset(tf.peak.links, sets = as.character(resolutions), keep.order = TRUE, nintersects = NA)


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
path <- "/g/scb/zaugg/zaugg_shared/data/10xMultiome/Daria/timecourse/GRaNIE/"

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
  labs(y = "count", title = "To which eGRNs do the shared links belong?", fill = "Resolutions", caption = "35 Different Links") + 
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
p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14

#%%%%%%%%%%%%%%%%%%%%%%#
# SPEARMAN CORRELATION #
#%%%%%%%%%%%%%%%%%%%%%%#


#%%%%%%%%%%%%%%%%%%%%%%#
# PEARSON VS. SPEARMAN #
#%%%%%%%%%%%%%%%%%%%%%%#

# Cell Type Pseudobulking - WNN - res: 0.5 - 16 clusters

# PEARSON
# load eGRN
GRN.p <- qread("/g/scb/zaugg/deuner/GRaNIE/outputdata/timecourse_celltype_wnn/output_pseudobulk_celltype_wnn_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/GRN.qs")
GRN.p

# build eGRN
GRN.p = build_eGRN_graph(GRN.p, forceRerun = TRUE)
GRN.p = visualizeGRN(GRN.p, plotAsPDF = FALSE)

# SPEARMAN
GRN.s <- qread("/g/scb/zaugg/deuner/GRaNIE/outputdata/timecourse_celltype_spearman_wnn/output_pseudobulk_celltype_wnn_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/GRN.qs")
GRN.s

# build eGRN
GRN.s = build_eGRN_graph(GRN.s, forceRerun = TRUE)
GRN.s = visualizeGRN(GRN.s, plotAsPDF = FALSE)

GRN.p@connections$all.filtered #more connections
GRN.s@connections$all.filtered 

# Set up the main directory
path = "/g/scb/zaugg/deuner/GRaNIE"
# Load the seurat object
seuratFile = "timecourse.pp.seuratObject.qs"
timecourse.s <- qread(paste0(path,"/tmp/",seuratFile))

## PEARSON
# TF-peak-gene connections here
GRN.p@connections$all.filtered[["0"]]
connections.p <- as.data.frame(GRN.p@connections$all.filtered[["0"]])

# SCT counts here
GRN.p@data[["RNA"]][["counts"]]
counts.p <- as.data.frame(GRN.p@data[["RNA"]][["counts"]])
counts.p$gene.ENSEMBL <- rownames(counts.p)
counts.p

# norm ATAC counts here
GRN.p@data[["peaks"]][["counts"]]
peaks.p <- as.data.frame(GRN.p@data[["peaks"]][["counts"]])
peaks.p$peak.ID <- rownames(peaks.p)
peaks.p

# connections + counts
joined.p <- connections.p %>% inner_join(counts.p, by = "gene.ENSEMBL") %>% inner_join(peaks.p, by = "peak.ID")


## SPEARMAN
# peak-gene connections here
GRN.s@connections$all.filtered[["0"]]
connections.s <- as.data.frame(GRN.s@connections$all.filtered[["0"]])

# SCT counts here
GRN.s@data[["RNA"]][["counts"]]
counts.s <- as.data.frame(GRN.s@data[["RNA"]][["counts"]])
counts.s$gene.ENSEMBL <- rownames(counts.s)
counts.s

# norm ATAC counts here
GRN.s@data[["peaks"]][["counts"]]
peaks.s <- as.data.frame(GRN.s@data[["peaks"]][["counts"]])
peaks.s$peak.ID <- rownames(peaks.s)
peaks.s

# connections + counts
joined.s <- connections.s %>% inner_join(counts.s, by = "gene.ENSEMBL") %>% inner_join(peaks.s, by = "peak.ID")

# check shared TF-peak-links
tf.peak.gene.p <- GRN.p@connections$all.filtered[["0"]]
tf.peak.gene.s <- GRN.s@connections$all.filtered[["0"]]
colnames(tf.peak.gene.p)

tf.shared.p.s <- inner_join(tf.peak.gene.p, tf.peak.gene.s, by = c("TF.ID")) 
gene.shared.p.s <- inner_join(tf.peak.gene.p, tf.peak.gene.s, by = c("gene.name"))
peak.shared.p.s <- inner_join(tf.peak.gene.p, tf.peak.gene.s, by = c("peak.ID"))

View(tf.peak.gene.p)
tf.peak.gene.s

tf.shared.p.s
gene.shared.p.s
peak.shared.p.s


# CHECK IF ASSUMPTIONS ARE HELD!
# Pearson:

# select top5 most significant peak-gene links
top5.p <- joined.p %>% top_n(5, peak_gene.r)
top5.p.1 <- top5.p[,1:32]
top5.p.2 <- top5.p[,c(1:16, 33:48)]
# long format
top5.p.long.1 <- reshape2::melt(top5.p.1 , id.vars = colnames(top5.p)[1:16], variable.name = "cluster.RNA", value.name = "value.RNA")
top5.p.long.2 <- reshape2::melt(top5.p.2, id.vars = colnames(top5.p)[1:16], variable.name = "cluster.ATAC", value.name = "value.ATAC")
top5.p.both <- cbind(top5.p.long.1, top5.p.long.2[,17:ncol(top5.p.long.2)])

# assign informative colors to the cell type clusters
cols <- c("#FF6600", "#FF0000", "#CC0033", "#99FFFF", "#0099FF", "#33FF33", "#FF00FF", "#990000", "#FF3300", "#FFCC66", "#660066", "#0000FF", "#FFFF33", "#FF6666", "#FF9933", "#66FF99")

ggplot(top5.p.both, aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Peak-gene link with highest correlation coefficient") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_label_repel(data=top5.p.both[top5.p.both$TF.name == "YY1",], aes(label=as.factor(TF.name), alpha=0.7), size=5, force=1.3, max.overlaps = 30) + 
  scale_color_manual(values = cols)

table(joined.p$TF.ID)


ggplot(top5.p.both, aes(x = value.RNA, y = value.ATAC, col = cluster.RNA, label=TF.name)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log 10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 10 most correlated peak-genes") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_text(hjust=0, vjust=0)

# set new ids
new.ids <- as.character(c(1:5))
top5.p.both$id <- new.ids
names(new.ids) <- levels(top5.p.both$id)


ggplot(top5.p.both, aes(x = log10(value.RNA), y = value.ATAC, col = id)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Id", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-genes") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5))


p1 <- ggplot(top5.p.both %>% dplyr::filter(id == "1"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "TTop 1 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p2 <-ggplot(top5.p.both %>% dplyr::filter(id == "2"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 2 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p3 <-  ggplot(top5.p.both %>% dplyr::filter(id == "3"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 3 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p4 <-ggplot(top5.p.both %>% dplyr::filter(id == "4"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 4 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p5 <-ggplot(top5.p.both %>% dplyr::filter(id == "5"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)

ggarrange(p1,p2,p3,p4,p5)
p1
p2

# Spearman:

# select top5 most significant peak-gene links
## add also the pearson coefficient as a columns
top5.s <- joined.s %>% top_n(5, peak_gene.r)
top5.s.1 <- top5.s[,1:32]
top5.s.2 <- top5.s[,c(1:16, 33:48)]
# long format
top5.s.long.1 <- reshape2::melt(top5.s.1 , id.vars = colnames(top5.s)[1:16], variable.name = "cluster.RNA", value.name = "value.RNA")
top5.s.long.2 <- reshape2::melt(top5.s.2, id.vars = colnames(top5.s)[1:16], variable.name = "cluster.ATAC", value.name = "value.ATAC")
top5.s.both <- cbind(top5.s.long.1, top5.s.long.2[,17:ncol(top5.s.long.2)])

# set new ids
new.ids <- as.character(c(1:5))
top5.s.both$id <- new.ids
names(new.ids) <- levels(top5.s.both$id)
new.ids

ggplot(top5.s.both, aes(x = log10(value.RNA), y = value.ATAC, col = id)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Id", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-genes") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5))

ggplot(top5.s.both, aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Peak-gene link with highest correlation coefficient") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)

p1 <- ggplot(top5.s.both %>% dplyr::filter(id == "1"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 1 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p2 <-ggplot(top5.s.both %>% dplyr::filter(id == "2"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 2 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p3 <-  ggplot(top5.s.both %>% dplyr::filter(id == "3"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 3 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p4 <-ggplot(top5.s.both %>% dplyr::filter(id == "4"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 4 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p5 <-ggplot(top5.s.both %>% dplyr::filter(id == "5"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)

ggarrange(p1,p2,p3,p4,p5)
p1
p4
p5

# check links overlaps

# merge links

# redoit but with all peak-gene links and putting pearson and spearman correlation coefficients in the same dataframe
# extract peak-genes from GRN ran with Spearman correlation 
pg.links <- GRN.s@connections$peak_genes$`0` # IMPORTANT LINE
colnames(pg.links)[4] <- "spearman.r" 
# add pearson correlation values
pg.links$pearson.r <- GRN.p@connections$peak_genes$`0`$peak_gene.r

# add id column
pg.links$id <- as.character(c(1:nrow(pg.links)))

# add gene expresion and peak accessibility counts for each connection
pg.joined.s <- pg.links %>% inner_join(counts.s, by = "gene.ENSEMBL") %>% inner_join(peaks.p, by = "peak.ID")

 ############
 # SPEARMAN # 
 ############

# select top5 most significant peak-gene links
## add also the pearson coefficient as a columns
pg.top5.s <- pg.joined.s %>% top_n(5, spearman.r)
pg.top5.s.1 <- pg.top5.s[,1:23]
pg.top5.s.2 <- pg.top5.s[,c(1:7, 24:39)]
# long format
pg.top5.s.long.1 <- reshape2::melt(pg.top5.s.1 , id.vars = colnames(pg.top5.s)[1:7], variable.name = "cluster.RNA", value.name = "value.RNA")
pg.top5.s.long.2 <- reshape2::melt(pg.top5.s.2, id.vars = colnames(pg.top5.s)[1:7], variable.name = "cluster.ATAC", value.name = "value.ATAC")
pg.top5.s.both <- cbind(pg.top5.s.long.1, pg.top5.s.long.2[,8:ncol(pg.top5.s.long.2)])
pg.top5.s.both

# set new ids
new.ids <- as.character(c(1:5))
names(new.ids) <- levels(pg.top5.s.both$id)
new.ids

ggplot(pg.top5.s.both, aes(x = log10(value.RNA), y = value.ATAC, col = id)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Id", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-genes") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5))


p1 <- ggplot(pg.top5.s.both %>% dplyr::filter(id == "316796"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 1 most correlated peak-gene link. Id: 316796") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p2 <-ggplot(pg.top5.s.both %>% dplyr::filter(id == "380465"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 2 most correlated peak-gene link. Id: 380465") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p3 <-  ggplot(pg.top5.s.both %>% dplyr::filter(id == "484736"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 3 most correlated peak-gene link. Id: 484736") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p4 <-ggplot(pg.top5.s.both %>% dplyr::filter(id == "714900"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 4 most correlated peak-gene link. Id: 714900") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p5 <-ggplot(pg.top5.s.both %>% dplyr::filter(id == "714906"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-gene link. Id: 714906") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)

ggarrange(p1,p2,p3,p4,p5)
p1

 ###########
 # PEARSON #
 ###########

# select top5 most significant peak-gene links
## add also the pearson coefficient as a columns
pg.top5.s <- pg.joined.s %>% top_n(5, pearson.r)
pg.top5.s.1 <- pg.top5.s[,1:23]
pg.top5.s.2 <- pg.top5.s[,c(1:7, 24:39)]
# long format
pg.top5.s.long.1 <- reshape2::melt(pg.top5.s.1 , id.vars = colnames(pg.top5.s)[1:7], variable.name = "cluster.RNA", value.name = "value.RNA")
pg.top5.s.long.2 <- reshape2::melt(pg.top5.s.2, id.vars = colnames(pg.top5.s)[1:7], variable.name = "cluster.ATAC", value.name = "value.ATAC")
pg.top5.s.both <- cbind(pg.top5.s.long.1, pg.top5.s.long.2[,8:ncol(pg.top5.s.long.2)])
pg.top5.s.both

# set new ids
new.ids <- as.character(c(1:5))
names(new.ids) <- levels(pg.top5.s.both$id)
new.ids

ggplot(pg.top5.s.both, aes(x = log10(value.RNA), y = value.ATAC, col = id)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Id", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-genes") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5))


p1 <- ggplot(pg.top5.s.both %>% dplyr::filter(id == "146686"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "TTop 1 most correlated peak-gene link. Id: 146686") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p2 <-ggplot(pg.top5.s.both %>% dplyr::filter(id == "18914"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 2 most correlated peak-gene link. Id: 18914") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p3 <-  ggplot(pg.top5.s.both %>% dplyr::filter(id == "441623"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 3 most correlated peak-gene link. Id: 441623") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p4 <-ggplot(pg.top5.s.both %>% dplyr::filter(id == "634342"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 4 most correlated peak-gene link. Id: 634342") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p5 <-ggplot(pg.top5.s.both %>% dplyr::filter(id == "757121"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-gene link. Id: 757121") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)

ggarrange(p1,p2,p3,p4,p5)
p1

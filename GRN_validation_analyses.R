#############################
# eGRNs VALIDATION ANALYSES #
#############################

## Experimental data used: pcHI-C, ChiP-Seq and eQTL

## input data format (one df per exp data)
## rows: links in GRN (ids) , columns: gene/TF, peak, resolution, validated, setting, dataset, corr.method,  
#         gene/TF         peak          resolution      validated            setting           dataset         corr.method     
# link1   FOX1    chr1:187183-187634        0.5           Yes           combined_spearman     combined          spearman  

# load libraries
library(dplyr)

# define datasets
datasets <- c("timecourse", "combined")
corr.methods <- c("pearson", "spearman")
val.data <- c("pcHi-C", "ChiP-seq", "eQTL")

# set validation data paths
pcHiC.dir <- "/g/scb/zaugg/deuner/valdata/pcHi-C/results/"
ChiP.dir <- "/g/scb/zaugg/deuner/valdata/ChiP-seq/results/"
eQTL.dir <- "/g/scb/zaugg/deuner/valdata/eQTL/results/"

# get pcHi-C data
pcHiC.df <- as.data.frame(matrix(nrow = 0, ncol = 5))
names(pcHiC.df) <- c("gene", "peak", "res", "validated", "setting", "dataset", "corr.method")

for (file in list.files(pcHiC.dir)){
  df <- read.csv(paste0(pcHiC.dir, file))
  df <- df %>%
    mutate(setting = rep(file, nrow(df)))
  df[c("dataset", "corr.method")] <- str_split_fixed(df$setting, "_", 2)
  pcHiC.df <- rbind(pcHiC.df, df)
}

# get ChiP-seq data
ChiP.df <- as.data.frame(matrix(nrow = 0, ncol = 5))
names(ChiP.df) <- c("TF", "peak", "res", "validated", "setting", "dataset", "corr.method")

for (file in list.files(ChiP.dir)){
  df <- read.csv(paste0(ChiP.dir, file))
  df <- df %>%
    mutate(setting = rep(file, nrow(df)))
  df[c("dataset", "corr.method")] <- str_split_fixed(df$setting, "_", 2)
  ChiP.df <- rbind(ChiP.df, df)
}

# get eQTL data
eQTL.df <- as.data.frame(matrix(nrow = 0, ncol = 5))
names(eQTL.df) <- c("gene", "peak", "res", "validated", "setting", "dataset", "corr.method")

for (file in list.files(eQTL.dir)){
  df <- read.csv(paste0(eQTL.dir, file))
  df <- df %>%
    mutate(setting = rep(file, nrow(df)))
  df[c("dataset", "corr.method")] <- str_split_fixed(df$setting, "_", 2)
  eQTL.df <- rbind(eQTL.df, df)
}

# Merge all 
val.df <- do.call("rbind", list(pcHiC.df, ChiP.df, eQTL.df))
val.df$validation <- rep(c("pcHi-C", "ChiP-seq", "eQTL"), times = c(nrow(pcHiC.df), nrow(ChiP.df), nrow(eQTL.df)))

#%%%%%%%#
# PLOTS #
#%%%%%%%#

# example data

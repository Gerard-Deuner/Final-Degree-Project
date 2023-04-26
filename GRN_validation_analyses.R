#############################
# eGRNs VALIDATION ANALYSES #
#############################

## Experimental data used: pcHI-C, ChiP-Seq and eQTL

#%%%%%%%#
# SETUP #
#%%%%%%%#

## input data format (one df per exp data)
## rows: links in GRN (ids) , columns: gene/TF, peak, resolution, validated, setting, dataset, corr.method,  
#         gene/TF         peak          resolution      validated            setting           dataset         corr.method     
# link1   FOX1    chr1:187183-187634        0.5           Yes           combined_spearman     combined          spearman  

# load libraries
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# define datasets
datasets <- c("timecourse", "combined")
corr.methods <- c("pearson", "spearman")
val.data <- c("pcHi-C", "ChiP-seq", "eQTL")

# set validation data paths
pcHiC.dir <- "/g/scb/zaugg/deuner/valdata/pcHi-C/validated_links/"
ChiP.dir <- "/g/scb/zaugg/deuner/valdata/ChiP-seq/validated_links/"
eQTL.dir <- "/g/scb/zaugg/deuner/valdata/eQTL/validated_links/"

# get pcHi-C data
pcHiC.df <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(pcHiC.df) <- c("gene", "peak", "res", "validated", "setting", "dataset", "corr.method")

for (file in list.files(pcHiC.dir)){
  print(file)
  df <- read.csv(paste0(pcHiC.dir, file), row.names = 1)
  dataset <- str_split(file, "_") %>% map(1) %>% as.character
  corr.method <- str_split(file, "_", 2) %>% map(2) %>% as.character %>% strsplit("[.]") %>% map(1) %>% as.character
  df <- df %>%
    tibble::add_column(setting = paste(dataset, corr.method, sep = "_"), 
                      dataset = dataset,
                      corr.method = corr.method)
  pcHiC.df <- rbind(pcHiC.df, df)
}

# get ChiP-seq data
ChiP.df <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(ChiP.df) <- c("TF", "peak", "res", "validated", "setting", "dataset", "corr.method")

for (file in list.files(ChiP.dir)){
  print(file)
  df <- read.csv(paste0(ChiP.dir, file), row.names = 1)
  dataset <- str_split(file, "_") %>% map(1) %>% as.character
  corr.method <- str_split(file, "_", 2) %>% map(2) %>% as.character %>% strsplit("[.]") %>% map(1) %>% as.character
  df <- df %>%
    tibble::add_column(setting = paste(dataset, corr.method, sep = "_"), 
                       dataset = dataset,
                       corr.method = corr.method)
  ChiP.df <- rbind(ChiP.df, df)
}

# get eQTL data
eQTL.df <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(eQTL.df) <- c("gene", "peak", "res", "validated", "setting", "dataset", "corr.method")

for (file in list.files(eQTL.dir)){
  print(file)
  df <- read.csv(paste0(eQTL.dir, file))
  dataset <- str_split(file, "_") %>% map(1) %>% as.character
  corr.method <- str_split(file, "_", 2) %>% map(2) %>% as.character %>% strsplit("[.]") %>% map(1) %>% as.character
  df <- df %>%
    tibble::add_column(setting = paste(dataset, corr.method, sep = "_"), 
                       dataset = dataset,
                       corr.method = corr.method)
  eQTL.df <- rbind(eQTL.df, df)
}

## Merge all 
# change TF column of ChiP-seq to gene so the names match 
names(ChiP.df)[1] <- "gene"
val.df <- do.call("rbind", list(pcHiC.df, ChiP.df)) # eQTL missing
val.df$validation <- rep(c("pcHi-C", "ChiP-seq"), times = c(nrow(pcHiC.df), nrow(ChiP.df))) # eQTL missing
val.df$resolution <- as.factor(val.df$resolution)
val.df$validated <- as.numeric(val.df$validated)

# add gene.ENSEMBL column
library(biomaRt)

ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)

ensembl.con <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)

annotLookup <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  #filter = "external_gene_name", # If I want the ENSEMBL IDs then it is ensembl_gene_id
  #values = val.df$gene,
  mart = ensembl.con
)

#%%%%%%%#
# PLOTS #
#%%%%%%%#




# TF recovery from ChiP-seq data
ggplot(val.df %>%                  # validated == 1 adds more occurrences, need to figure out why
         filter(validation == "ChiP-seq") %>% # want to measure frequency of recoveries #, grepl("nomicro", setting)
         distinct(gene, resolution, setting, .keep_all = TRUE) %>% # just consider the TF once per GRN, it if appears in the eGRN multiple times it has a recovery of 1
         dplyr::select(resolution, setting, validated) %>%  # remove useless columns
         group_by(setting, resolution) %>% # group the data by resolution and setting
         summarise(recoveredTFs = sum(validated), .groups = "drop"), # count number of TFs recovered per each resolution and setting
       aes(resolution, y = recoveredTFs , group = setting, col = setting)) + 
  geom_line(size = 1) +
  labs(y = "# TFs recovered", x = "Cluster Resolutions", title = "TF Recovery") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme_bw() 

# works well no!
val.df %>% 
  filter(validation == "ChiP-seq", resolution == 3, setting == "timecourse_pearson") %>% distinct(gene, resolution, setting, .keep_all = TRUE) %>%
  dplyr::select(validated) %>% table()

 # CHECK FIRST GROUPING BY SETTING LOOK THE SUMS AND THEN BY SETTING AND RES ANS SEE IF ALLL THE RES SUMS EQUAL THE SETTING SUM

ggplot(val.df %>% 
         filter(validation == "ChiP-seq") %>% # want to measure frequency of recoveries
         distinct(gene, resolution, setting, .keep_all = TRUE) %>% # just consider the TF once per GRN, it if appears in the eGRN multiple times it has a recovery of 1
         dplyr::select(resolution, setting, validated) %>%  # remove useless columns
         group_by(resolution, setting) %>% # group the data by resolution and setting
         summarise(recoveredTFs = sum(validated)), # count number of TFs recovered per each resolution and setting
       aes(resolution, y = recoveredTFs , group = setting, fill = setting)) + 
  geom_bar(stat = "identity") + 
  facet_grid("setting") + 
  labs(y = "# TFs recovered", x = "Cluster Resolutions", title = "TF Recovery") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.background = element_blank(),
        strip.text.y = element_blank())


# Peak-Gene links recovery from pcHi-C
ggplot(val.df %>%
         filter(validation == "pcHi-C"))

# Peaks recovery from pcHi-C

# Genes  recovery from pcHi-C

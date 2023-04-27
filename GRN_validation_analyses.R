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
library(ggpubr)
library(RColorBrewer)

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

# define colours for the 4 nomicro settings
colours <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")

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
a <- ggplot(val.df %>%                 
         filter(validation == "pcHi-C") %>% # want to measure frequency of recoveries 
         filter(grepl("nomicro", setting)) %>% # just consider GRNs built without microglial cells
         #dplyr::select(gene, peak, setting, resolution, validated, corr.method) %>%  # remove useless columns
         group_by(setting, resolution) %>% # group the data by resolution and setting
         summarise(recoveredTFs = sum(validated)), #%>% # count number of peak-gene links recovered per each resolution and setting
         #as.data.frame() %>%
         #separate(setting, c("dataset", "corr.method", "micro"), sep = "_", remove = FALSE), # create corresponding dataset and corr.method columns to be used in aes()
       aes(resolution, y = recoveredTFs , group = setting, col = setting)) + 
  geom_smooth(size = 1, se = FALSE) +
  labs(y = "# Peak-Gene Links Recovered", x = "Cluster Resolutions", title = "Peak-Gene Connections Recovery") +
  scale_color_manual(values = colours) +  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme_bw() 

b <- ggplot(val.df %>%                 
         filter(validation == "pcHi-C") %>% # want to measure frequency of recoveries 
         filter(grepl("nomicro", setting)) %>% # just consider GRNs built without microglial cells
         group_by(setting, resolution) %>% # group the data by resolution and setting
         summarise(recoveredTFs = sum(validated)), #%>% # count number of peak-gene links recovered per each resolution and setting
       aes(resolution, y = recoveredTFs , group = setting, fill = setting)) + 
  geom_bar(stat = "identity") + 
  facet_grid("setting") +
  labs(y = "# Peak-Gene Links Recovered", x = "Cluster Resolutions", title = "Recovery distribution across resulutions") +
  scale_fill_manual(values = colours) +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.background = element_blank(),
        strip.text.y = element_blank())
  
c <- ggplot(val.df %>%                 
         filter(validation == "pcHi-C") %>% # want to measure frequency of recoveries 
         filter(grepl("nomicro", setting)) %>% # just consider GRNs built without microglial cells
         group_by(setting, resolution) %>% # group the data by resolution and setting
         summarise(recovered = sum(validated), GRN_links = length(validated)) %>% #%>% # count number of peak-gene links recovered per each resolution and setting
         mutate(non_recovered = GRN_links - recovered) %>%
         pivot_longer(c(recovered, non_recovered), names_to = "links", values_to = "value"),
         aes(resolution, y = value , group = setting, fill = links)) + 
  geom_bar(stat = "identity") + 
  facet_grid("setting") +
  labs(y = "# Peak-Gene Links", x = "Cluster Resolutions", title = "Recovered vs. Non-recovered Links") +
  scale_fill_manual(values = c("grey", "#FFCC00")) +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.text.y = element_blank(),
        strip.background = element_rect(fill = colours))

d <- ggplot(val.df %>%                 
         filter(validation == "pcHi-C") %>% # want to measure frequency of recoveries 
         filter(grepl("nomicro", setting)) %>% # just consider GRNs built without microglial cells
         group_by(setting, resolution) %>% # group the data by resolution and setting
         summarise(recovered = sum(validated), GRN_links = length(validated)) %>% #%>% # count number of peak-gene links recovered per each resolution and setting
         mutate(non_recovered = GRN_links - recovered, ratio = recovered / non_recovered),
       aes(resolution, y = ratio , group = setting, col = setting)) + 
  geom_smooth(size = 1, se = FALSE) +
  labs(y = "# Peak-Gene Links Recovery Ratio", x = "Cluster Resolutions", title = "Recovery Ratio") +
  scale_color_manual(values = colours) +  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme_bw() 

pg <- ggarrange(a, b, c, d, ncol = 2, nrow = 2, common.legend = TRUE, labels = "AUTO")
pg
ggsave("/g/scb/zaugg/deuner/valdata/figures/PeakGeneRecovery.pdf", pg)

# genes per regulon
# regions per regulon


################################
# GENERAL GRN STATISTICS PLOTS #
################################

# path of setting
path <- "/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/combined_batch_mode_spearman_nomicro/Batch_Mode_Outputs/"

# define resolutions 
resolutions <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))

# read GRNs
i <- 1
for (res in resolutions){
  print(res)
  grn <- qread(paste0(path, "output_pseudobulk_clusterRes", res, "_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/GRN.qs"))
  grn <- grn@connections$all.filtered$`0`
  grn$resolution <- rep(resolutions[i], nrow(grn))
  if (i == 1) {
    all.grn.links <- grn
  } else {
    all.grn.links <- rbind(all.grn.links, grn)
  }
  i <- i + 1
  
} 

# set colours
coul <- brewer.pal(4, "YlOrRd") 

# Add more colors to this palette :
coul <- colorRampPalette(coul)(19)

# convert resolutions to factor
all.grn.links$resolution <- as.factor(all.grn.links$resolution)

# number of TFs (unique)
p1 <- ggplot(all.grn.links %>%
         dplyr::select(TF.name, TF.ENSEMBL, peak.ID, gene.name, gene.ENSEMBL, resolution) %>% # select important columms
         group_by(resolution) %>%
         summarise(TFcount = n_distinct(TF.ENSEMBL)),
       aes(x = resolution, y = TFcount, group = resolution, fill = resolution)) + 
  geom_bar(stat = "identity", col = "black", fill = coul) + 
  labs(x = "Resolution", y = "TF count", title = "# TFs per resolution") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme_bw()

# number of genes
p2 <- ggplot(all.grn.links %>%
         dplyr::select(TF.name, TF.ENSEMBL, peak.ID, gene.name, gene.ENSEMBL, resolution) %>% # select important columms
         group_by(resolution) %>%
         summarise(GENEcount = n_distinct(gene.ENSEMBL)),
       aes(x = resolution, y = GENEcount, group = resolution, fill = resolution)) + 
  geom_bar(stat = "identity", col = "black", fill = coul) + 
  labs(x = "Resolution", y = "Gene count", title = "# Genes per resolution") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme_bw()

# number of peaks
p3 <- ggplot(all.grn.links %>%
         dplyr::select(TF.name, TF.ENSEMBL, peak.ID, gene.name, gene.ENSEMBL, resolution) %>% # select important columms
         group_by(resolution) %>%
         summarise(PEAKcount = n_distinct(peak.ID)),
       aes(x = resolution, y = PEAKcount, group = resolution, fill = resolution)) + 
  geom_bar(stat = "identity", col = "black", fill = coul) + 
  labs(x = "Resolution", y = "Peak count", title = "# Peaks per resolution") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme_bw()

# number total links
p4 <- ggplot(all.grn.links %>%
         dplyr::select(TF.name, TF.ENSEMBL, peak.ID, gene.name, gene.ENSEMBL, resolution) %>% # select important columms
         group_by(resolution) %>%
         tally(),
       aes(x = resolution, y = n, group = resolution, fill = resolution)) + 
  geom_bar(stat = "identity", col = "black", fill = coul) + 
  labs(x = "Resolution", y = "Link count", title = "# Links per resolution") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme_bw()


# genes per regulon distribution
p5 <- ggplot(all.grn.links %>%
         dplyr::select(TF.name, TF.ENSEMBL, peak.ID, gene.name, gene.ENSEMBL, resolution) %>% # select important columms
         group_by(resolution,TF.ENSEMBL) %>%
         summarise(GENEcount = n_distinct(gene.ENSEMBL)),
         #filter(n < 200), # remove outliers
       aes(x = resolution, y = GENEcount, group = resolution, fill = resolution)) +
  geom_boxplot()

# regions per regulon distribution
p6 <- ggplot(all.grn.links %>%
         dplyr::select(TF.name, TF.ENSEMBL, peak.ID, gene.name, gene.ENSEMBL, resolution) %>% # select important columms
         group_by(resolution,TF.ENSEMBL) %>%
         summarise(PEAKcount = n_distinct(peak.ID)),
       #filter(n < 200), # remove outliers
       aes(x = resolution, y = PEAKcount, group = resolution, fill = resolution)) +
  geom_boxplot()

pgrn <- ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 3, common.legend = TRUE, labels = "AUTO")
pgrn
ggsave("/g/scb/zaugg/deuner/valdata/figures/GRNstats_comb_sp_nm.pdf", pg)


egrn <- qread("/g/scb/zaugg/deuner/GRaNIE/outputdata/combined_spearman_nomciro_Res_0.25/output_pseudobulk_wsnn_res.0.25_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/GRN.qs")
egrn@connections$all.filtered


#############################
# eGRNs VALIDATION ANALYSES #
#############################

## Experimental data used: pcHI-C, ChiP-Seq and eQTL

## input data format 
## rows: links in GRN (ids) , columns: gene/TF, peak, resolution, validated, setting, dataset, corr.method,  
#         gene/TF         peak          resolution      validated            setting           dataset         corr.method     
# link1   FOX1    chr1:187183-187634        0.5           1           combined_spearman     combined          spearman  

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

####################################
# ALL SIGNIFICANT LINKS VALIDATION #
####################################

#%%%%%%%#
# SETUP #
#%%%%%%%#

# get pcHi-C data
pcHiC.df <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(pcHiC.df) <- c("gene", "peak", "res", "validated", "setting", "dataset", "corr.method")

for (file in grep(list.files(pcHiC.dir, pattern = "nomicro_"), pattern = 'filtered', invert=TRUE, value=TRUE)){
  print(file)
  df <- read.csv(paste0(pcHiC.dir, file), row.names = 1)
  dataset <- str_split(file, "_") %>% map(1) %>% as.character
  corr.method <- str_split(file, "_", 2) %>% map(2) %>% as.character %>% strsplit("[.]") %>% map(1) %>% as.character
  df <- df %>%
    tibble::add_column(setting = paste(dataset, corr.method, sep = "_"), 
                      dataset = dataset,
                      corr.method = corr.method) %>%
    dplyr::select(-gene.ENSEMBL)
  pcHiC.df <- rbind(pcHiC.df, df)
}

# get ChiP-seq data
ChiP.df <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(ChiP.df) <- c("TF", "peak", "res", "validated", "setting", "dataset", "corr.method")

for (file in grep(grep(list.files(ChiP.dir, pattern = c("nomicro", "all")), pattern = 'filtered', invert=TRUE, value=TRUE), pattern = 'spearman_nomicro.csv', invert=TRUE, value=TRUE)){
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

for (file in grep(list.files(eQTL.dir, pattern = "nomicro_all"), pattern = 'filtered', invert=TRUE, value=TRUE)){
  print(file)
  df <- read.csv(paste0(eQTL.dir, file), row.names = 1)
  dataset <- str_split(file, "_") %>% map(1) %>% as.character
  corr.method <- str_split(file, "_", 2) %>% map(2) %>% as.character %>% strsplit("[.]") %>% map(1) %>% as.character
  df <- df %>%
    tibble::add_column(setting = paste(dataset, corr.method, sep = "_"), 
                       dataset = dataset,
                       corr.method = corr.method) %>%
    distinct() # the validated links are duplicated many times
  eQTL.df <- rbind(eQTL.df, df)
}

## Merge all 
# change TF column of ChiP-seq to gene so the names match 
names(ChiP.df)[1] <- "gene"
val.df <- do.call("rbind", list(pcHiC.df, ChiP.df, eQTL.df)) # eQTL missing
val.df$validation <- rep(c("pcHi-C", "ChiP-seq", "eQTL"), times = c(nrow(pcHiC.df), nrow(ChiP.df), nrow(eQTL.df)))
val.df$resolution <- as.factor(val.df$resolution)
val.df$validated <- as.numeric(val.df$validated)

# # add gene.ENSEMBL column
# library(biomaRt)
# 
# ensembl <- useEnsembl(biomart = "genes")
# datasets <- listDatasets(ensembl)
# 
# ensembl.con <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 
# attr <- listAttributes(ensembl.con)
# filters <- listFilters(ensembl.con)
# 
# annotLookup <- getBM(
#   attributes = c("ensembl_gene_id", "external_gene_name"),
#   #filter = "external_gene_name", # If I want the ENSEMBL IDs then it is ensembl_gene_id
#   #values = val.df$gene,
#   mart = ensembl.con
# )

#%%%%%%%#
# PLOTS #
#%%%%%%%#

# define colours for the 4 nomicro settings
colours <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")
labels <- c("Combined - Pearson", "Combined - Spearman", "Timecourse - Pearson", "Timecourse - Spearman")

# TFs for which I have Chip-seq data
chip.TFs <- val.df %>%
  dplyr::filter(validation == "ChiP-seq", validated == 1) %>%
  pull(gene) %>%
  unique()

# TF recovery from ChiP-seq data (unique Tfs)
tf1 <- ggplot(val.df %>% 
         dplyr::filter(validation == "ChiP-seq") %>% # want to measure frequency of recoveries, 
         distinct(gene, resolution, setting, .keep_all = TRUE) %>% # just consider the TF once per GRN, it if appears in the eGRN multiple times it has a recovery of 1
         dplyr::select(resolution, setting, validated) %>%  # remove useless columns
         group_by(resolution, setting) %>% # group the data by resolution and setting
         summarise(recoveredTFs = sum(validated)), # count number of TFs recovered per each resolution and setting
       aes(resolution, y = recoveredTFs , group = setting, fill = setting)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = colours, labels = labels) +  
  labs(y = "# TFs recovered", x = "Cluster Resolutions", title = "TF Recovery", caption = "Unique TFs") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid("setting")

# TF recovery from ChiP-seq data (allow duplicate TFs)
tf2<- ggplot(val.df %>% 
         dplyr::filter(validation == "ChiP-seq") %>% # want to measure frequency of recoveries
         #distinct(gene, resolution, setting, .keep_all = TRUE) %>% # just consider the TF once per GRN, it if appears in the eGRN multiple times it has a recovery of 1
         dplyr::select(resolution, setting, validated) %>%  # remove useless columns
         group_by(resolution, setting) %>% # group the data by resolution and setting
         summarise(recoveredTFs = sum(validated)), # count number of TFs recovered per each resolution and setting
       aes(resolution, y = recoveredTFs , group = setting, fill = setting)) + 
  geom_bar(stat = "identity") + 
  labs(y = "# TFs recovered", x = "Cluster Resolutions", title = "TF Recovery", caption = "(allow duplicates)") +
  scale_fill_manual(values = colours, labels = labels) +  
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid("setting")

tfplots <- ggarrange(tf1, tf2, nrow = 1, ncol = 2, common.legend = TRUE) 
tfplots
ggsave("/g/scb/zaugg/deuner/valdata/figures/TFPeakRecovery.pdf", tfplots)

# Check TFs
val.TFs <- val.df %>% 
  dplyr::filter(validation == "ChiP-seq", validated == 1) %>% 
  distinct(gene, resolution, setting, .keep_all = TRUE) %>% 
  group_by(setting, resolution)

for (set in unique(val.df$setting)){
  val.TFs %>%
    dplyr::filter(setting == set) %>%
    dplyr::select(resolution, gene) %>% 
    print(n = 100)
} 

# TF-peak links recovery from ChiP-seq data
tfp1 <- ggplot(val.df %>% 
                dplyr::filter(validation == "ChiP-seq") %>% # want to measure frequency of recoveries, Also restrain background to Chip-seq TFs
                dplyr::select(resolution, setting, validated) %>%  # remove useless columns
                group_by(setting, resolution) %>% # group the data by resolution and setting
                summarise(recoveredTFPEAKs = sum(validated)), # count number of TFs recovered per each resolution and setting
              aes(resolution, y = recoveredTFPEAKs , group = setting, fill = setting)) + 
  geom_bar(stat = "identity") + 
  labs(y = "# TF-peak links recovered", x = "Cluster Resolutions", title = "TF-peak links Recovery") +
  scale_fill_manual(values = colours, labels = labels) +  
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid("setting") 


tfp3 <- ggplot(val.df %>%                 
              dplyr::filter(validation == "ChiP-seq", gene %in% chip.TFs) %>% # want to measure frequency of recoveries, Also restrain background to Chip-seq TFs
              group_by(setting, resolution) %>% # group the data by resolution and setting
              summarise(recovered = sum(validated), GRN_links = length(validated)) %>% #%>% # count number of peak-gene links recovered per each resolution and setting
              mutate(non_recovered = GRN_links - recovered) %>%
              pivot_longer(c(recovered, non_recovered), names_to = "links", values_to = "value"),
            aes(resolution, y = value , group = setting, fill = links)) + 
  geom_bar(stat = "identity") + 
  facet_grid("setting") +
  labs(y = "# TF-PEAK Links", x = "Cluster Resolutions", title = "Recovered vs. Non-recovered Links") +
  scale_fill_manual(values = c("grey", "#FFCC00")) +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.text.y = element_blank(),
        strip.background = element_rect(fill = colours))

tfp4 <- ggplot(val.df %>%                 
              dplyr::filter(validation == "ChiP-seq", gene %in% chip.TFs) %>% # want to measure frequency of recoveries, Also restrain background to Chip-seq TFs
              group_by(setting, resolution) %>% # group the data by resolution and setting
              summarise(recovered = sum(validated), GRN_links = length(validated)) %>% #%>% # count number of peak-gene links recovered per each resolution and setting
              mutate(non_recovered = GRN_links - recovered, ratio = recovered / non_recovered),
            aes(resolution, y = ratio , group = setting, col = setting)) + 
  geom_smooth(size = 1, se = FALSE) +
  labs(y = "# TF-Peak Links Recovery Ratio", x = "Cluster Resolutions", title = "Recovery Ratio") +
  scale_color_manual(values = colours, labels = labels) +  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme_bw() 

tfp <- ggarrange(tfp1, tfp3, tfp4, ncol = 2, nrow = 2, common.legend = TRUE, labels = "AUTO")
tfp
ggsave("/g/scb/zaugg/deuner/valdata/figures/TFPeakRecovery.pdf", tfp)


# Peak-Gene links recovery from pcHi-C
a <- ggplot(val.df %>%                 
         dplyr::filter(validation == "pcHi-C") %>% # want to measure frequency of recoveries 
         #dplyr::select(gene, peak, setting, resolution, validated, corr.method) %>%  # remove useless columns
         group_by(setting, resolution) %>% # group the data by resolution and setting
         summarise(recoveredTFs = sum(validated)), #%>% # count number of peak-gene links recovered per each resolution and setting
         #as.data.frame() %>%
         #separate(setting, c("dataset", "corr.method", "micro"), sep = "_", remove = FALSE), # create corresponding dataset and corr.method columns to be used in aes()
       aes(resolution, y = recoveredTFs , group = setting, col = setting)) + 
  geom_smooth(size = 1, se = FALSE) +
  labs(y = "# Peak-Gene Links Recovered", x = "Cluster Resolutions", title = "Peak-Gene Connections Recovery") +
  scale_color_manual(values = colours, labels = labels) +  
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme_bw() 

a2 <- ggplot(val.df %>%                 
              dplyr::filter(validation == "pcHi-C", validated == 1) %>% # want to measure frequency of recoveries 
              #dplyr::select(gene, peak, setting, resolution, validated, corr.method) %>%  # remove useless columns
              group_by(setting, resolution) %>% # group the data by resolution and setting
              tally(), #%>% # count number of peak-gene links recovered per each resolution and setting
            #as.data.frame() %>%
            #separate(setting, c("dataset", "corr.method", "micro"), sep = "_", remove = FALSE), # create corresponding dataset and corr.method columns to be used in aes()
            aes(resolution, y = n , group = setting, col = setting)) + 
  geom_smooth(size = 1, se = FALSE) +
  labs(y = "# Peak-Gene Links Recovered", x = "Cluster Resolutions", title = "Peak-Gene Connections Recovery") +
  scale_color_manual(values = colours, labels = labels) +  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme_bw() 

b <- ggplot(val.df %>%                 
         dplyr::filter(validation == "pcHi-C") %>% # want to measure frequency of recoveries 
         group_by(setting, resolution) %>% # group the data by resolution and setting
         summarise(recoveredTFs = sum(validated)), #%>% # count number of peak-gene links recovered per each resolution and setting
       aes(resolution, y = recoveredTFs , group = setting, fill = setting)) + 
  geom_bar(stat = "identity") + 
  facet_grid("setting") +
  labs(y = "# Peak-Gene Links Recovered", x = "Cluster Resolutions", title = "Recovery distribution across resulutions") +
  scale_fill_manual(values = colours, labels = labels) +  
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank())
  
c <- ggplot(val.df %>%                 
         dplyr::filter(validation == "pcHi-C") %>% # want to measure frequency of recoveries 
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
         dplyr::filter(validation == "pcHi-C") %>% # want to measure frequency of recoveries 
         group_by(setting, resolution) %>% # group the data by resolution and setting
         summarise(recovered = sum(validated), GRN_links = length(validated)) %>% #%>% # count number of peak-gene links recovered per each resolution and setting
         mutate(non_recovered = GRN_links - recovered, ratio = recovered / non_recovered),
       aes(resolution, y = ratio , group = setting, col = setting)) + 
  geom_smooth(size = 1, se = FALSE) +
  labs(y = "# Peak-Gene Links Recovery Ratio", x = "Cluster Resolutions", title = "Recovery Ratio") +
  scale_color_manual(values = colours, labels = labels) +  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme_bw() 

pg <- ggarrange(a, b, c, d, ncol = 2, nrow = 2, common.legend = TRUE, labels = "AUTO")
pg
ggsave("/g/scb/zaugg/deuner/valdata/figures/PeakGeneRecovery.pdf", pg)

# inspect res 18 of combined pearson nomicro
val.df %>%
  dplyr::filter(setting == "combined_pearson_nomicro") %>%
  group_by(resolution) %>%
  tally()

# check GRN object
grn18 <- qread(paste0("/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/combined_batch_mode_pearson_nomicro/Batch_Mode_Outputs/output_pseudobulk_clusterRes18_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/GRN.qs"))
grn18@connections$peak_genes$`0` %>% dplyr::filter(peak_gene.p_raw < 0.05) %>% nrow
grn18@connections$all.filtered$`0` %>% nrow
grn20 <- qread(paste0("/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/combined_batch_mode_pearson_nomicro/Batch_Mode_Outputs/output_pseudobulk_clusterRes20_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/GRN.qs"))
grn20@connections$peak_genes$`0` %>% dplyr::filter(peak_gene.p_raw < 0.05) %>% nrow
grn20@connections$all.filtered$`0` %>% nrow

# Check GRNs files
GRN18 <- fread(paste0("/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/combined_batch_mode_pearson_nomicro/Batch_Mode_Outputs/output_pseudobulk_clusterRes18_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/connections_TFPeak0.2_peakGene0.1.tsv.gz"))
nrow(GRN18)
GRN16 <- fread(paste0("/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/combined_batch_mode_pearson_nomicro/Batch_Mode_Outputs/output_pseudobulk_clusterRes16_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/connections_TFPeak0.2_peakGene0.1.tsv.gz"))
nrow(GRN16)
GRN20 <- fread(paste0("/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/combined_batch_mode_pearson_nomicro/Batch_Mode_Outputs/output_pseudobulk_clusterRes20_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/connections_TFPeak0.2_peakGene0.1.tsv.gz"))
nrow(GRN20)
# This is correct

# Check pcHi-C validated links
pchic <- fread("/g/scb/zaugg/deuner/valdata/pcHi-C/validated_links/combined_pearson_nomicro.csv")
pchic <- dplyr::select(pchic, -V1)

pchic %>%
  group_by(resolution) %>%
  tally()

pchic2 <- fread("/g/scb/zaugg/deuner/valdata/pcHi-C/validated_links/combined_pearson_nomicro_all.csv")
pchic2 <- dplyr::select(pchic2, -V1)

pchic2 %>%
  group_by(resolution) %>%
  tally()
# Problem is here - check pchic validation file

# Peak - Gene validation from eQTL data

colours <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")

i <- ggplot(val.df %>%                 
         dplyr::filter(validation == "eQTL") %>% # want to measure frequency of recoveries 
         #filter(grepl("nomicro", setting)) %>% # just consider GRNs built without microglial cells
         #dplyr::select(gene, peak, setting, resolution, validated, corr.method) %>%  # remove useless columns
         group_by(setting, resolution) %>% # group the data by resolution and setting
         summarise(recoveredLINKs = sum(validated)), #%>% # count number of peak-gene links recovered per each resolution and setting
       #as.data.frame() %>%
       #separate(setting, c("dataset", "corr.method", "micro"), sep = "_", remove = FALSE), # create corresponding dataset and corr.method columns to be used in aes()
       aes(resolution, y = recoveredLINKs , group = setting, col = setting)) + 
  geom_smooth(size = 1, se = FALSE) +
  labs(y = "# Peak-Gene Links Recovered", x = "Cluster Resolutions", title = "Peak-Gene Connections Recovery") +
  scale_color_manual(values = colours, labels = labels) +  
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme_bw()
i
ii <- ggplot(val.df %>%                 
              dplyr::filter(validation == "eQTL") %>% # want to measure frequency of recoveries 
              #filter(grepl("nomicro", setting)) %>% # just consider GRNs built without microglial cells
              group_by(setting, resolution) %>% # group the data by resolution and setting
              summarise(recoveredLINKs = sum(validated)), #%>% # count number of peak-gene links recovered per each resolution and setting
            aes(resolution, y = recoveredLINKs , group = setting, fill = setting)) + 
  geom_bar(stat = "identity") + 
  facet_grid("setting") +
  labs(y = "# Peak-Gene Links Recovered", x = "Cluster Resolutions", title = "Recovery distribution across resulutions") +
  scale_fill_manual(values = colours, labels = labels) +  
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank())

iii <- ggplot(val.df %>%                 
              dplyr::filter(validation == "eQTL") %>% # want to measure frequency of recoveries 
              #filter(grepl("nomicro", setting)) %>% # just consider GRNs built without microglial cells
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

iv <- ggplot(val.df %>%                 
              dplyr::filter(validation == "eQTL") %>% # want to measure frequency of recoveries 
              #filter(grepl("nomicro", setting)) %>% # just consider GRNs built without microglial cells
              group_by(setting, resolution) %>% # group the data by resolution and setting
              summarise(recovered = sum(validated), GRN_links = length(validated)) %>% #%>% # count number of peak-gene links recovered per each resolution and setting
              mutate(non_recovered = GRN_links - recovered, ratio = recovered / non_recovered),
            aes(resolution, y = ratio , group = setting, col = setting)) + 
  geom_smooth(size = 1, se = FALSE) +
  labs(y = "# Peak-Gene Links Recovery Ratio", x = "Cluster Resolutions", title = "Recovery Ratio") +
  scale_color_manual(values = colours, labels = labels) +  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme_bw() 

pge <- ggarrange(i, ii, iii, iv, ncol = 2, nrow = 2, common.legend = TRUE, labels = "AUTO")
pge
ggsave("/g/scb/zaugg/deuner/valdata/figures/PeakGeneRecovery_eQTL.png", pge, device = "png")

# validated links with pcHi-C and eQTL data show a similar pattern
# check if links with pcHi-C data are the same as links validated with eQTL data

# paper figure plots #
a <- a + theme(legend.title=element_blank(),
               axis.title.x=element_blank(),
               axis.title.y=element_blank(),
               axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + labs(title = "")
i <- i + theme(legend.title=element_blank(), 
               axis.title.x=element_blank(),
               axis.title.y=element_blank(),
               axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + labs(title = "")
tfp1 <- tfp1 + theme(legend.title=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank()) + labs(title = "")

fp <- ggarrange(a, i, ggplot() + theme_transparent(), tfp1, common.legend = T, ncol = 2, nrow = 2)

tiff(paste0("/g/scb/zaugg/deuner/figs/validation_plots", ".pdf"), units="in", width=6.4, height=4.8, res=300, type = "cairo")
fp
dev.off()


# Define the colors for the lines
colours <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")

# peak-gene recovery distribution plot
# Combine the data for both validation types
combined_data <- rbind(
  val.df %>%                 
    dplyr::filter(validation == "pcHi-C") %>%
    group_by(setting, resolution) %>%
    summarise(recoveredLINKs = sum(validated)),
  val.df %>%                 
    dplyr::filter(validation == "eQTL") %>%
    group_by(setting, resolution) %>%
    summarise(recoveredLINKs = sum(validated))
)
combined_data$validation <- c(rep("pcHi-C", nrow(combined_data)/2), rep("eQTL", nrow(combined_data)/2))

new_settings <- c()
for (i in 1:nrow(combined_data)){
  split <- str_split(combined_data$setting[i], "_")
if (split[[1]][1] == "combined"){
  new_setting <- ifelse(split[[1]][2] == "spearman", "Combined - Spearman", "Combined - Pearson") 
} else {
  new_setting <- ifelse(split[[1]][2] == "spearman", "Timecourse - Spearman", "Timecourse - Pearson") 
}
  new_settings <- c(new_settings, new_setting)
}
combined_data$setting <- new_settings


# Create the plot
peakgene_dist <- ggplot(combined_data, aes(resolution, y = recoveredLINKs, group = interaction(setting, validation), col = setting, linetype = validation)) +
  geom_smooth(size = 1, se = FALSE) +
  labs(y = "No. Recovered Links", x = "Cluster Resolutions", title = "Peak-Gene Recovery") +
  scale_color_manual(values = colours) +
  scale_linetype_manual(values = c("pcHi-C" = "solid", "eQTL" = "dashed")) +
  theme_bw() + 
  theme(legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5))

# Display the plot
peakgene_dist


# peak-gene ratio dist
# Combine the data for both validation types
combined_data_ratio <- rbind(
  val.df %>%                 
    dplyr::filter(validation == "pcHi-C") %>%
    group_by(setting, resolution) %>%
    summarise(recovered = sum(validated), GRN_links = length(validated)) %>% 
    mutate(non_recovered = GRN_links - recovered, ratio = recovered / non_recovered),
  val.df %>%                 
    dplyr::filter(validation == "eQTL") %>%
    group_by(setting, resolution) %>%
    summarise(recovered = sum(validated), GRN_links = length(validated)) %>% 
    mutate(non_recovered = GRN_links - recovered, ratio = recovered / non_recovered)
)
combined_data_ratio$validation <- c(rep("pcHi-C", nrow(combined_data)/2), rep("eQTL", nrow(combined_data)/2))

new_settings <- c()
for (i in 1:nrow(combined_data_ratio)){
  split <- str_split(combined_data_ratio$setting[i], "_")
  if (split[[1]][1] == "combined"){
    new_setting <- ifelse(split[[1]][2] == "spearman", "Combined - Spearman", "Combined - Pearson") 
  } else {
    new_setting <- ifelse(split[[1]][2] == "spearman", "Timecourse - Spearman", "Timecourse - Pearson") 
  }
  new_settings <- c(new_settings, new_setting)
}
combined_data_ratio$setting <- new_settings


# Create the plot
peakgene_ratio <- ggplot(combined_data_ratio, aes(resolution, y = ratio, group = interaction(setting, validation), col = setting, linetype = validation)) +
  geom_smooth(size = 1, se = FALSE) +
  labs(y = "Ratio", x = "Cluster Resolutions", title = "") +
  scale_color_manual(values = colours) +
  scale_linetype_manual(values = c("pcHi-C" = "solid", "eQTL" = "dashed")) +
  theme_bw() + 
  theme(legend.title = element_blank())

# Display the plot
peakgene_ratio


# tf recovery plot #
tfrecovery <- ggplot(val.df %>% 
         dplyr::filter(validation == "ChiP-seq") %>% # want to measure frequency of recoveries, Also restrain background to Chip-seq TFs
         dplyr::select(resolution, setting, validated) %>%  # remove useless columns
         group_by(setting, resolution) %>% # group the data by resolution and setting
         summarise(recoveredTFPEAKs = sum(validated)), # count number of TFs recovered per each resolution and setting
       aes(resolution, y = recoveredTFPEAKs , group = setting, fill = setting)) + 
  geom_bar(stat = "identity", col = "black") + 
  labs(y = "", x = "Cluster Resolutions", title = "TF-peak Recovery") +
  scale_fill_manual(values = colours, labels = labels) +  
  theme_bw() + 
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid("setting") 
  

all_plots <- ggarrange(peakgene_dist, tfrecovery, peakgene_ratio, common.legend = T, labels = c("A", "B", "C")) 

tiff(paste0("/g/scb/zaugg/deuner/figs/validation_plots_unfiltered", ".tiff"), units="in", width=6.4, height=4.8, res=300, type = "cairo")
all_plots
dev.off()

#############################
# FILTERED LINKS VALIDATION #
#############################

#%%%%%%%#
# SETUP #
#%%%%%%%#

# get pcHi-C data
pcHiC.df <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(pcHiC.df) <- c("gene", "peak", "res", "validated", "setting", "dataset", "corr.method")

for (file in list.files(pcHiC.dir, pattern = "nomicro_filtered")){
  print(file)
  df <- read.csv(paste0(pcHiC.dir, file), row.names = 1)
  dataset <- str_split(file, "_") %>% map(1) %>% as.character
  corr.method <- str_split(file, "_", 2) %>% map(2) %>% as.character %>% strsplit("[.]") %>% map(1) %>% as.character
  df <- df %>%
    tibble::add_column(setting = paste(dataset, corr.method, sep = "_"), 
                       dataset = dataset,
                       corr.method = corr.method) %>%
    dplyr::select(-gene.ENSEMBL)
  pcHiC.df <- rbind(pcHiC.df, df)
}

# get ChiP-seq data
ChiP.df <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(ChiP.df) <- c("TF", "peak", "res", "validated", "setting", "dataset", "corr.method")

for (file in list.files(ChiP.dir, pattern = "nomicro_filtered")){
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

for (file in list.files(eQTL.dir, pattern = "nomicro_filtered")){
  print(file)
  df <- read.csv(paste0(eQTL.dir, file), row.names = 1)
  dataset <- str_split(file, "_") %>% map(1) %>% as.character
  corr.method <- str_split(file, "_", 2) %>% map(2) %>% as.character %>% strsplit("[.]") %>% map(1) %>% as.character
  df <- df %>%
    tibble::add_column(setting = paste(dataset, corr.method, sep = "_"), 
                       dataset = dataset,
                       corr.method = corr.method) %>%
    distinct() # the validated links are duplicated many times
  eQTL.df <- rbind(eQTL.df, df)
}

## Merge all 
# change TF column of ChiP-seq to gene so the names match 
names(ChiP.df)[1] <- "gene"
val.df <- do.call("rbind", list(pcHiC.df, ChiP.df, eQTL.df)) # eQTL missing
val.df$validation <- rep(c("pcHi-C", "ChiP-seq", "eQTL"), times = c(nrow(pcHiC.df), nrow(ChiP.df), nrow(eQTL.df)))
val.df$resolution <- as.factor(val.df$resolution)
val.df$validated <- as.numeric(val.df$validated)

#%%%%%%%#
# PLOTS #
#%%%%%%%#

# define colours for the 4 nomicro settings
colours <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")

# TF recovery from ChiP-seq data
# TF recovery from ChiP-seq data (unique Tfs)
tf1 <- ggplot(val.df %>% 
                dplyr::filter(validation == "ChiP-seq") %>% # want to measure frequency of recoveries
                distinct(gene, resolution, setting, .keep_all = TRUE) %>% # just consider the TF once per GRN, it if appears in the eGRN multiple times it has a recovery of 1
                dplyr::select(resolution, setting, validated) %>%  # remove useless columns
                group_by(resolution, setting) %>% # group the data by resolution and setting
                summarise(recoveredTFs = sum(validated)), # count number of TFs recovered per each resolution and setting
              aes(resolution, y = recoveredTFs , group = setting, fill = setting)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = colours) +  
  labs(y = "# TFs recovered", x = "Cluster Resolutions", title = "TF Recovery", caption = "Unique TFs") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid("setting")

# TF recovery from ChiP-seq data (allow duplicate TFs)
tf2<- ggplot(val.df %>% 
               dplyr::filter(validation == "ChiP-seq") %>% # want to measure frequency of recoveries
               #distinct(gene, resolution, setting, .keep_all = TRUE) %>% # just consider the TF once per GRN, it if appears in the eGRN multiple times it has a recovery of 1
               dplyr::select(resolution, setting, validated) %>%  # remove useless columns
               group_by(resolution, setting) %>% # group the data by resolution and setting
               summarise(recoveredTFs = sum(validated)), # count number of TFs recovered per each resolution and setting
             aes(resolution, y = recoveredTFs , group = setting, fill = setting)) + 
  geom_bar(stat = "identity") + 
  labs(y = "# TFs recovered", x = "Cluster Resolutions", title = "TF Recovery", caption = "(allow duplicates)") +
  scale_fill_manual(values = colours) +  
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid("setting")

tfplots <- ggarrange(tf1, tf2, nrow = 1, ncol = 2, common.legend = TRUE)
tfplots
ggsave("/g/scb/zaugg/deuner/valdata/figures/TFPeakRecovery_filtered.pdf", tfplots)

# TF-peak links recovery from ChiP-seq data
tfp1 <- ggplot(val.df %>% 
                 dplyr::filter(validation == "ChiP-seq") %>% # want to measure frequency of recoveries
                 dplyr::select(resolution, setting, validated) %>%  # remove useless columns
                 group_by(setting, resolution) %>% # group the data by resolution and setting
                 summarise(recoveredTFPEAKs = sum(validated)), # count number of TFs recovered per each resolution and setting
               aes(resolution, y = recoveredTFPEAKs , group = setting, fill = setting)) + 
  geom_bar(stat = "identity") + 
  labs(y = "# TF-peak links recovered", x = "Cluster Resolutions", title = "TF-peak links Recovery") +
  scale_fill_manual(values = colours, labels = labels) +  
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid("setting") 


tfp3 <- ggplot(val.df %>%                 
                 dplyr::filter(validation == "ChiP-seq", gene %in% chip.TFs) %>% # want to measure frequency of recoveries, Also restrain background to Chip-seq TFs 
                 group_by(setting, resolution) %>% # group the data by resolution and setting
                 summarise(recovered = sum(validated), GRN_links = length(validated)) %>% #%>% # count number of peak-gene links recovered per each resolution and setting
                 mutate(non_recovered = GRN_links - recovered) %>%
                 pivot_longer(c(recovered, non_recovered), names_to = "links", values_to = "value"),
               aes(resolution, y = value , group = setting, fill = links)) + 
  geom_bar(stat = "identity") + 
  facet_grid("setting") +
  labs(y = "# TF-PEAK Links", x = "Cluster Resolutions", title = "Recovered vs. Non-recovered Links") +
  scale_fill_manual(values = c("grey", "#FFCC00")) +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.text.y = element_blank(),
        strip.background = element_rect(fill = colours))

tfp4 <- ggplot(val.df %>%                 
                 dplyr::filter(validation == "ChiP-seq", gene %in% chip.TFs) %>% # want to measure frequency of recoveries, Also restrain background to Chip-seq TFs 
                 group_by(setting, resolution) %>% # group the data by resolution and setting
                 summarise(recovered = sum(validated), GRN_links = length(validated)) %>% #%>% # count number of peak-gene links recovered per each resolution and setting
                 mutate(non_recovered = GRN_links - recovered, ratio = recovered / non_recovered),
               aes(resolution, y = ratio , group = setting, col = setting)) + 
  geom_smooth(size = 1, se = FALSE) +
  labs(y = "# TF-Peak Links Recovery Ratio", x = "Cluster Resolutions", title = "Recovery Ratio") +
  scale_color_manual(values = colours) +  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme_bw() 

tfp <- ggarrange(tfp1, tfp3, tfp4, ncol = 2, nrow = 2, common.legend = TRUE, labels = "AUTO")
tfp
ggsave("/g/scb/zaugg/deuner/valdata/figures/TFPeakRecovery_filtered.pdf", pg)

# Check TFs
val.df %>%                 
  dplyr::filter(validation == "ChiP-seq",
                validated == 1,
                setting == "combined_spearman_nomicro_filtered") %>% 
  distinct() %>% 
  pull(gene) %>%
  unique()

# Peak-Gene links recovery from pcHi-C
a <- ggplot(val.df %>%                 
              dplyr::filter(validation == "pcHi-C") %>% # want to measure frequency of recoveries 
              #dplyr::select(gene, peak, setting, resolution, validated, corr.method) %>%  # remove useless columns
              group_by(setting, resolution) %>% # group the data by resolution and setting
              summarise(recoveredTFs = sum(validated)), #%>% # count number of peak-gene links recovered per each resolution and setting
            #as.data.frame() %>%
            #separate(setting, c("dataset", "corr.method", "micro"), sep = "_", remove = FALSE), # create corresponding dataset and corr.method columns to be used in aes()
            aes(resolution, y = recoveredTFs , group = setting, col = setting)) + 
  geom_smooth(size = 1, se = FALSE) +
  labs(y = "# Peak-Gene Links Recovered", x = "Cluster Resolutions", title = "Peak-Gene Connections Recovery") +
  scale_color_manual(values = colours, labels = labels) +  
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank())

b <- ggplot(val.df %>%                 
              dplyr::filter(validation == "pcHi-C") %>% # want to measure frequency of recoveries 
              group_by(setting, resolution) %>% # group the data by resolution and setting
              summarise(recoveredTFs = sum(validated)), #%>% # count number of peak-gene links recovered per each resolution and setting
            aes(resolution, y = recoveredTFs , group = setting, fill = setting)) + 
  geom_bar(stat = "identity") + 
  facet_grid("setting") +
  labs(y = "# Peak-Gene Links Recovered", x = "Cluster Resolutions", title = "Recovery distribution across resolutions") +
  scale_fill_manual(values = colours) +  
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank())

c <- ggplot(val.df %>%                 
              dplyr::filter(validation == "pcHi-C") %>% # want to measure frequency of recoveries 
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
              dplyr::filter(validation == "pcHi-C") %>% # want to measure frequency of recoveries 
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
ggsave("/g/scb/zaugg/deuner/valdata/figures/PeakGeneRecovery_filtered.pdf", pg)


# Peak - Gene validation from eQTL data

colours <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "purple")

i <- ggplot(val.df %>%                 
              dplyr::filter(validation == "eQTL") %>% # want to measure frequency of recoveries 
              #filter(grepl("nomicro", setting)) %>% # just consider GRNs built without microglial cells
              #dplyr::select(gene, peak, setting, resolution, validated, corr.method) %>%  # remove useless columns
              group_by(setting, resolution) %>% # group the data by resolution and setting
              summarise(recoveredLINKs = sum(validated)), #%>% # count number of peak-gene links recovered per each resolution and setting
            #as.data.frame() %>%
            #separate(setting, c("dataset", "corr.method", "micro"), sep = "_", remove = FALSE), # create corresponding dataset and corr.method columns to be used in aes()
            aes(resolution, y = recoveredLINKs , group = setting, col = setting)) + 
  geom_smooth(size = 1, se = FALSE) +
  labs(y = "# Peak-Gene Links Recovered", x = "Cluster Resolutions", title = "Peak-Gene Connections Recovery") +
  scale_color_manual(values = colours, labels = labels) +  
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme_bw()
i
ii <- ggplot(val.df %>%                 
               dplyr::filter(validation == "eQTL") %>% # want to measure frequency of recoveries 
               #filter(grepl("nomicro", setting)) %>% # just consider GRNs built without microglial cells
               group_by(setting, resolution) %>% # group the data by resolution and setting
               summarise(recoveredLINKs = sum(validated)), #%>% # count number of peak-gene links recovered per each resolution and setting
             aes(resolution, y = recoveredLINKs , group = setting, fill = setting)) + 
  geom_bar(stat = "identity") + 
  facet_grid("setting") +
  labs(y = "# Peak-Gene Links Recovered", x = "Cluster Resolutions", title = "Recovery distribution across resulutions") +
  scale_fill_manual(values = colours) +  
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank())

iii <- ggplot(val.df %>%                 
                dplyr::filter(validation == "eQTL") %>% # want to measure frequency of recoveries 
                #filter(grepl("nomicro", setting)) %>% # just consider GRNs built without microglial cells
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

iv <- ggplot(val.df %>%                 
               dplyr::filter(validation == "eQTL") %>% # want to measure frequency of recoveries 
               #filter(grepl("nomicro", setting)) %>% # just consider GRNs built without microglial cells
               group_by(setting, resolution) %>% # group the data by resolution and setting
               summarise(recovered = sum(validated), GRN_links = length(validated)) %>% #%>% # count number of peak-gene links recovered per each resolution and setting
               mutate(non_recovered = GRN_links - recovered, ratio = recovered / non_recovered),
             aes(resolution, y = ratio , group = setting, col = setting)) + 
  geom_smooth(size = 1, se = FALSE) +
  labs(y = "# Peak-Gene Links Recovery Ratio", x = "Cluster Resolutions", title = "Recovery Ratio") +
  scale_color_manual(values = colours) +  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme_bw() 

pge <- ggarrange(i, ii, iii, iv, ncol = 2, nrow = 2, common.legend = TRUE, labels = "AUTO")
pge
ggsave("/g/scb/zaugg/deuner/valdata/figures/PeakGeneRecovery_eQTL_filtered.png", pge, device = "png")

# paper figure plots #

# Define the colors for the lines
colours <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")

# peak-gene recovery distribution plot
# Combine the data for both validation types
combined_data <- rbind(
  val.df %>%                 
    dplyr::filter(validation == "pcHi-C") %>%
    group_by(setting, resolution) %>%
    summarise(recoveredLINKs = sum(validated)),
  val.df %>%                 
    dplyr::filter(validation == "eQTL") %>%
    group_by(setting, resolution) %>%
    summarise(recoveredLINKs = sum(validated))
)
combined_data$validation <- c(rep("pcHi-C", nrow(combined_data)/2), rep("eQTL", nrow(combined_data)/2))

new_settings <- c()
for (i in 1:nrow(combined_data)){
  split <- str_split(combined_data$setting[i], "_")
  if (split[[1]][1] == "combined"){
    new_setting <- ifelse(split[[1]][2] == "spearman", "Combined - Spearman", "Combined - Pearson") 
  } else {
    new_setting <- ifelse(split[[1]][2] == "spearman", "Timecourse - Spearman", "Timecourse - Pearson") 
  }
  new_settings <- c(new_settings, new_setting)
}
combined_data$setting <- new_settings


# Create the plot
peakgene_dist <- ggplot(combined_data, aes(resolution, y = recoveredLINKs, group = interaction(setting, validation), col = setting, linetype = validation)) +
  geom_smooth(size = 1, se = FALSE) +
  labs(y = "No. Recovered Links", x = "Cluster Resolutions", title = "Peak-Gene Recovery") +
  scale_color_manual(values = colours) +
  scale_linetype_manual(values = c("pcHi-C" = "solid", "eQTL" = "dashed")) +
  theme_bw() + 
  theme(legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5))

# Display the plot
peakgene_dist


# peak-gene ratio dist
# Combine the data for both validation types
combined_data_ratio <- rbind(
  val.df %>%                 
    dplyr::filter(validation == "pcHi-C") %>%
    group_by(setting, resolution) %>%
    summarise(recovered = sum(validated), GRN_links = length(validated)) %>% 
    mutate(non_recovered = GRN_links - recovered, ratio = recovered / non_recovered),
  val.df %>%                 
    dplyr::filter(validation == "eQTL") %>%
    group_by(setting, resolution) %>%
    summarise(recovered = sum(validated), GRN_links = length(validated)) %>% 
    mutate(non_recovered = GRN_links - recovered, ratio = recovered / non_recovered)
)
combined_data_ratio$validation <- c(rep("pcHi-C", nrow(combined_data)/2), rep("eQTL", nrow(combined_data)/2))

new_settings <- c()
for (i in 1:nrow(combined_data_ratio)){
  split <- str_split(combined_data_ratio$setting[i], "_")
  if (split[[1]][1] == "combined"){
    new_setting <- ifelse(split[[1]][2] == "spearman", "Combined - Spearman", "Combined - Pearson") 
  } else {
    new_setting <- ifelse(split[[1]][2] == "spearman", "Timecourse - Spearman", "Timecourse - Pearson") 
  }
  new_settings <- c(new_settings, new_setting)
}
combined_data_ratio$setting <- new_settings


# Create the plot
peakgene_ratio <- ggplot(combined_data_ratio, aes(resolution, y = ratio, group = interaction(setting, validation), col = setting, linetype = validation)) +
  geom_smooth(size = 1, se = FALSE) +
  labs(y = "Ratio", x = "Cluster Resolutions", title = "") +
  scale_color_manual(values = colours) +
  scale_linetype_manual(values = c("pcHi-C" = "solid", "eQTL" = "dashed")) +
  theme_bw() + 
  theme(legend.title = element_blank())

# Display the plot
peakgene_ratio


# tf recovery plot #
tfrecovery <- ggplot(val.df %>% 
                       dplyr::filter(validation == "ChiP-seq") %>% # want to measure frequency of recoveries, Also restrain background to Chip-seq TFs
                       dplyr::select(resolution, setting, validated) %>%  # remove useless columns
                       group_by(setting, resolution) %>% # group the data by resolution and setting
                       summarise(recoveredTFPEAKs = sum(validated)), # count number of TFs recovered per each resolution and setting
                     aes(resolution, y = recoveredTFPEAKs , group = setting, fill = setting)) + 
  geom_bar(stat = "identity", col = "black") + 
  labs(y = "", x = "Cluster Resolutions", title = "TF-peak Recovery") +
  scale_fill_manual(values = colours, labels = labels) +  
  theme_bw() + 
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid("setting") 


all_plots <- ggarrange(peakgene_dist, tfrecovery, peakgene_ratio, common.legend = T, labels = c("E", "F", "G")) 

tiff(paste0("/g/scb/zaugg/deuner/figs/validation_plots_filtered", ".tiff"), units="in", width=6.4, height=4.8, res=300, type = "cairo")
all_plots
dev.off()

################################
# GENERAL GRN STATISTICS PLOTS #
################################

# path of setting
path <- "/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/combined_batch_mode_spearman_nomicro/Batch_Mode_Outputs/"

# define resolutions 
resolutions <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))

# df were all links are going to be stores
all.grn.links <- as.data.frame(matrix(nrow = 0, ncol = 16))
names(all.grn.links) <- c("TF.ID",                  "TF.name",                "TF.ENSEMBL",             "TF_peak.r_bin",          "TF_peak.r",             
                          "TF_peak.fdr",            "TF_peak.fdr_direction",  "TF_peak.connectionType", "peak.ID",                "peak_gene.distance",    
                          "peak_gene.r",            "peak_gene.p_raw",        "peak_gene.p_adj",        "gene.ENSEMBL",           "gene.name",             
                          "gene.type")
# read GRNs
for (res in resolutions){
  print(res)
  grn <- qread(paste0(path, "output_pseudobulk_clusterRes", res, "_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/GRN.qs"))
  grn <- grn@connections$all.filtered$`0`
  names(grn)
    
  grn$resolution <- rep(res, nrow(grn))
  
  all.grn.links <- rbind(all.grn.links, grn)
  # if (res == 0.1) {
  #   all.grn.links <- grn
  # } else {
  #   all.grn.links <- rbind(all.grn.links, grn)
  # }
  
} 

# convert resolutions to factor
all.grn.links$resolution <- as.factor(all.grn.links$resolution)

# number of TFs (unique)
p1 <- ggplot(all.grn.links %>%
         dplyr::select(TF.name, TF.ENSEMBL, peak.ID, gene.name, gene.ENSEMBL, resolution) %>% # select important columms
         group_by(resolution) %>%
         summarise(TFcount = n_distinct(TF.ENSEMBL)),
       aes(x = resolution, y = TFcount, group = resolution, fill = TFcount)) + 
  geom_bar(stat = "identity", col = "black") + 
  labs(x = "Resolution", y = "TF count", title = "# TFs per resolution") + 
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5), legend.position="none") + 
  scale_fill_gradient(low = "#132B43",
                      high = "#56B1F7")

# number of genes
p2 <- ggplot(all.grn.links %>%
         dplyr::select(TF.name, TF.ENSEMBL, peak.ID, gene.name, gene.ENSEMBL, resolution) %>% # select important columms
         group_by(resolution) %>%
         summarise(GENEcount = n_distinct(gene.ENSEMBL)),
       aes(x = resolution, y = GENEcount, group = resolution, fill = GENEcount)) + 
  geom_bar(stat = "identity", col = "black") + 
  labs(x = "Resolution", y = "Gene count", title = "# Genes per resolution") + 
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5), legend.position="none") + 
  scale_fill_gradient(low = "#006633",
                      high = "#99FFCC")

# number of peaks
p3 <- ggplot(all.grn.links %>%
         dplyr::select(TF.name, TF.ENSEMBL, peak.ID, gene.name, gene.ENSEMBL, resolution) %>% # select important columms
         group_by(resolution) %>%
         summarise(PEAKcount = n_distinct(peak.ID)),
       aes(x = resolution, y = PEAKcount, group = resolution, fill = PEAKcount)) + 
  geom_bar(stat = "identity", col = "black") + 
  labs(x = "Resolution", y = "Peak count", title = "# Peaks per resolution") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust=1), plot.title = element_text(hjust = 0.5), legend.position="none") + 
  scale_fill_gradient(low = "#FF9900",
                      high = "#FFFF99")

# number total links
p4 <- ggplot(all.grn.links %>%
         dplyr::select(TF.name, TF.ENSEMBL, peak.ID, gene.name, gene.ENSEMBL, resolution) %>% # select important columms
         group_by(resolution) %>%
         tally(),
       aes(x = resolution, y = n, group = resolution, fill = n)) + 
  geom_bar(stat = "identity", col = "black") + 
  labs(x = "Resolution", y = "Link count", title = "# Links per resolution") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust=1), plot.title = element_text(hjust = 0.5), legend.position="none") + 
  scale_fill_gradient(low = "#9966CC",
                      high = "#FFCCFF")
  

pgrn <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
pgrn
ggsave("/g/scb/zaugg/deuner/valdata/figures/GRNstats_comb_sp_nm.png", pgrn, device = "png")

tiff(paste0("/g/scb/zaugg/deuner/figs/GRNstats_comb_sp_nm", ".tiff"), units="in", width=6.4, height=4.8, res=300, type = "cairo")
pgrn
dev.off()

# genes per regulon distribution
p5 <- ggplot(all.grn.links %>%
         dplyr::select(TF.name, TF.ENSEMBL, peak.ID, gene.name, gene.ENSEMBL, resolution) %>% # select important columms
         group_by(resolution,TF.ENSEMBL) %>%
         summarise(GENEcount = n_distinct(gene.ENSEMBL)),
         #filter(n < 200), # remove outliers
       aes(x = resolution, y = GENEcount, group = resolution, fill = resolution)) +
  geom_boxplot() +
  labs(x = "Resolution", y = "# Genes per regulon", title = "Genes per Regulon Distribution") + 
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.position="none")

p6 <- ggplot(all.grn.links %>%
               dplyr::select(TF.name, TF.ENSEMBL, peak.ID, gene.name, gene.ENSEMBL, resolution) %>% # select important columms
               group_by(resolution,TF.ENSEMBL) %>%
               summarise(GENEcount = n_distinct(gene.ENSEMBL)) %>%
               dplyr::filter(GENEcount < 20), # remove outliers
             aes(x = resolution, y = GENEcount, group = resolution, fill = resolution)) +
  geom_boxplot() + 
  labs(x = "Resolution", y = "# Genes per regulon", title = "Genes per Regulon Distribution") + 
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.position="none")

# regions per regulon distribution
p7 <- ggplot(all.grn.links %>%
         dplyr::select(TF.name, TF.ENSEMBL, peak.ID, gene.name, gene.ENSEMBL, resolution) %>% # select important columms
         group_by(resolution,TF.ENSEMBL) %>%
         summarise(PEAKcount = n_distinct(peak.ID)),
       #filter(n < 200), # remove outliers
       aes(x = resolution, y = PEAKcount, group = resolution, fill = resolution)) +
  geom_boxplot() + 
  labs(x = "Resolution", y = "# Peaks per regulon", title = "Peaks per Regulon Distribution") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust=1), plot.title = element_text(hjust = 0.5),
        legend.position="none")


p8 <- ggplot(all.grn.links %>%
               dplyr::select(TF.name, TF.ENSEMBL, peak.ID, gene.name, gene.ENSEMBL, resolution) %>% # select important columms
               group_by(resolution,TF.ENSEMBL) %>%
               summarise(PEAKcount = n_distinct(peak.ID)) %>%
               dplyr::filter(PEAKcount < 20), # remove outliers
             aes(x = resolution, y = PEAKcount, group = resolution, fill = resolution)) +
  geom_boxplot() + 
  labs(x = "Resolution", y = "# Peaks per regulon", title = "Peaks per Regulon Distribution") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust=1), plot.title = element_text(hjust = 0.5),
        legend.position="none")

prgn <- ggarrange(p5, p6, p7, p8, ncol = 2, nrow = 2, labels = "AUTO")
prgn
ggsave("/g/scb/zaugg/deuner/valdata/figures/RegulonsStats_comb_sp_nm.png", prgn, device = "png")


# genes per regulon distribution (log)
p9 <- ggplot(all.grn.links %>%
               dplyr::select(TF.name, TF.ENSEMBL, peak.ID, gene.name, gene.ENSEMBL, resolution) %>% # select important columms
               group_by(resolution,TF.ENSEMBL) %>%
               summarise(GENEcount = n_distinct(gene.ENSEMBL)),
             #filter(n < 200), # remove outliers
             aes(x = resolution, y = log(GENEcount), group = resolution, fill = resolution)) +
  geom_boxplot() +
  labs(x = "Resolution", y = "# Genes per regulon (log)", title = "Genes per Regulon Distribution") + 
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.position="none")

# peaks per regulon distribution (log)
p10 <- ggplot(all.grn.links %>%
               dplyr::select(TF.name, TF.ENSEMBL, peak.ID, gene.name, gene.ENSEMBL, resolution) %>% # select important columms
               group_by(resolution,TF.ENSEMBL) %>%
               summarise(PEAKcount = n_distinct(peak.ID)),
             #filter(n < 200), # remove outliers
             aes(x = resolution, y = log(PEAKcount), group = resolution, fill = resolution)) +
  geom_boxplot() + 
  labs(x = "Resolution", y = "# Peaks per regulon (log)", title = "Peaks per Regulon Distribution") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust=1), plot.title = element_text(hjust = 0.5),
        legend.position="none")

prgn.log <- ggarrange(p9, p10, nrow = 2)
prgn.log
ggsave("/g/scb/zaugg/deuner/valdata/figures/RegulonsStats_log_comb_sp_nm.png", prgn.log, device = "png")

tiff(paste0("/g/scb/zaugg/deuner/figs/RegulonsStats_log_comb_sp_nm", ".tiff"), units="in", width=6.4, height=4.8, res=300, type = "cairo")
prgn.log
dev.off()

# arrange paper plots
tiff(paste0("/g/scb/zaugg/deuner/figs/eGRN_stats_regulons_cb_sp_nolabels", ".tiff"), units="in", width=12.8, height=4.8, res=300, type = "cairo")
ggarrange(pgrn, prgn.log, ncol = 2)
dev.off()

ggarrange(pgrn, ggplot() + theme_void(), prgn.log, ncol = 3, widths = c(1, 0.2, 1))

cb <- ggarrange(pgrn, prgn.log, ncol = 2)

################################################
# Score Resolutions based on previous analyses #
################################################

# create function to evaluate which is the best resolution based on validation scores
evaluateRes <- function(data,                  # df with resolutions as rows and validation stats as columns. (data.frame)
                        res,                   # resolution to evaluate. (int)
                        qc,                    # parameter that indicates if that resolution has passed the peak-gene QC filters. 0 = No. 1 = Yes. (int)
                        pcHiC.weight = 0.5,    # weight given to the pcHi-C validation. (int). Default = 0.5. 
                        ChiP.weight = 0.1,     # weight given to the ChiP-seq validation. (int). Default = 0.1.
                        eQTL.weight = 0.4      # weight given to the eQTL validation. (int). Default = 0.4 
) {
  
  # data format:
  #   res  pcHiC.val   ChiP.val    eQTL.val    peak.gene.links    tf.peak.links
  #   0.1     130         3           85          1320                  754
  #   ...     ...        ...         ...          ...                   ...

  # or 
  
  #         gene/TF         peak          resolution      validated            setting           dataset         corr.method     
  # link1   FOX1    chr1:187183-187634        0.5           1           combined_spearman     combined          spearman  
  # At the moment it the function is implemented for this format
  
  
  # Scoring function: (0-1)
  #     f(eGRN) = x(0.5y + 0.1z + 0.4t)
  
  
  # get number of validated TF-peak links from ChiP-seq
  ChiP.val <- val.df %>% 
    dplyr::filter(validation == "ChiP-seq") %>% # want to measure frequency of recoveries
    dplyr::select(resolution, setting, validated) %>%  # remove useless columns
    group_by(setting, resolution) %>% # group the data by resolution and setting
    summarise(recoveredTFPEAKs = sum(validated)) %>% # count number of TFs recovered per each resolution and setting
    dplyr::filter(resolution == res, setting == "combined_spearman_nomicro") %>% # filter for resolution and setting
    pull(recoveredTFPEAKs) # get number of TF-peak links

  # get number of validated peak-gene links from pcHi-C
  pcHiC.val <- val.df %>% 
    dplyr::filter(validation == "pcHi-C") %>% # want to measure frequency of recoveries
    dplyr::select(resolution, setting, validated) %>%  # remove useless columns
    group_by(setting, resolution) %>% # group the data by resolution and setting
    summarise(recoveredLINKs = sum(validated)) %>% # count number of TFs recovered per each resolution and setting
    dplyr::filter(resolution == res, setting == "combined_spearman_nomicro") %>% # filter for resolution and setting
    pull(recoveredLINKs) # get number of TF-peak links
  
  # get number of validated peak-gene links from eQTL
  # eQTL.val <- val.df %>% 
  #   dplyr::filter(validation == "eQTL") %>% # want to measure frequency of recoveries
  #   dplyr::select(resolution, setting, validated) %>%  # remove useless columns
  #   group_by(setting, resolution) %>% # group the data by resolution and setting
  #   summarise(recoveredLINKs = sum(validated)) %>% # count number of TFs recovered per each resolution and setting
  #   dplyr::filter(resolution == res, setting == "combined_spearman_nomicro") %>% # filter for resolution and setting
  #   pull(recoveredLINKs) # get number of TF-peak links
  eQTL.val <- 0
  
  # get total number of inferred TF-peak links
  tf.peak.links <- val.df %>%                 
    dplyr::filter(validation == "ChiP-seq") %>% # want to measure frequency of recoveries 
    group_by(setting, resolution) %>% # group the data by resolution and setting
    summarise(recovered = sum(validated), GRN_links = length(validated)) %>% # count number of peak-gene links recovered per each resolution and setting
    dplyr::filter(resolution == res, setting == "combined_spearman_nomicro") %>% # filter for resolution and setting
    pull(GRN_links) # get number of total TF-peak links
  
  # get total number of inferred peak-gene links
  peak.gene.links <- val.df %>%                 
    dplyr::filter(validation == "pcHi-C") %>% # want to measure frequency of recoveries 
    group_by(setting, resolution) %>% # group the data by resolution and setting
    summarise(recovered = sum(validated), GRN_links = length(validated)) %>% # count number of peak-gene links recovered per each resolution and setting
    dplyr::filter(resolution == res, setting == "combined_spearman_nomicro") %>% # filter for resolution and setting
    pull(GRN_links) # get number of total peak-gene links
  
  # peak-gene QC score
  x <- qc
  
  # pcHi-C validation score
  y <- pcHiC.val**2 / peak.gene.links
  
  # ChiP-seq validation score
  z <- ChiP.val**2 / tf.peak.links
  
  # eQTL validation score
  t <- eQTL.val / peak.gene.links
  
  # compute function
  out <- x * ( (pcHiC.weight * y) + (ChiP.weight * z) + (eQTL.weight * t) )
  
  print(paste("Res", res, "Evaluation Score:", out, sep = " "))
}


# set qc scores vector (combined spearman)
qc.vec <- c(rep(0, 4), rep(1, 15))

i <- 1
for (r in resolutions){
  evaluateRes(data = val.df, res = r, qc = qc.vec[i])
  i <- i + 1
  }


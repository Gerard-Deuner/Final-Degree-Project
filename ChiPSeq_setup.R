##################
# ChiP-Seq SETUP #
##################

# load libraries
library(dplyr)
library(tidyr)
library(Seurat)
library(qs)
library(GenomicRanges)
library(stringr)

# read dataset from command line
args <- commandArgs(trailingOnly = TRUE)

# define dataset
dataset <- args[1] # timecourse | combined

# set path to ChiP-seq data
path <- "/g/scb2/zaugg/deuner/ChipSeqData/"

# cell types for which I have Chip-Seq data
celltypes <- c("cortical-interneuron", "dopaminergic-neuron", "hiPSC", "neural-progenitor", "neuron", "neuron-progenitor", "NPC")

# dataframe where all the filtered chip-seq data will be stored
chip <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(chip) <- c("V1", "V2")

# get gene names (hence TF names) present in the dataset
s.obj <- qread(paste0("/g/scb/zaugg/deuner/GRaNIE/tmp/", dataset, ".pp.seuratObject.qs"))
DefaultAssay(s.obj) <- "RNA"
tc.gene.names <- rownames(s.obj)

# iterate over bed files and filter for TFs present in the dataset
for (celltype in celltypes){
  bed <- read.csv(paste0(path, "remap2022_", celltype, "_all_macs2_hg38_v1_0.bed.gz"), sep = ".", header = FALSE)
  bed <- bed %>% dplyr::filter(V2 %in% tc.gene.names) %>% dplyr::select(c("V1", "V2")) 
  bed <- bed %>% mutate(V3 = rep(celltype, nrow(bed)))
  print(celltype)
  chip <- rbind(chip, bed)
}

# split first column and add appropriate column names 
chip <- chip %>% separate(V1, c("chr", "chromStart", "chromEnd", "exp"), "\t")
colnames(chip)[5:6] <- c("TF", "celltype")
head(chip)

# filter for peaks found in the dataset
DefaultAssay(s.obj) <- "ATAC"
tc.peaks <- rownames(s.obj) %>% as.data.frame()
names(tc.peaks) <- "peak"
tc.peaks[c("chr", "chromStart", "chromEnd")] <- str_split_fixed(tc.peaks$peak, "-", 3)

# use findOverlaps method from GenomicRanges
ref <-  makeGRangesFromDataFrame(tc.peaks, keep.extra.columns = TRUE) 
qry <-  makeGRangesFromDataFrame(chip)
ovlp <- findOverlaps(qry, ref) # returns indexes of intersecting regions
head(ovlp)
ovlp.indx <- ovlp@from

# subset baits overlapping dataset gene promoters
chip <-  chip[ovlp.indx,]

# add peak column
chip$peak <- tc.peaks[ovlp@to, "peak"]
head(chip)

# important columns: TF and peak. Keep also celltype column
chip <- chip %>%
  select(TF, peak, celltype)

# save file
write.csv(chip, file = gzfile(paste0("/g/scb/zaugg/deuner/valdata/ChiP-seq/filtered_chipseq_", dataset, ".tsv.gz")))


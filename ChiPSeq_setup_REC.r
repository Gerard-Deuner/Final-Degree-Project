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

# define dataset
dataset <- "timeocourse"

# set path to ChiP-seq data
path <- "/g/scb2/zaugg/deuner/ChipSeqData/"

# cell types for which I have Chip-Seq data
celltypes <- c("cortical-interneuron", "dopaminergic-neuron", "hiPSC", "neural-progenitor", "neuron", "neuron-progenitor", "NPC")

# dataframe where all the filtered chip-seq data will be stored
chip <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(chip) <- c("V1", "V2")

# get gene names (hence TF names) present in timecourse dataset
timecourse.s <- qread("/g/scb/zaugg/deuner/GRaNIE/tmp/timecourse.pp.seuratObject.qs")
DefaultAssay(timecourse.s) <- "RNA"
tc.gene.names <- rownames(timecourse.s)

# iterate over bed files and filter for TFs present in the timecourse dataset
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

# save file
write.csv(chip, file = gzfile("/g/scb2/zaugg/deuner/ChipSeqData/filtered_TFs_chipseq_timecourse.csv.gz"))


# filter for peaks found in the timecourse dataset
DefaultAssay(timecourse.s) <- "ATAC"
tc.peaks <- rownames(timecourse.s)
tc.peaks <- as.data.frame(str_split_fixed(tc.peaks, "-", 3))
colnames(tc.peaks) <- c("chr", "chromStart", "chromEnd")
# use findOverlaps method from GenomicRanges
ref <-  makeGRangesFromDataFrame(tc.peaks, keep.extra.columns = TRUE) 
qry <-  makeGRangesFromDataFrame(chip)
ovlp <- findOverlaps(qry, ref) # returns indexes of intersecting regions
head(ovlp)
ovlp.indx <- unique(ovlp@to)
# subset baits overlapping timecourse gene promoters
chip <-  chip[ovlp.indx,]
chip <- chip[,2:ncol(chip)]
head(chip)

# save file
write.csv(chip, file = gzfile("/g/scb2/zaugg/deuner/ChipSeqData/filtered_chipseq_timecourse.csv.gz"))

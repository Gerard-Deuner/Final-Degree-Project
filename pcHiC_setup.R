################
# PCHi-C SETUP #
################

# load libraries
library(dplyr)
library(Seurat)
library(qs)
library(GenomicRanges)
library(tidyr)
library(stringr)
library(purrr)

# read dataset from command line
args <- commandArgs(trailingOnly = TRUE)

# define dataset
dataset <- args[1] # timecourse | combined

#%%%%%%%%#
# NEURON # 
#%%%%%%%%#

celltype <- "neuron"

# pull the promoter coordinates for the genes in the dataset (either timecourse or combined)
## take the TSS and extend 2kb upstream if strand +, downstream if strand -
# load the respective dataset seurat object to get the gene names present in the dataset
s.obj <- qread(paste0("/g/scb/zaugg/deuner/GRaNIE/tmp/", dataset, ".pp.seuratObject.qs"))
DefaultAssay(s.obj) <- "RNA"
ds.gene.names <- rownames(s.obj)

# get gene strand and position details from gtf file
gtf.file = "/g/scb/zaugg/marttine/RNA_ATAC_integration/annotation/gencode.v41.chr_patch_hapl_scaff.annotation.gtf.gz"
gtf <- rtracklayer::import(gtf.file)
gtf = as.data.frame(gtf)

# subset genes from tgf file that are present in the dataset in bed format
ds.genes.ann <- gtf %>% dplyr::filter((gene_name %in% ds.gene.names) & (type == "gene")) %>% dplyr::select(seqnames, start, end, strand, gene_name, gene_id)
colnames(ds.genes.ann)[1] <- "chr"

# get the promoter regions from each gene
prom.cords <- ds.genes.ann
prom.cords$start <- ifelse(prom.cords$strand == "+", ds.genes.ann$start-2000, ds.genes.ann$end)
prom.cords$end <- ifelse(prom.cords$strand == "+", ds.genes.ann$start, ds.genes.ann$end+2000)
prom.cords
# keep baits that intersect with promoter coordinates 
baits.file <- "/g/scb/zaugg/marttine/RNA_ATAC_integration/CRISPR/GRN/baits.csv"
baits <- read.csv(baits.file, row.names = 1)
# use findOverlaps method from GenomicRanges
qry <-  makeGRangesFromDataFrame(prom.cords, keep.extra.columns = TRUE) 
ref <-  makeGRangesFromDataFrame(baits)
ovlp <- findOverlaps(qry, ref) # returns indexes of intersecting regions
ovlp
ovlp.indx <- ovlp@to
# subset baits overlapping timecourse gene promoters
ds.baits <-  baits[ovlp.indx,]
head(ds.baits)
# keep the gene names to which the baits map to
ds.baits <- ds.baits %>% cbind(gene = prom.cords$gene_name[ovlp@from])
ds.baits
# look at all regions linked to these baits from links file, they essentially represent enhancers
links.file <- "/g/scb/zaugg/marttine/RNA_ATAC_integration/CRISPR/GRN/pchic.i3Neuron.Song2019.hg38.csv"
links <- read.csv(links.file, sep = " ")
# get links from baits found in timecourse data
ds.baits_ehcs <- ds.baits %>% right_join(links, by = c("chr", "start", "end"), multiple = "all")

# filter for these enhancers that are present in my peakset
## (expand peak +/- 1kb and check for any overlap with pchic enhancers)
ehcs <- data.frame(chr = ds.baits_ehcs$chr.1, start = ds.baits_ehcs$start.1, end = ds.baits_ehcs$end.1) 
# get peaks present in timecourse data
DefaultAssay(s.obj) <- "ATAC"
peaks.vect <- rownames(s.obj)
peaks <- data.frame(do.call(rbind, strsplit(peaks.vect, split = "-")))
colnames(peaks) <- c("chr", "start", "end")
# expand peaks +-1Kb
peaks.ext <- peaks
peaks.ext$start <- as.numeric(peaks$start) - 1000
peaks.ext$end <- as.numeric(peaks$end) + 1000
peaks.ext$peak <- peaks.vect

# find overlap between the dataset and HiC enhancers
qry.peaks <-  makeGRangesFromDataFrame(ehcs)
ref.peaks <-  makeGRangesFromDataFrame(peaks.ext, keep.extra.columns = TRUE)
ovlp.peaks <- findOverlaps(qry.peaks, ref.peaks) # returns indexes of intersecting regions
ovlp.peaks.indx <- ovlp.peaks@from
peak.names <- peaks.vect[ovlp.peaks@to]
ds.baits_ds.ehcs <- ds.baits_ehcs[ovlp.peaks.indx,]
ds.baits_ds.ehcs <- ds.baits_ds.ehcs %>% cbind(peak = peak.names)
#View(ds.baits_ds.ehcs)

# keep gene and peak columns
ds.baits_ds.ehcs <- ds.baits_ds.ehcs %>%
  select(gene, peak)

# save links into a csv file 
write.csv(ds.baits_ds.ehcs, paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/links/", dataset, "_", celltype, "_pchic_links.tsv"))


#%%%%%#
# NPC # 
#%%%%%#

# This dataset does not contain a baits file and has two links files, a promoter-promoter interactions file and a promoter-other 
#  interactions files

celltype <- "NPC"

# load links files
# pp.links <- read.csv("/g/scb/zaugg/deuner/valdata/pcHi-C/NPC/GSE86189_npc.pp.all.txt.bz2", sep= "\t") # promoter | promoter links
links <- read.csv("/g/scb/zaugg/deuner/valdata/pcHi-C/NPC/GSE86189_npc.po.all.txt.bz2", sep = "\t") # promoter | other links -> interested in this ones as other may represent enhancers
head(links)

# change format of frag1 and frag2 columns to chr:start-end
links$frag1 <- sub("[.]", ":", links$frag1) # replaces the first match
links$frag1 <- sub("[.]", "-", links$frag1)
links$frag2 <- sub("[.]", ":", links$frag2) # replaces the first match
links$frag2 <- sub("[.]", "-", links$frag2)
head(links)

# pcHiC read data was mapped to human genome hg19 so we need to adapt the hg19 annotations to hg38 reference genome
# save coordinates in separated files to convert them using UCSC liftOver webpage
write.table(links$frag1, paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/", celltype, "/prom_cords_hg19.csv"), row.names=FALSE, col.names = FALSE, quote = FALSE)
write.table(links$frag2, paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/", celltype, "/link_cords_hg19.csv"), row.names=FALSE, col.names = FALSE, quote = FALSE)

# adapt coordinates from hg19 to hg38
# using UCSC liftover: https://genome.ucsc.edu/cgi-bin/hgLiftOver
# Load converted regions as well as non intersected regions to eliminate them
proms.hg38 <- scan(paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/", celltype, "/prom_cords_hg38.bed"), character(), quote = "\t")
proms.hg38.err <- scan(paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/", celltype, "/prom_cords_hg38_error.txt"), character(), quote = "")
links.hg38 <- scan(paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/", celltype, "/link_cords_hg38.bed"), character(), quote = "")
links.hg38.err <- scan(paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/", celltype, "/link_cords_hg38_error.txt"), character(), quote = "")

# filter out links containing an unmatched region and add new regions
links1 <- links %>%
  filter(!frag1 %in% proms.hg38.err) %>%
  mutate(new1 = proms.hg38)

# check if unmatched regions have been removed effectively
nrow(links1) == length(proms.hg38)

# filter out links containing an unmatched region and add new regions
links2 <- links %>%
  filter(!frag2 %in% links.hg38.err) %>%
  mutate(new2 = links.hg38)

# check if unamatched regions have been removed effectively
nrow(links2) == length(links.hg38)

# keep regions involving coords well matched 
links3 <- inner_join(links1, links2 %>% select(frag1, frag2, new2), by = c("frag1", "frag2"))
head(links3)

# check if there is a decrease of links
nrow(links)
nrow(links3)
nrow(links3) <= nrow(links1) & nrow(links2)

links <- links3
rm(links3)

# replace old cords (hg19) with new cords (hg38) 
links$frag1 <- links$new1
links$frag2 <- links$new2
links <- links %>% select(-c("new1", "new2"))
head(links)


# create individual columns for chr, start, end respectively 
# 1. promoter
# 2. "enhancer"
links[c("chr.1", "start.1")] <- str_split_fixed(links$frag1, ":", 2) 
links[c("start.1", "end.1")] <- str_split_fixed(links$start.1, "-", 2) 
links[c("chr.2", "start.2")] <- str_split_fixed(links$frag2, ":", 2)
links[c("start.2", "end.2")] <- str_split_fixed(links$start.2, "-", 2) 
head(links)

# create ranges column 
links$ranges1 <- strsplit(links$frag1, ":", 2) %>% map(2)
links$ranges2 <- strsplit(links$frag2, ":", 2) %>% map(2)
head(links)

# create GRanges objects 
links1 <- makeGRangesFromDataFrame(links, seqnames.field = "chr.1", start.field = "start.1", end.field = "end.1", keep.extra.columns = TRUE) %>%
  as.data.frame() %>%
  rename(seqnames = "chr.1", start = "start.1", end = "end.1") %>%
  dplyr::select(c("chr.1", "start.1", "end.1", "frag1", "frag2"))

links2 <- makeGRangesFromDataFrame(links, seqnames.field = "chr.2", start.field = "start.2", end.field = "end.2", keep.extra.columns = TRUE) %>%
  as.data.frame() %>% 
  rename(seqnames = "chr.2", start = "start.2", end = "end.2") %>%
  dplyr::select(c("chr.2", "start.2", "end.2", "frag1", "frag2"))

nrow(links1)
nrow(links2)

# pull the promoter coordinates for the genes in the dataset (either timecourse or combined)
## take the TSS and extend 2kb upstream if strand +, downstream if strand -
# load the respective dataset seurat object to get the gene names present in the dataset
s.obj <- qread(paste0("/g/scb/zaugg/deuner/GRaNIE/tmp/", dataset, ".pp.seuratObject.qs"))
DefaultAssay(s.obj) <- "RNA"
ds.gene.names <- rownames(s.obj)

# get gene strand and position details from gtf file
gtf.file = "/g/scb/zaugg/marttine/RNA_ATAC_integration/annotation/gencode.v41.chr_patch_hapl_scaff.annotation.gtf.gz"
gtf <- rtracklayer::import(gtf.file)
gtf = as.data.frame(gtf)

# subset genes from tgf file that are present in the dataset in bed format
ds.genes.ann <- gtf %>% 
  dplyr::filter((gene_name %in% ds.gene.names) & (type == "gene")) %>% 
  dplyr::select(seqnames, start, end, strand, gene_name, gene_id)
colnames(ds.genes.ann)[1] <- "chr"

# get the promoter regions from each gene
prom.cords <- ds.genes.ann
prom.cords$start <- ifelse(prom.cords$strand == "+", ds.genes.ann$start-2000, ds.genes.ann$end)
prom.cords$end <- ifelse(prom.cords$strand == "+", ds.genes.ann$start, ds.genes.ann$end+2000)
head(prom.cords)

# format links in bed file format
# treat frag1 as baits (promoter regions)
baits <- links
baits[c("chr", "start")] <- str_split_fixed(baits$frag1, ":", 2)
baits[c("start", "end")] <- str_split_fixed(baits$start, "-", 2) 
baits <- baits %>%
  dplyr::select(c("chr", "start", "end", "frag2"))
head(baits)  

# keep promoter regions from links that intersect with promoter coordinates
# use findOverlaps method from GenomicRanges
qry <-  makeGRangesFromDataFrame(prom.cords, keep.extra.columns = TRUE) 
ref <-  makeGRangesFromDataFrame(baits)
ovlp <- findOverlaps(qry, ref) # returns indexes of intersecting regions
ovlp
ovlp.indx <- ovlp@to

# subset baits overlapping timecourse gene promoters
ds.baits <-  baits[ovlp.indx,]
head(ds.baits)
# keep the gene names to which the baits map to
ds.baits <- ds.baits %>% cbind(gene = prom.cords$gene_name[ovlp@from])
head(ds.baits)

# format baits in dataset so they match the links format
ds.baits$frag1 <- do.call(paste, c(ds.baits[1:2], sep= ":"))
ds.baits$frag1 <- do.call(paste, c(ds.baits[c("frag1", "end")], sep= "-"))
ds.baits <- ds.baits %>% select(-c("chr", "start", "end"))
head(ds.baits)

# look at all regions linked to these baits from links file, they essentially represent enhancers
# get links from baits found in timecourse data
ds.baits_ehcs <- ds.baits #frag2 represent the enhancers
ds.baits_ehcs[c("chr", "start")] <- str_split_fixed(ds.baits_ehcs$frag2, ":", 2)
ds.baits_ehcs[c("start", "end")] <- str_split_fixed(ds.baits_ehcs$start, "-", 2) 
head(ds.baits_ehcs)

# filter for these enhancers that are present in my peakset
## (expand peak +/- 1kb and check for any overlap with pchic enhancers)
ehcs <- data.frame(chr = ds.baits_ehcs$chr, start = ds.baits_ehcs$start, end = ds.baits_ehcs$end) 
# get peaks present in timecourse data
DefaultAssay(s.obj) <- "ATAC"
peaks.vect <- rownames(s.obj)
peaks <- data.frame(do.call(rbind, strsplit(peaks.vect, split = "-")))
colnames(peaks) <- c("chr", "start", "end")
head(peaks)
# expand peaks +-1Kb
peaks.ext <- peaks
peaks.ext$start <- as.numeric(peaks$start) - 1000
peaks.ext$end <- as.numeric(peaks$end) + 1000
peaks.ext$peak <- peaks.vect
head(peaks.ext)

# find overlap between the dataset and HiC enhancers
qry.peaks <-  makeGRangesFromDataFrame(ehcs, seqnames.field = "chr", start.field = "start", end.field = "end")
ref.peaks <-  makeGRangesFromDataFrame(peaks.ext, seqnames.field = "chr", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
ovlp.peaks <- findOverlaps(qry.peaks, ref.peaks) # returns indexes of intersecting regions
ovlp.peaks.indx <- ovlp.peaks@from
peak.names <- peaks.vect[ovlp.peaks@to]
ds.baits_ds.ehcs <- ds.baits_ehcs[ovlp.peaks.indx,]
ds.baits_ds.ehcs <- ds.baits_ds.ehcs %>% cbind(peak = peak.names)
head(ds.baits_ds.ehcs)
#View(ds.baits_ds.ehcs)

# clean it
ds.baits_ds.ehcs <- ds.baits_ds.ehcs %>%
  select(gene, peak)

# save links into a csv file 
# IMPORTANT: 2 mandatory columns: peak and gene
write.csv(ds.baits_ds.ehcs, paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/links/", dataset, "_", celltype, "_pchic_links.tsv"))

#%%%%%%#
# iPSC # 
#%%%%%%#

celltype <- "iPSC"

# (same procedure as for NPC)
# ALSO HG19

# pcHI-C files
# capt-CHiCAGO_interactions-iPSC-1.ibed
# capt-CHiCAGO_interactions-iPSC-2.ibed
# capt-CHiCAGO_interactions-iPSC-3.ibed

# load links files
links.file1 <- read.csv("/g/scb/zaugg/deuner/valdata/pcHi-C/iPSC/CM+iPSC_Hi-C/capt-CHiCAGO_interactions-iPSC-1.ibed", sep = "\t")
links.file2 <- read.csv("/g/scb/zaugg/deuner/valdata/pcHi-C/iPSC/CM+iPSC_Hi-C/capt-CHiCAGO_interactions-iPSC-2.ibed", sep = "\t")
links.file3 <- read.csv("/g/scb/zaugg/deuner/valdata/pcHi-C/iPSC/CM+iPSC_Hi-C/capt-CHiCAGO_interactions-iPSC-3.ibed", sep = "\t")

# merge them
links <- rbind(links.file1, links.file2, links.file3)
nrow(links)
head(links)

# create frag1 (promoter) and frag2(enhancer) columns so the subsequent procedure is more straightforward. format chr:start-end
links$frag1 <- do.call(paste, c(links[c("bait_chr", "bait_start")], sep= ":"))
links$frag1 <- do.call(paste, c(links[c("frag1", "bait_end")], sep= "-"))
links$frag2 <- do.call(paste, c(links[c("otherEnd_chr", "otherEnd_start")], sep= ":"))
links$frag2 <- do.call(paste, c(links[c("frag2", "otherEnd_end")], sep= "-"))
head(links)

# pcHiC read data was mapped to human genome hg19 so we need to adapt the hg19 annotations to hg38 reference genome
# save coordinates in separated files to convert them using UCSC liftOver webpage
write.table(links$frag1, paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/", celltype, "/prom_cords_hg19.csv"), row.names=FALSE, col.names = FALSE, quote = FALSE)
write.table(links$frag2, paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/", celltype, "/link_cords_hg19.csv"), row.names=FALSE, col.names = FALSE, quote = FALSE)

# adapt coordinates from hg19 to hg38
# using UCSC liftover: https://genome.ucsc.edu/cgi-bin/hgLiftOver
# Load converted regions as well as non intersected regions to eliminate them
proms.hg38 <- scan(paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/", celltype, "/prom_cords_hg38.bed"), character(), quote = "\t")
proms.hg38.err <- scan(paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/", celltype, "/prom_cords_hg38_error.txt"), character(), quote = "")
links.hg38 <- scan(paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/", celltype, "/link_cords_hg38.bed"), character(), quote = "")
links.hg38.err <- scan(paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/", celltype, "/link_cords_hg38_error.txt"), character(), quote = "")

# filter out links containing an unmatched region and add new regions
links1 <- links %>%
  filter(!frag1 %in% proms.hg38.err) %>%
  mutate(new1 = proms.hg38)

# check if unmatched regions have been removed effectively
nrow(links1) == length(proms.hg38)

# filter out links containing an unmatched region and add new regions
links2 <- links %>%
  filter(!frag2 %in% links.hg38.err) %>%
  mutate(new2 = links.hg38)

# check if unamatched regions have been removed effectively
nrow(links2) == length(links.hg38)

# keep regions involving coords well matched 
links3 <- inner_join(links1, links2 %>% select(frag1, frag2, new2), by = c("frag1", "frag2"), multiple = "any")
head(links3)

# check if there is a decrease of links
nrow(links)
nrow(links3)
nrow(links3) <= nrow(links1) & nrow(links2)

links <- links3
rm(links3)

# replace old cords (hg19) with new cords (hg38) 
links$frag1 <- links$new1
links$frag2 <- links$new2
links <- links %>% select(-c("new1", "new2"))
head(links)


# create individual columns for chr, start, end respectively 
# 1. promoter
# 2. "enhancer"
links[c("chr.1", "start.1")] <- str_split_fixed(links$frag1, ":", 2) 
links[c("start.1", "end.1")] <- str_split_fixed(links$start.1, "-", 2) 
links[c("chr.2", "start.2")] <- str_split_fixed(links$frag2, ":", 2)
links[c("start.2", "end.2")] <- str_split_fixed(links$start.2, "-", 2) 
head(links)

# create ranges column 
links$ranges1 <- strsplit(links$frag1, ":", 2) %>% map(2)
links$ranges2 <- strsplit(links$frag2, ":", 2) %>% map(2)
head(links)

# create GRanges objects 
links1 <- makeGRangesFromDataFrame(links, seqnames.field = "chr.1", start.field = "start.1", end.field = "end.1", keep.extra.columns = TRUE) %>%
  as.data.frame() %>%
  rename(seqnames = "chr.1", start = "start.1", end = "end.1") %>%
  dplyr::select(c("chr.1", "start.1", "end.1", "frag1", "frag2"))

links2 <- makeGRangesFromDataFrame(links, seqnames.field = "chr.2", start.field = "start.2", end.field = "end.2", keep.extra.columns = TRUE) %>%
  as.data.frame() %>% 
  rename(seqnames = "chr.2", start = "start.2", end = "end.2") %>%
  dplyr::select(c("chr.2", "start.2", "end.2", "frag1", "frag2"))

nrow(links1)
nrow(links2)

# pull the promoter coordinates for the genes in the dataset (either timecourse or combined)
## take the TSS and extend 2kb upstream if strand +, downstream if strand -
# load the respective dataset seurat object to get the gene names present in the dataset
s.obj <- qread(paste0("/g/scb/zaugg/deuner/GRaNIE/tmp/", dataset, ".pp.seuratObject.qs"))
DefaultAssay(s.obj) <- "RNA"
ds.gene.names <- rownames(s.obj)

# get gene strand and position details from gtf file
gtf.file = "/g/scb/zaugg/marttine/RNA_ATAC_integration/annotation/gencode.v41.chr_patch_hapl_scaff.annotation.gtf.gz"
gtf <- rtracklayer::import(gtf.file)
gtf = as.data.frame(gtf)

# subset genes from tgf file that are present in the dataset in bed format
ds.genes.ann <- gtf %>% 
  dplyr::filter((gene_name %in% ds.gene.names) & (type == "gene")) %>% 
  dplyr::select(seqnames, start, end, strand, gene_name, gene_id)
colnames(ds.genes.ann)[1] <- "chr"

# get the promoter regions from each gene
prom.cords <- ds.genes.ann
prom.cords$start <- ifelse(prom.cords$strand == "+", ds.genes.ann$start-2000, ds.genes.ann$end)
prom.cords$end <- ifelse(prom.cords$strand == "+", ds.genes.ann$start, ds.genes.ann$end+2000)
head(prom.cords)

# format links in bed file format
# treat frag1 as baits (promoter regions)
baits <- links
baits[c("chr", "start")] <- str_split_fixed(baits$frag1, ":", 2)
baits[c("start", "end")] <- str_split_fixed(baits$start, "-", 2) 
baits <- baits %>%
  dplyr::select(c("chr", "start", "end", "frag2"))
head(baits)  

# keep promoter regions from links that intersect with promoter coordinates
# use findOverlaps method from GenomicRanges
qry <-  makeGRangesFromDataFrame(prom.cords, keep.extra.columns = TRUE) 
ref <-  makeGRangesFromDataFrame(baits)
ovlp <- findOverlaps(qry, ref) # returns indexes of intersecting regions
ovlp
ovlp.indx <- ovlp@to

# subset baits overlapping timecourse gene promoters
ds.baits <-  baits[ovlp.indx,]
head(ds.baits)
# keep the gene names to which the baits map to
ds.baits <- ds.baits %>% cbind(gene = prom.cords$gene_name[ovlp@from])
head(ds.baits)

# format baits in dataset so they match the links format
ds.baits$frag1 <- do.call(paste, c(ds.baits[1:2], sep= ":"))
ds.baits$frag1 <- do.call(paste, c(ds.baits[c("frag1", "end")], sep= "-"))
ds.baits <- ds.baits %>% select(-c("chr", "start", "end"))
head(ds.baits)

# look at all regions linked to these baits from links file, they essentially represent enhancers
# get links from baits found in timecourse data
ds.baits_ehcs <- ds.baits #frag2 represent the enhancers
ds.baits_ehcs[c("chr", "start")] <- str_split_fixed(ds.baits_ehcs$frag2, ":", 2)
ds.baits_ehcs[c("start", "end")] <- str_split_fixed(ds.baits_ehcs$start, "-", 2) 
head(ds.baits_ehcs)

# filter for these enhancers that are present in my peakset
## (expand peak +/- 1kb and check for any overlap with pchic enhancers)
ehcs <- data.frame(chr = ds.baits_ehcs$chr, start = ds.baits_ehcs$start, end = ds.baits_ehcs$end) 
# get peaks present in timecourse data
DefaultAssay(s.obj) <- "ATAC"
peaks.vect <- rownames(s.obj)
peaks <- data.frame(do.call(rbind, strsplit(peaks.vect, split = "-")))
colnames(peaks) <- c("chr", "start", "end")
head(peaks)
# expand peaks +-1Kb
peaks.ext <- peaks
peaks.ext$start <- as.numeric(peaks$start) - 1000
peaks.ext$end <- as.numeric(peaks$end) + 1000
peaks.ext$peak <- peaks.vect
head(peaks.ext)

# find overlap between the dataset and HiC enhancers
qry.peaks <-  makeGRangesFromDataFrame(ehcs, seqnames.field = "chr", start.field = "start", end.field = "end")
ref.peaks <-  makeGRangesFromDataFrame(peaks.ext, seqnames.field = "chr", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
ovlp.peaks <- findOverlaps(qry.peaks, ref.peaks) # returns indexes of intersecting regions
ovlp.peaks.indx <- ovlp.peaks@from
peak.names <- peaks.vect[ovlp.peaks@to]
ds.baits_ds.ehcs <- ds.baits_ehcs[ovlp.peaks.indx,]
ds.baits_ds.ehcs <- ds.baits_ds.ehcs %>% cbind(peak = peak.names)
head(ds.baits_ds.ehcs)
#View(ds.baits_ds.ehcs)

# clean it
ds.baits_ds.ehcs <- ds.baits_ds.ehcs %>%
  select(gene, peak)

# save links into a csv file 
# IMPORTANT: 2 mandatory columns: peak and gene
write.csv(ds.baits_ds.ehcs, paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/links/", dataset, "_", celltype, "_pchic_links.tsv"))

#%%%%%#
# ALL #
#%%%%%#

## merge all cell types pcHi-C links
# read them
neuron <- read.csv(paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/links/", dataset, "_neuron_pchic_links.tsv"), sep = "\t")
npc <- read.csv(paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/links/", dataset, "_NPC_pchic_links.tsv"), sep = "\t")
ipsc <- read.csv(paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/links/", dataset, "_iPSC_pchic_links.tsv"), sep = "\t")

# add cell type id 
neuron$celltype <- rep("Neuron", nrow(neuron))
npc$celltype <- rep("NPC", nrow(npc))
ipsc$celltype <- rep("iPSC", nrow(iPSC))

# combine them all
all <- rbind(neuron, npc, ipsc)

# save links
write.csv(ds.baits_ds.ehcs, paste0("/g/scb/zaugg/deuner/valdata/pcHi-C/links/", dataset, "_all_pchic_links.tsv"))

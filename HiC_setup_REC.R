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

# define dataset
dataset <- "timecourse"

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
View(ds.baits_ds.ehcs)

# save links into a csv file 
write.csv(ds.baits_ds.ehcs, paste0("/g/scb/zaugg/deuner/GRaNIE/validationdata/", dataset, "_", celltype, "_pchic_links.csv"))


#%%%%%#
# NPC # 
#%%%%%#

# This dataset does not contain a baits file and has two links files, a promoter-promoter interactions file and a promoter-other 
#  interactions files

celltype <- "NPC"

# load links files
pp.links <- read.csv("/g/scb/zaugg/deuner/valdata/pcHi-C/NPC/GSE86189_npc.pp.all.txt.bz2", sep= "\t")
po.links <- read.csv("/g/scb/zaugg/deuner/valdata/pcHi-C/NPC/GSE86189_npc.po.all.txt.bz2", sep = "\t")

# merge them 
links <- rbind(pp.links, po.links)

# create individual columns for chr, start, end
links[c("chr.1", "start.1", "end.1")] <- str_split_fixed(links$frag1, "[.]", 3)
links[c("chr.2", "start.2", "end.2")] <- str_split_fixed(links$frag2, "[.]", 3)
head(links)

# change format of frag1 and frag2 columns to chr:start-end
links$frag1 <- sub("[.]", ":", links$frag1) # replaces the first match
links$frag1 <- sub("[.]", "-", links$frag1)
links$frag2 <- sub("[.]", ":", links$frag2) # replaces the first match
links$frag2 <- sub("[.]", "-", links$frag2)
head(links)

# pcHiC read data was mapped to human genome hg19 so we need to adapt the hg19 annotations to hg38 reference genome
library(rtracklayer)
chain <- import.chain("/g/scb/zaugg/deuner/valdata/hg19ToHg38.over.chain")

# create ranges column 
links$ranges1 <- strsplit(links$frag1, ":", 2) %>% map(2)
links$ranges2 <- strsplit(links$frag2, ":", 2) %>% map(2)
head(links)

# save coordinates in separated files to convert them using liftOver webpage
write.csv(links$frag1, "/g/scb/zaugg/deuner/valdata/prom_cords.csv", row.names=FALSE)
write.csv(links$frag2, "/g/scb/zaugg/deuner/valdata/link_cords.csv", row.names=FALSE)

# create GRanges objects, 
links1 <- makeGRangesFromDataFrame(links, seqnames.field = "chr.1", start.field = "start.1", end.field = "end.1", keep.extra.columns = TRUE)   
links2 <- makeGRangesFromDataFrame(links, seqnames.field = "chr.2", start.field = "start.2", end.field = "end.2", keep.extra.columns = TRUE)

# adapt coordinates
l1 <- liftOver(links1, chain) #11016843

l1 <- unlist(l1) # 12563604
ll1 <- data.frame(l1) #12563604
nrow(ll1 %>% distinct())

# HI HA DUPLICATES EN LES COLUMNES FRAG1 I FRAG2, MIRAR ELS DUPLICATES I AGAFAR EL MIN DELS STARTS I EL MAX DELS ENDS!!!

class(links1[1][1,ranges()])
ranges(links1)[[1]][1]
liftover()  


  rename(seqnames = "chr.1", start = "start.1", end = "end.1") %>%
  dplyr::select(c("chr.1", "start.1", "end.1", "frag1", "frag2"))
links2 <- makeGRangesFromDataFrame(links2, keep.extra.columns = TRUE) %>%
  liftOver(chain) %>%
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
baits <- links %>% separate(frag1, c("chr", "start", "end"), "[.]")
baits <- baits %>% dplyr::select(c("chr", "start", "end", "frag2"))

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
ds.baits <- ds.baits %>% unite("frag1", 1:3, sep = ".")
head(ds.baits)

# look at all regions linked to these baits from links file, they essentially represent enhancers
# get links from baits found in timecourse data
ds.baits_ehcs <- ds.baits #frag2 represent the enhancers
ds.baits_ehcs <- ds.baits_ehcs %>% separate(frag2, c("chr", "start", "end"), "[.]")

# filter for these enhancers that are present in my peakset
## (expand peak +/- 1kb and check for any overlap with pchic enhancers)
ehcs <- data.frame(chr = ds.baits_ehcs$chr, start = ds.baits_ehcs$start, end = ds.baits_ehcs$end) 
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
View(ds.baits_ds.ehcs)

# save links into a csv file 
# IMPORTANT: 2 mandatory columns: peak and gene
write.csv(ds.baits_ds.ehcs, paste0("/g/scb/zaugg/deuner/GRaNIE/validationdata/", dataset, "_", celltype, "_pchic_links.csv"))

#%%%%%#
# PSC # 
#%%%%%#

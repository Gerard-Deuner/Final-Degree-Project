##############
# HI-C SETUP #
##############

# load libraries
library(dplyr)
library(Seurat)
library(qs)
library(GenomicRanges)

# pull the promoter coordinates for the genes in the timecourse dataset
## take the TSS and extend 2kb upstream if strand +, downstream if strand -
# load timecourse seurat object to get the gene names present in the dataset
timecourse.s <- qread("/g/scb/zaugg/deuner/GRaNIE/tmp/timecourse.pp.seuratObject.qs")
DefaultAssay(timecourse.s) <- "RNA"
tc.gene.names <- rownames(timecourse.s)

# get gene strand and position details from gtf file
gtf.file = "/g/scb/zaugg/marttine/RNA_ATAC_integration/annotation/gencode.v41.chr_patch_hapl_scaff.annotation.gtf.gz"
gtf <- rtracklayer::import(gtf.file)
gtf = as.data.frame(gtf)

# subset genes from tgf file that are present in timecourse dataset in bed format
tc.genes.ann <- gtf %>% filter((gene_name %in% tc.gene.names) & (type == "gene")) %>% select(seqnames, start, end, strand, gene_name)
colnames(tc.genes.ann)[1] <- "chr"

# get the promoter regions from each gene
prom.cords <- tc.genes.ann
prom.cords$start <- ifelse(prom.cords$strand == "+", tc.genes.ann$start-2000, tc.genes.ann$end)
prom.cords$end <- ifelse(prom.cords$strand == "+", tc.genes.ann$start, tc.genes.ann$end+2000)

# keep baits that intersect with promoter coordinates 
baits.file <- "/g/scb/zaugg/marttine/RNA_ATAC_integration/CRISPR/GRN/baits.csv"
baits <- read.csv(baits.file, row.names = 1)
# use findOverlaps method from GenomicRanges
qry <-  makeGRangesFromDataFrame(prom.cords, keep.extra.columns = TRUE) 
ref <-  makeGRangesFromDataFrame(baits)
ovlp <- findOverlaps(qry, ref) # returns indexes of intersecting regions
ovlp.indx <- ovlp@to
# subset baits overlapping timecourse gene promoters
tc.baits <-  baits[ovlp.indx,]
head(tc.baits)
# keep the gene names to which the baits map to
tc.baits <- tc.baits %>% cbind(gene = prom.cords$gene_name[ovlp@from])

# look at all regions linked to these baits from links file, they essentially represent enhancers
links.file <- "/g/scb/zaugg/marttine/RNA_ATAC_integration/CRISPR/GRN/pchic.i3Neuron.Song2019.hg38.csv"
links <- read.csv(links.file, sep = " ")
# get links from baits found in timecourse data
tc.baits_ehcs <- tc.baits %>% right_join(links, by = c("chr", "start", "end"), multiple = "all")

# filter for these enhancers that are present in my peakset
## (expand peak +/- 1kb and check for any overlap with pchic enhancers)
ehcs <- data.frame(chr = tc.baits_ehcs$chr.1, start = tc.baits_ehcs$start.1, end = tc.baits_ehcs$end.1) 
# get peaks present in timecourse data
DefaultAssay(timecourse.s) <- "ATAC"
peaks.vect <- rownames(timecourse.s)
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
tc.baits_tc.ehcs <- tc.baits_ehcs[ovlp.peaks.indx,]
tc.baits_tc.ehcs <- tc.baits_tc.ehcs %>% cbind(peak = peak.names)
View(tc.baits_tc.ehcs)

# save links
write.csv(tc.baits_tc.ehcs, "/g/scb/zaugg/deuner/GRaNIE/code/timecourse_pchic_links.csv")

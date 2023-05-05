############################
# eQTL ENRICHMENT ANALYSIS #
############################

# load libraries
library(dplyr)
library(Seurat)
library(qs)
library(GenomicRanges)
library(tidyr)
library(plyranges)
library(data.table)
library(purrr)
library(stringr)
library(ggplot2)

# read arguments from command line
args <- commandArgs(trailingOnly = TRUE)

# define dataset
dataset <- args[1] # timecourse | combined

# define correlation method
corr.method <- args[2] # pearson | spearman

# define nature of the eQTL data 
nature <- args[3] # metabrain | hipsci | all

# validate links against all significat peak-gene links or peak-gene links from filtered GRN 
links.to.validate <- args[4] # all | filtered

# get list of eqtl files
if (nature == "all"){
  
  file.names.brain <- list.files("/g/scb/zaugg/deuner/valdata/eQTL/metabrain/inputdata/")
  file.names.brain <- grep(".txt.gz", file.names.brain, value = TRUE)
  file.names.hipsci <- list.files("/g/scb/zaugg/deuner/valdata/eQTL/hipsci/inputdata/")
  file.names.hipsci <- grep(".txt.gz", file.names.hipsci, value = TRUE)
  file.names.hipsci <- grep("full", file.names.hipsci, value = TRUE)
  file.names <- c(file.names.brain, file.names.hipsci)

  } else {
  
  # set eQTL data directory
  data.dir <- paste0("/g/scb/zaugg/deuner/valdata/eQTL/", nature, "/inputdata/")
  # set output directory
  out.dir <- paste0("/g/scb/zaugg/deuner/valdata/eQTL/", nature, "/setup_outputs/")
  # list of files
  file.names <- list.files(data.dir) 
  # just take into account files containing the correct content
  file.names <- grep(".txt.gz", file.names, value = TRUE)

  }


# dataframe where all eQTLs are going to be stored
eqtl <- data.frame(matrix(ncol = 5, nrow = 0))
eqtl.names <- c("Gene", "SNP", "SNPChr", "SNPPos", "MetaP")
colnames(eqtl) <- eqtl.names

# convert all columns to proper class so bind_rows works
eqtl <- eqtl %>%
  mutate(Gene = as.character(Gene),
         SNP = as.character(SNP),
         SNPChr = as.integer(SNPChr),
         SNPPos = as.integer(SNPPos),
         MetaP = as.numeric(MetaP))

## store the eQTL data into a single dataframe
# iterate over the eQTL files
for (file in file.names){

  print(file)
  # read files one by one

  if (nature == "all"){
    # read MetaBrain files
    if (grepl("cortex", file, fixed = TRUE)){
      data <- fread(paste0("/g/scb/zaugg/deuner/valdata/eQTL/metabrain/inputdata/", file), sep = "\t")
    # read hiPSCs files
    } else {
    data <- fread(paste0("/g/scb/zaugg/deuner/valdata/eQTL/hipsci/inputdata/", file), sep = "\t")
    }
    
    # remove uninformative columns filter for the eQTL significance threshold to reduce the size
    col.names <- if (grepl("cortex", file, fixed = TRUE)){
      c("Gene", "SNP", "SNPChr", "SNPPos", "MetaP")
    } else {
      c("gene", "alt", "chr", "snp_pos", "pvalue")
    }
    
    # select important columns
    data <- data %>% dplyr::select(all_of(col.names))

  } else {
    data <- fread(paste0(data.dir, file), sep = "\t")
  }
  
  # rename column names (so they match with the global eQTL df's names)
  names(data) <- eqtl.names
  
  # convert all columns to proper class so bind_rows works
  data <- data %>%
    mutate(Gene = as.character(Gene),
           SNP = as.character(SNP),
           SNPChr = as.integer(SNPChr),
           SNPPos = as.integer(SNPPos),
           MetaP = as.numeric(MetaP))
  
  # concatenate the file with the rest of eQTL files
  eqtl <- bind_rows(list(eqtl, data)) #bind_rows much faster than rbind

}


head(eqtl)

# set threshold for significant eQTLs
eQTL.FDR <- 0.3

# number of eQTLs before filtering
print(paste("number of eQTLs before filtering", nrow(eqtl)))

# remove gene ensembl versions from the ensembl id
eqtl$Gene <- strsplit(eqtl$Gene, "[.]") %>%
  map(1) %>%
  unlist()

# make a genomic ranges list of significant eQTLs
gr.eqtl <- eqtl %>%
  dplyr::select(Gene, SNP, SNPChr, SNPPos, MetaP) %>%
  mutate(fdr = p.adjust(MetaP, method = "fdr")) %>%
  dplyr::filter(fdr < eQTL.FDR) %>%
  dplyr::rename(eqtl.gene = Gene) %>%
  makeGRangesFromDataFrame(., keep.extra.columns=T, seqnames.field = 'SNPChr', start.field = 'SNPPos', end.field = 'SNPPos')

# number of eQTLs after filtering
print(paste("number of eQTLs before filtering", length(gr.eqtl)))

# create df to store all links
all.links <- data.frame(matrix(nrow = 0, ncol = 4))
names(all.links) <- c("gene", "peak", "resolution", "validated")
# convert all columns to proper class so bind_rows works
all.links <- all.links %>%
  mutate(gene = as.character(gene),
         peak = as.character(peak),
         resolution = as.numeric(resolution),
         validated = as.integer(validated))

# set cluster resolutions to test
resolutions <- c(c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2)))

# iterate over GRNs
for (res in resolutions){
  
  # read GRN (NOT THE GRN ITSELF BUT THE LINKS TABLE)
  grn <- qread(paste0("/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/", dataset, "_batch_mode_", corr.method, "/Batch_Mode_Outputs/output_pseudobulk_clusterRes", res, "_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/GRN.qs"))
  
  # take desired links
  grn <- if(links.to.validate == "all"){
    grn@connections$peak_genes$`0`
    } else { 
      grn@connections$all.filtered$`0`
    }
  
  if (nrow(grn) > 0){
    grn <- as.data.frame(grn)
    names(grn)[5] <- "peak_gene.p_adj"
    
    # read genes and their positions from Ensembl:
    genes <- fread("/g/scb/zaugg/claringb/eQTL_overlap_GRN/input/ENSG_genes_biomart_GRCh38_20220519.txt")
    
    # read GRN object file to get genes that were used to create this GRN
    grn_genes <- qread(paste0("/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/", dataset, "_batch_mode_", corr.method, "/Batch_Mode_Outputs/output_pseudobulk_clusterRes", res, "_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/GRN.qs"))
    
    ## Filter GRN
    # Set an FDR threshold to include peak-gene links below a certain significance threshold.
    
    # Set threshold for significant GRN links
    GRN_FDR <- 0.05
    
    #Filter based on FDR and include only positive links
    grn <- grn %>%
      dplyr::filter(peak_gene.r > 0) %>%
      dplyr::filter(peak_gene.p_adj < GRN_FDR)
    
    # filter for GRN enhancers
    # Include the GRN enhancers that overlap with at least one significant eQTL SNP
    # keep distinct enhancers in this GRN (GRanges object)
    gr.grn_enhancers <- grn %>%
      dplyr::select(peak.ID) %>%
      tidyr::separate(peak.ID, into = c("peak.chr", "peak.start", "peak.end"), remove = F) %>%
      mutate(peak.chr = str_replace_all(peak.chr, "chr", ""),
             peak.start.pos = peak.start) %>%
      distinct() %>%  
      makeGRangesFromDataFrame(., keep.extra.columns=T, seqnames.field = 'peak.chr', start.field = 'peak.start', end.field = 'peak.end')
    
    # overlap GRN enhancers with eQTLs to keep enhancers that contain at least one significant eQTL SNP
    gr.filt_enhancers <- join_overlap_intersect(gr.eqtl, gr.grn_enhancers)
    
    # make into data frame
    filt_enhancers <- as.data.frame(gr.filt_enhancers) %>%
      mutate(p_beta = as.character(MetaP),
             fdr = as.numeric(fdr))
    
    nr_bf_filt <- length(unique(gr.grn_enhancers$peak.ID))
    nr_aft_fil <- length(unique(filt_enhancers$peak.ID))
    
    # Number of enhancers in GRN before filtering:
    nr_bf_filt
    # Number of enhancers in GRN after filtering:
    nr_aft_fil
    
    # Validate enhancer-gene links based on the GRN
    # Keep distinct enhancer-gene links from the GRN and test how many of the GRN enhancer-gene links overlap with any eQTLs by annotating 
    #  each link as `TRUE` (if any eQTL SNP is located in the enhancer and its target gene matches) or `FALSE` 
    #  (eQTL SNP is located in the enhancer but none of its target genes overlaps with the GRN one)
    # select enhancers - target gene links from the GRN
    grn <- grn %>%
      #select(starts_with("peak") | starts_with("gene")) %>%
      mutate(peak_gene.distance_bin = case_when(peak_gene.distance <= 50000 ~ "0-50kb",
                                                peak_gene.distance > 50000 & peak_gene.distance <= 100000 ~ "50-100kb",
                                                peak_gene.distance > 100000 & peak_gene.distance <= 150000 ~ "100-150kb",
                                                peak_gene.distance > 150000 & peak_gene.distance <= 200000 ~ "150-200kb",
                                                peak_gene.distance > 200000 & peak_gene.distance <= 250000 ~ "200-250kb",
                                                peak_gene.distance > 250000 ~ ">250kb")) %>%
      mutate(peak_gene.distance_bin = factor(peak_gene.distance_bin, levels = c("0-50kb", "50-100kb", "100-150kb", "150-200kb", "200-250kb"))) %>%
      dplyr::select(peak.ID, gene.ENSEMBL, peak_gene.distance, peak_gene.distance_bin) %>%
      distinct() %>%
      dplyr::filter(peak.ID %in% filt_enhancers$peak.ID)
    
    # Number of enhancer-gene links in GRN:
    nrow(grn)
    
    # classify links as positive or negative
    grn_links_validated <- grn %>%
      left_join(filt_enhancers, by = "peak.ID", multiple = "all") %>%
      dplyr::rename(grn.gene = gene.ENSEMBL) %>%
      dplyr::select(peak.ID, grn.gene, eqtl.gene, SNP, fdr) %>%
      group_by(peak.ID, grn.gene) %>%
      mutate(link_validate = ifelse(grn.gene == eqtl.gene, TRUE, FALSE))  # case_when(any(grn.gene == eqtl.gene, na.rm = T) ~ TRUE, TRUE ~ FALSE))
    
    # show the GRN links that can be tested and whether they are validated or not
    data.table(grn_links_validated, rownames = F, filter = "top", options = list(pageLength = 5, scrollX = T))
    
    # format validated links to be stored
    links <- grn_links_validated %>%
      tibble::add_column(resolution = res) %>% 
      dplyr::select(gene = grn.gene, peak = peak.ID, resolution, validated = link_validate) %>%
      distinct()
    
    # format links so it has the same classes as all.links
    links <- links %>%
      mutate(gene = as.character(gene),
             peak = as.character(peak),
             resolution = as.numeric(resolution),
             validated = as.integer(validated))
    
    # add links data for this res to the global links df
    all.links <- bind_rows(all.links, links)  
    
    #----------------------------------------------------------
    
    # # Background validation #
    # 
    # # Validate enhancer-gene links based on the closest gene
    # # As a background, produce enhancer-gene links for the same enhancers by selecting the closest gene to the peak location.
    # # Test how many of the GRN enhancer-gene links overlap with eQTLs by annotating each link as `TRUE` (if any eQTL SNP is located in
    # #  the enhancer and its target gene matches) or `FALSE` (eQTL SNP is located in the enhancer but none of its target genes overlaps with the GRN one).
    # 
    # # select only genes that were considered in making this GRN
    # genes_in_grn <- grn_genes@connections$peak_genes$`0` %>%
    #   dplyr::select(gene.ENSEMBL) %>%
    #   inner_join(genes, by = c("gene.ENSEMBL" = "ensembl_gene_id"))
    # 
    # # make ensembl genes into Genomic Ranges dataset
    # gr.genes <- genes_in_grn %>%
    #   mutate(gene.start = start_position) %>%
    #   makeGRangesFromDataFrame(., keep.extra.columns=T, seqnames.field = 'chromosome_name', start.field = 'start_position', end.field = 'end_position')
    # 
    # # find the nearest gene for each enhancer
    # # make into dataframe
    # # classify links as positive or negative
    # nearest_links_validated <- join_nearest(gr.filt_enhancers, gr.genes) %>%
    #   as.data.frame() %>%
    #   dplyr::rename(nearest.gene = gene.ENSEMBL) %>%
    #   dplyr::select(peak.ID, nearest.gene, eqtl.gene, SNP, fdr) %>%
    #   group_by(peak.ID, nearest.gene) %>%
    #   mutate(link_validate = case_when(any(nearest.gene == eqtl.gene, na.rm = T) ~ TRUE,
    #                                    TRUE ~ FALSE))
    # 
    # # Number of enhancer-gene links to be tested
    # nr_links <- nearest_links_validated %>% dplyr::select(peak.ID, nearest.gene) %>% distinct() %>% nrow()
    # 
    # # Number of enhancer-gene links based on nearest gene:
    # nr_links
    # 
    # # show the nearest gene links that can be tested and whether they are validated or not
    # data.table(nearest_links_validated, rownames = F, filter = "top", options = list(pageLength = 5, scrollX = T))
    # 
    # 
    # # Validate enhancer-gene links based on the randomly sampled distance-matched genes
    # # As a background, produce enhancer-gene links for the same enhancers by selecting a gene with the same distance to the peak location.
    # # Test how many of the GRN enhancer-gene links overlap with eQTLs by annotating each link as `TRUE` (if any eQTL SNP is located in the enhancer and its target gene matches) or `FALSE` (eQTL SNP is located in the enhancer but none of its target genes overlaps with the GRN one).
    # # Expand ranges for enhancers to 250kb on either side
    # # Use the gr.grn_enhancers GRanges object for this because that includes the enhancer info as ranges
    # gr.grn_enhancer_windows <- resize(gr.grn_enhancers, width = width(gr.grn_enhancers)+(500000*2), fix = "start")
    # 
    # # Intersect regions around enhancers with all genes
    # gr.grn_enhancer_windows <- join_overlap_inner(gr.grn_enhancer_windows, gr.genes) #keep all genes that fall within the extended enhancer regions
    # 
    # # Filter for enhancers that overlap with an eQTL SNP
    # gr.filt_enhancer_windows <- gr.grn_enhancer_windows[which(elementMetadata(gr.grn_enhancer_windows)[,1] %in% gr.filt_enhancers$peak.ID)]
    # 
    # # Add enhancer - gene distance bin
    # gr.filt_enhancer_windows$peak_gene.distance <- abs(as.numeric(gr.filt_enhancer_windows$peak.start.pos) - gr.filt_enhancer_windows$gene.start)
    # gr.filt_enhancer_windows$peak_gene.distance_bin <- cut(gr.filt_enhancer_windows$peak_gene.distance,
    #                                                        breaks = c(0, 50000, 100000, 150000, 200000, 250000, 500000),
    #                                                        labels = c("0-50kb", "50-100kb", "100-150kb", "150-200kb", "200-250kb", ">250kb"))
    # 
    # # Calculate number of observations that need to be sampled from each bin size
    # sample_freq <- grn %>%
    #   group_by(peak_gene.distance_bin) %>%
    #   tally(name = 'freq')
    # 
    # # Add frequency metadata #NB GRN had quite a few gene-peak links with >250kb, which are removed at this step since the sampling should reflect the enhancers
    # filt_enhancer_windows <- merge(gr.filt_enhancer_windows, sample_freq, by = "peak_gene.distance_bin")
    # 
    # # Sample using the frequency rate
    # # Add eQTL information
    # # classify links as positive or negative
    # 
    # #Replicate the sampling 20 times
    # random_links_validated_all <- lapply(1:20, function(rep) {
    #   group_by(filt_enhancer_windows, peak_gene.distance_bin) %>%
    #     sample_n(freq[1]) %>%
    #     ungroup() %>%
    #     inner_join(filt_enhancers, by = "peak.ID") %>%
    #     dplyr::rename(random.gene = gene.ENSEMBL) %>%
    #     dplyr::select(peak.ID, random.gene, eqtl.gene, SNP, fdr) %>%
    #     #  select(peak.ID, random.gene, eqtl.gene, variant, fdr, peak_gene.distance) %>%
    #     group_by(peak.ID, random.gene) %>%
    #     mutate(link_validate = case_when(any(random.gene == eqtl.gene, na.rm = T) ~ TRUE,
    #                                      TRUE ~ FALSE))
    # })
    # 
    # #Number of enhancer-gene links to be tested
    # nr_links_rand <- random_links_validated_all[[1]] %>% dplyr::select(peak.ID, random.gene) %>% distinct() %>% nrow()
    # 
    # # Make figures output side by side
    # par(mar = c(4, 4, .1, .1))
    # 
    # #Quick sanity check if the gene-distance sampling went okay
    # # Removed it after implementing re-sampling to avoid having to add peak_gene.distance column
    # ggplot(grn) +
    #   geom_histogram(aes(x=peak_gene.distance), fill = "white", col = "black", binwidth = 25000) +
    #   theme_classic() +
    #   ggtitle("Peak gene distance of GRN enhancer-gene pairs") +
    #   xlab("Peak gene distance")
    # 
    # ggplot(random_links_validated_all) +
    #   geom_histogram(aes(x=peak_gene.distance), fill = "white", col = "black", binwidth = 25000) +
    #   theme_classic() +
    #   ggtitle("Peak gene distance of subsampled enhancer-gene pairs") +
    #   xlab("Peak gene distance")
    # 
    # 
    # # Number of enhancer-gene links based on random gene:
    # nr_links_rand
    # 
    # # Compare validation of GRN links with validation of nearest gene links
    # # Sum up positive and negative links GRN:
    # # function to make a table of the results
    # make_table <- function(df){
    #   colnames(df) <- c("peak.ID", "gene", "eqtl.gene", "SNP", "fdr", "link_validate")
    #   table <- df %>%
    #     dplyr::select(peak.ID, gene, link_validate) %>%
    #     distinct() %>%
    #     group_by(link_validate) %>%
    #     tally() %>%
    #     mutate(percentage = 100*n/sum(n))
    #   return(table)
    # }
    # 
    # t_grn <- make_table(grn_links_validated)
    # t_grn
    # 
    # # Sum up positive and negative links nearest gene:
    # t_near <- make_table(nearest_links_validated)
    # t_near
    # 
    # # Sum up positive and negative links random distance-matched gene (one example):
    # t_rand <- lapply(random_links_validated_all, make_table)
    # t_rand[1]
    # 
    # # Is the percentage of GRN enhancer-gene links validated with an eQTL higher than enhancer-nearest gene links?
    # t_near %>% filter(link_validate == TRUE) %>% dplyr::select(percentage) < t_grn %>% filter(link_validate == TRUE) %>% dplyr::select(percentage)
    # 
    # # Is the percentage of GRN enhancer-gene links validated with an eQTL higher than enhancer-gene links with a similar distance?
    # t_rand[[5]] %>% filter(link_validate == TRUE) %>% select(percentage) < t_grn %>% filter(link_validate == TRUE) %>% select(percentage)
    # 
    # odds_grn <- t_grn[t_grn$link_validate == T,]$n / t_grn[t_grn$link_validate == F,]$n
    # odds_rand <- list()
    # or <- list()
    # 
    # for (i in 1:20){
    #   odds_rand[[i]] <- t_rand[[i]][t_rand[[i]]$link_validate == T,]$n / t_rand[[i]][t_rand[[i]]$link_validate == F,]$n
    #   or[[i]] <- odds_grn / odds_rand[[i]]
    # }
    # 
    # mean_or <- mean(unlist(or))
    # min_or <- min(unlist(or))
    # max_or <- max(unlist(or))
    # 
    # # filename <- paste0("/g/scb/zaugg/deuner/valdata/eQTL/", nature, "/setup_outputs/or/",  dataset, "_", corr.method, "-FDR", eQTL.FDR, "_in_",
    # #                    res, ".GRN", "-peak-gene-FDR", GRN_FDR,".txt")
    # # write.table(unlist(or), filename, sep = "\t", quote = F, col.names = F, row.names = F)
  }
}

# write file (and create a dataset and corr.method -specific folder if it does not exist yet)
write.csv(all.links, paste0("/g/scb/zaugg/deuner/valdata/eQTL/validated_links/", dataset, "_", corr.method, "_", links.to.validate, "_eQTL_links.tsv"))
  
# The odds ratio of the GRN links over the random distance-matched links being validated by eQTLs is:
round(mean_or, digits = 3) 
(round(min_or, digits = 3) - round(max_or, digits = 3))



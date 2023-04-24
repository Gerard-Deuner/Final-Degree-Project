#######################
# ChiP-Seq Validation #
#######################

# load libraries
library(dplyr)
library(qs)
library(Seurat)
library(ggplot2)
library(pROC)

# define datasets
datasets <- c("timecourse", "combined")

# define correlation methods
corr.methods <- c("pearson", "spearman")

# define GRN method
#methods <- c("GRaNIE", "SCENIC+", "Pando")

# set vector with the resolutions tested
resolutions <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))

# dataframes to store TPs and FPs
TPs <- data.frame(matrix(ncol = length(datasets) * length(corr.methods), nrow = length(resolutions)), row.names = resolutions)
FPs <- data.frame(matrix(ncol = length(datasets) * length(corr.methods), nrow = length(resolutions)), row.names = resolutions)
it <- 1
for (dataset in datasets){
  for (corr.method in corr.methods){
    names(TPs)[it] <- paste(dataset, corr.method, sep = "_")
    names(FPs)[it] <- paste(dataset, corr.method, sep = "_")
    it <- it + 1
  }
}

# list of dfs where all links are going to be stored independent of resolution
TP.all.list <- list() 

# iterate over datasets
for (dataset in datasets){
  # iterate over correlation methods
  for (corr.method in corr.methods){
    # set output dir
    out.dir <- paste0("/g/scb/zaugg/deuner/valdata/ChiP-seq/results/", dataset, "_", corr.method, "/")
    dir.create(out.dir) # create dir if it doesn't exist already
    
    # path of GRNs
    GRNs.dir <- paste0("/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/", dataset, "_batch_mode_", corr.method, "/Batch_Mode_Outputs/")
    
    # ChiP-Seq links path
    links.dir <- paste0("/g/scb/zaugg/deuner/valdata/ChiP-seq/", "filtered_chipseq_", dataset,".tsv.gz")
    links <- read.csv(links.dir) 
    
    # subset just useful columns and unique links
    links <- links %>%
      dplyr::select(TF, peak) %>%
      distinct(TF, peak, .keep_all = T)
    head(links)
    
    # set colors for ROC curves
    colors <- c("#FFC107", "#2196F3", "#4CAF50", "#FF5722", "#9C27B0", "#3F51B5", "#FFEB3B", "#8BC34A", "#E91E63", 
                         "#673AB7", "#00BCD4", "#FF9800", "#CDDC39", "#795548", "#607D8B", "#9E9E9E", "#F44336", "#00ACC1")
                         
    # set vecor of auc positioning in y axis
    aucy <- seq(10, 110, by = 5)
    
    # iterate over GRNs
    j = 0
    
    # boolean for first plot
    first = TRUE
    
    # dataframe where all peak-gene links will be stored (independent of the resolution)
    TP.all <- data.frame()
    
    # vectors of number of TPs and FPs per dataset and corr.method
    TP.vec <- c()
    FP.vec <- c()
    
    for (res in resolutions){
      # set index
      j <- j + 1
      
      # set up GRN directory
      GRN.dir <- paste0(GRNs.dir, "output_pseudobulk_clusterRes", res, "_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/")
      # read GRN object
      GRN <- qread(paste0(GRN.dir, "GRN.qs"))
      
      # get TF-peak connections from GRN
      GRN.links <- GRN@connections$TF_peaks$`0`$main %>% dplyr::filter(TF_peak.fdr < 0.2) %>% as.data.frame() %>% dplyr::select(peak.ID, TF.ID, TF_peak.r, TF_peak.fdr, TF_peak.fdr_direction)
      
      # adapt peaks format to links format
      form.peak <- rep("", nrow(GRN.links))
      for (i in 1:nrow(GRN.links)){
        peak <- as.character(GRN.links$peak.ID[i])
        split.peak <- strsplit(peak, split = ":")
        new.peak <- paste(split.peak[[1]][1], split.peak[[1]][2], sep = "-")
        form.peak[i] <- new.peak
      }
      GRN.links$id <- seq(1:nrow(GRN.links))
      GRN.links$peak.ID <- form.peak 
      names(GRN.links)[1:2] <- c("peak", "TF")
      # just keep the TF name
      GRN.links$TF <- as.character(GRN.links$TF)
      GRN.links$TF <- strsplit(GRN.links$TF, "[.]") %>% map(1) %>% as.character()
      
      # Binary classification
      # get rows that match in GRN links and HiC links by gene name and peak id
      TP <- inner_join(GRN.links, links, by = c("peak", "TF"), multiple = "all")
      TP$res <- rep(resolutions[j], nrow(TP))
      TP.all <- rbind(TP.all, TP)
      # set 1 if they are in ChiP-seq data and 0 if not
      bc <- ifelse(GRN.links$id %in% TP$id, yes = 1, no = 0)
      # add column to GRN_links with that information
      GRN.links$bc <- bc
      # create an confusion matrix 
      conf.m <- table(bc)
      
      # store GRN links data in a variable
      assign(paste("GRN_links", res, sep = "."), GRN.links)
      
      # Compute ROC curve
      par(pty = "s")
      
      # print metrics
      print("")
      print(paste0("Resolution: ", res))
      print(paste0("GRN links: ", nrow(GRN.links)))
      print(paste0("TP: ", nrow(TP)))
      FP <- nrow(GRN.links) - nrow(TP)
      print(paste0("FP: ", FP))
      print("")
      
      # Store TPs and FPs
      TP.vec <- c(TP.vec, nrow(TP))
      FP.vec <- c(FP.vec, FP)
      
      # if ((first == TRUE) & (sum(bc) > 0)){
      #   ROC <- roc(GRN.links$bc, GRN.links$TF_peak.r, plot = TRUE, legacy.axes = TRUE, percent = TRUE,
      #              xlab="False Positive Percentage | 1 - Specificity", ylab="True Positive Percentage | Sensitivity",
      #              col = colors[j], lwd=2, print.auc = TRUE, main = "ROC", add = FALSE, print.auc.x = 115, print.auc.y = aucy[j])
      #   legend("bottomright",  legend = resolutions, col = colors, ncol = 2)
      #   first = FALSE
      #   print(ci.auc(ROC))
      # }
      # if ((first == FALSE) & (sum(bc) > 0)){ 
      #   ROC <- roc(GRN.links$bc, GRN.links$TF_peak.r, plot = TRUE, legacy.axes = TRUE, percent = TRUE,
      #              xlab="False Positive Percentage | 1 - Specificity", ylab="True Positive Percentage | Sensitivity",
      #              col = colors[j], lwd=2, print.auc = TRUE, main = "ROC", add = TRUE, print.auc.x = 115, print.auc.y = aucy[j])
      #   legend("bottomright", legend = resolutions, col = colors, lwd = 4, ncol = 2)
      #   print(ci.auc(ROC))
      # }
      
      
    }
    # save TPs and FPs data
    print(TP.vec)
    print(FP.vec)
    TPs[, paste(dataset, corr.method, sep = "_")] <- TP.vec
    FPs[, paste(dataset, corr.method, sep = "_")] <- FP.vec
    TP.all.list <- append(TP.all.list, TP.all)
  }
}

TPs$resolution <- as.factor(resolutions)

# TF recovery plot
ggplot(TPs, aes(resolution)) + 
  geom_line(aes(y = timecourse_pearson, group = 1, colour = "timecourse_pearson", linetype = "timecourse_pearson")) + 
  geom_line(aes(y = timecourse_spearman, group = 1, colour = "timecourse_spearman", linetype = "timecourse_spearman")) + 
  geom_line(aes(y = timecourse_spearman, group = 1, colour = "combined_pearson", linetype = "combined_pearson")) + 
  geom_line(aes(y = timecourse_spearman, group = 1, colour = "combined_spearman", linetype = "combined_spearman")) + 
  labs(title = "TF Recovery", x = "Resolutions", y = "# Recovered TFs", colour = "Dataset", linetype = "Correlation Method") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_bw()
  
  
  

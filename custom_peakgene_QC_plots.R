###########################################################
# CUSTOM PEAK-GENE QUALITY CONTROL PLOTS (MULTIPLE EGRNS) #
###########################################################

library(GRaNIE)
library(dplyr)
library(Seurat)
library(qs)

# DISCLAIMER: This is a modification of GRaNIE original code

# List of GRNs
# /g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/combined_batch_mode_spearman_nomicro/Batch_Mode_Outputs/output_pseudobulk_clusterRes0.1_RNA_limma_quantile_ATAC_DESeq2_sizeFactors

resolutions <- c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))#c(20, 18, 16, 14, 12, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0.75, 0.5, 0.25)

density_plots <- c()
ratio_plots <- c()
rcorr_plots <- c()

ii <- 1
for (res in resolutions){
  
  GRN <- qread(paste0("/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/combined_batch_mode_spearman_nomicro/Batch_Mode_Outputs/output_pseudobulk_clusterRes", res, "_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/GRN.qs"))
  
  gene.types <- list(c("all")) #, c("protein_coding"))
  useFiltered <- FALSE
  plotDetails <- FALSE
  plotPerTF <- FALSE
  fileBase <- NULL
  pdf_width <- 12
  pdf_height <- 12
  pages <- NULL
  forceRerun <- FALSE
  
  
  cols_keep = c("peak.ID", "gene.ENSEMBL", "peak_gene.r", "peak_gene.p_raw", "peak_gene.distance")
  
  # Change from 10 to 5 
  nCategoriesBinning = 5
  probs = seq(0,1,1/nCategoriesBinning)
  
  
  range = GRN@config$parameters$promoterRange
  class_levels = c(paste0("real_",range), paste0("random_",range))
  networkType_details = c(paste0("real_",range), paste0("random_",range))
  
  colors_vec = c("black", "darkgray")
  networkType_vec = c("real", "background")
  names(colors_vec) = names(networkType_vec) = networkType_details
  
  options(dplyr.summarise.inform = FALSE) 
  
  includeRobustCols = FALSE
  #Either take all or the filtered set of connections
  if (useFiltered) {
    
    # TODO: Robust columns filter
    peakGeneCorrelations.all = 
      rbind(dplyr::select(GRN@connections$all.filtered[["0"]], dplyr::everything()) %>%
              dplyr::mutate(class = paste0("real_",range))) %>%
      rbind(dplyr::select(GRN@connections$all.filtered[["1"]],  dplyr::everything()) %>% 
              dplyr::mutate(class = paste0("random_",range)))
    
  } else {
    
    robustColumns = c("peak_gene.p_raw.robust", "peak_gene.bias_M_p.raw", "peak_gene.bias_LS_p.raw", "peak_gene.r_robust")
    if (all(robustColumns %in% colnames(GRN@connections$peak_genes[["0"]]))) {
      includeRobustCols = TRUE
      cols_keep = c(cols_keep, robustColumns)
    }
    
    
    
    if (!"peak.GC.perc" %in% colnames(GRN@annotation$peaks)) {
      GRN@annotation$peaks$peak.GC.perc = NA
    }
    
    peakGeneCorrelations.all = 
      rbind(
        dplyr::select(GRN@connections$peak_genes[["0"]], tidyselect::all_of(cols_keep)) %>%
          dplyr::mutate(class = factor(paste0("real_",range), levels = class_levels)),
        dplyr::select(GRN@connections$peak_genes[["1"]],  tidyselect::all_of(cols_keep)) %>% 
          dplyr::mutate(class = factor(paste0("random_",range), levels = class_levels))) %>%
      dplyr::left_join(dplyr::select(GRN@annotation$genes, "gene.ENSEMBL", "gene.type", 
                                     "gene.mean", "gene.median", "gene.CV"), by = "gene.ENSEMBL", multiple = "all") %>%
      dplyr::left_join(GRN@annotation$peaks %>% 
                         dplyr::select(-dplyr::starts_with("peak.gene."), -"peak.GC.perc"), by = "peak.ID") %>%
      dplyr::select(-"gene.ENSEMBL")
    
  }
  
  if (!plotPerTF) {
    peakGeneCorrelations.all = dplyr::select(peakGeneCorrelations.all, -"peak.ID")
  }
  
  nClasses_distance  = 10
  
  peakGeneCorrelations.all = peakGeneCorrelations.all %>%
    dplyr::mutate(r_positive = .data$peak_gene.r > 0,
                  peak_gene.distance_class = 
                    forcats::fct_explicit_na(addNA(cut(.data$peak_gene.distance, breaks = nClasses_distance, include.lowest = TRUE)), "random"),
                  peak_gene.distance_class_abs = forcats::fct_explicit_na(addNA(cut(abs(.data$peak_gene.distance), 
                                                                                    breaks = nClasses_distance, include.lowest = TRUE, ordered_result = TRUE)), "random"),
                  peak_gene.p.raw.class = cut(.data$peak_gene.p_raw, breaks = seq(0,1,0.05), include.lowest = TRUE, ordered_result = TRUE),
                  peak_gene.r.class = cut(.data$peak_gene.r, breaks = seq(-1,1,0.05), include.lowest = TRUE, ordered_result = TRUE)) %>%
    dplyr::filter(!is.na(.data$peak_gene.r)) # Eliminate rows with NA for peak_gene.r. This can happen if the normalized gene counts are identical across ALL samples, and due to the lack of any variation, the correlation cannot be computed
  
  
  # Oddity of cut: When breaks is specified as a single number, the range of the data is divided into breaks pieces of equal length, and then the outer limits are moved away by 0.1% of the range to ensure that the extreme values both fall within the break intervals. 
  levels(peakGeneCorrelations.all$peak_gene.distance_class_abs)[1] = 
    gsub("(-\\d+)", "0", levels(peakGeneCorrelations.all$peak_gene.distance_class_abs)[1], perl = TRUE)
  
  
  if (includeRobustCols) {
    peakGeneCorrelations.all = peakGeneCorrelations.all %>%
      dplyr::mutate(peak_gene.p_raw.robust.class = 
                      cut(.data$peak_gene.p_raw.robust, breaks = seq(0,1,0.05), include.lowest = TRUE, ordered_result = TRUE))
  }
  
  
  # Prepare plots #
  #informative_palette <- colorRampPalette(c("#ff0000", "#0000ff"))
  #informative_colors <- informative_palette(19)
  informative_colors <- c("#FF0000", "#FF3300", "#FF6600", "#FF9900", "#FFCC00", "#FFFF00", "#CCFF00", "#99FF00", "#66FF00", "#33FF00", "#00FF00", "#00FF33", "#00FF66", "#00FF99", "#00FFCC", "#00FFFF", "#00CCFF", "#0099FF", "#0066FF")
  background_colors <- c("#FFD9D9", "#FFCEC0", "#FFC9AF", "#FFC3AD", "#FFBD9B", "#FFF4AD", "#D9FFAD", "#C4FFAD", "#A9FFAD", "#8EFFAD", "#ADFFAD", "#ADFFCE", "#ADFFD9", "#ADFFE6", "#ADFFEF", "#ADEDFF", "#ADDEFF", "#ADD3EB", "#ADCEE6")
  
  # Background Color Gradient
  #background_palette <- colorRampPalette(c("#000000", "#999999"))
  #background_colors <- background_palette(19)
  
  colors_class = c("black", "black")
  names(colors_class) = unique(peakGeneCorrelations.all$class)
  colors_class[which(grepl("random", names(colors_class)))] = "darkgray"
    
  r_pos_class = c(informative_colors[[ii]], background_colors[[ii]]) #c("black", "darkgray")
  names(r_pos_class) = c("TRUE", "FALSE")
  
  dist_class = c(informative_colors[[ii]], background_colors[[ii]]) #c("dark red", "#fc9c9c") 
  names(dist_class) = class_levels
  
  freqs = table(peakGeneCorrelations.all$class)
  freq_class = paste0(gsub(names(freqs), pattern = "(.+)(_.*)", replacement = "\\1"), " (n=", .prettyNum(freqs) , ")")
  # Change upstream and go with "background" everywhere
  freq_class = gsub(freq_class, pattern = "random", replacement = "background")
  names(freq_class) <- names(freqs)
  
  xlabels_peakGene_r.class = levels(peakGeneCorrelations.all$peak_gene.r.class)
  nCur = length(xlabels_peakGene_r.class)
  xlabels_peakGene_r.class[setdiff(seq_len(nCur), c(1, floor(nCur/2), nCur))] <- ""
  
  # For the last plot, which is wider, we label a few more
  xlabels_peakGene_r.class2 = levels(peakGeneCorrelations.all$peak_gene.r.class)
  nCur = length(xlabels_peakGene_r.class2)
  xlabels_peakGene_r.class2[setdiff(seq_len(nCur), c(1, floor(nCur/4), floor(nCur/2), floor(nCur/4*3), nCur))] <- ""
  
  xlabels_peakGene_praw.class = levels(peakGeneCorrelations.all$peak_gene.p.raw.class)
  nCur = length(xlabels_peakGene_praw.class)
  xlabels_peakGene_praw.class[setdiff(seq_len(nCur), c(1, floor(nCur/2), nCur))] <- ""
  
  #
  # ITERATE THROUGH ALL GENE TYPES, WITH ONE PER PLOT
  #
  for (geneTypesSelected in gene.types) {
    
    # Reset page counter for each PDF anew
    pageCounter = 1
    
    futile.logger::flog.info(paste0(" Gene type ", paste0(geneTypesSelected, collapse = "+")))
    
    if ("all" %in% geneTypesSelected) {
      indexCur = seq_len(nrow(peakGeneCorrelations.all))
    } else {
      indexCur = which(peakGeneCorrelations.all$gene.type %in% geneTypesSelected)
    }
    
    
    # START PLOTTING #
    
    if (!is.null(fileBase)) {
      filenameCur = paste0(fileBase, paste0(geneTypesSelected, collapse = "+"), filteredStr, ".pdf")
      .checkOutputFile(filenameCur)
      grDevices::pdf(file = filenameCur, width = pdf_width, height = pdf_height)
      
      futile.logger::flog.info(paste0(" Plotting to file ", filenameCur))
      
    }
    
    if (plotPerTF) {
      
      TF.nRows = rep(-1, length(GRN@config$allTF))
      #TF.nRows = rep(-1, 10)
      TF.peaks = list()
      names(TF.nRows) = GRN@config$allTF
      for (TFCur in GRN@config$allTF) {
        TF.peaks[[TFCur]] = names(which(GRN@data$TFs$TF_peak_overlap[,TFCur] == 1))
        TF.nRows[TFCur] = peakGeneCorrelations.all[indexCur,] %>% dplyr::filter(.data$peak.ID %in% TF.peaks[[TFCur]]) %>% nrow()
      }
      
      TFs_sorted = names(sort(TF.nRows, decreasing = TRUE))
      
      allTF = c("all", TFs_sorted)
      
    } else {
      allTF = "all"
    }
    
    counter = 0
    for (TFCur in allTF) {
      
      counter = counter + 1
      
      if (is.null(pages) | (!is.null(pages) && pageCounter %in% pages)) {
        
        if (length(allTF) > 1) {
          futile.logger::flog.info(paste0(" QC plots for TF ", TFCur, " (", counter, " of ", length(allTF), ")"))
        }
        
        
        if ("all" %in% geneTypesSelected) {
          indexCur = seq_len(nrow(peakGeneCorrelations.all))
        } else {
          indexCur = which(peakGeneCorrelations.all$gene.type %in% geneTypesSelected)
        }
        
        
        if (TFCur != "all") {
          indexCur = intersect(indexCur, which(peakGeneCorrelations.all$peak.ID %in% TF.peaks[[TFCur]]))
        }
        
        if (length(indexCur) == 0) {
          message = " No connections left after filtering and intersection, skip plots"
          .checkAndLogWarningsAndErrors(NULL, message, isWarning = TRUE)
          next
        }
        
        # Get subset also for just the real data
        indexCurReal = intersect(indexCur, which(peakGeneCorrelations.all$class == names(dist_class)[1]))
        
        
        xlabel = paste0("Correlation raw p-value")
        
        # DENSITY PLOTS P VALUE
        
        # TODO: Densities as ratio
        # https://stackoverflow.com/questions/58629959/how-can-i-extract-data-from-a-kernel-density-function-in-r-for-many-samples-at-o#:~:text=To%20compute%20the%20density%20you,can%20use%20the%20package%20spatstat%20.
        
        # Produce the labels for the class-specific subtitles
        customLabel_class = .customLabeler(table(peakGeneCorrelations.all[indexCur,]$class))
        
        r_pos_freq = table(peakGeneCorrelations.all[indexCur,]$r_positive)
        labeler_r_pos = ggplot2::labeller(r_positive = c("TRUE"  = paste0("r positive (", r_pos_freq["TRUE"], ")"), 
                                                         "FALSE" = paste0("r negative (", r_pos_freq["FALSE"], ")")) )
        theme_main = ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, size = ggplot2::rel(0.8)),
                                    axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.8)),
                                    panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
        
        # dist_class2 <- dist_class
        # names(dist_class2) <- as.logical(c("TRUE", "FALSE"))

        # peakGeneCorrelations.all2 <- peakGeneCorrelations.all2 %>%
        #   mutate(r_pos_class = r_pos_class[peakGeneCorrelations.all2$r_positive])
        # 
        peakGeneCorrelations.all$r_pos_class = unname(r_pos_class[as.character(peakGeneCorrelations.all$r_positive)])
        peakGeneCorrelations.all$dist_class = unname(dist_class[as.character(peakGeneCorrelations.all$class)])
        
        names(peakGeneCorrelations.all$r_pos_class)
        
        ## p-val density curves stratified by real vs background ##
        
        if (res == 0.1) {
          gA2 <- ggplot2::ggplot(peakGeneCorrelations.all[indexCur,], ggplot2::aes(peak_gene.p_raw)) +
            ggplot2::geom_density(ggplot2::aes(color = dist_class, linetype = r_pos_class)) +
            ggplot2::facet_wrap(~ class, labeller = ggplot2::labeller(class = freq_class)) +
            ggplot2::xlab(xlabel) +
            ggplot2::ylab("Density") +
            ggplot2::theme_bw() +
            ggplot2::theme(legend.position = "none", 
                           axis.text = ggplot2::element_text(size = ggplot2::rel(0.6)), 
                           strip.text.x = ggplot2::element_text(size = ggplot2::rel(0.8))) +
            ggplot2::scale_color_identity()
        } else {
          gA2 <- gA2 + ggplot2::geom_density(data = peakGeneCorrelations.all[indexCur,],
                                             mapping = ggplot2::aes(peak_gene.p_raw, color = dist_class, linetype = r_pos_class),
                                             show.legend = FALSE) +
            ggplot2::geom_density(data = peakGeneCorrelations.all[indexCur,],
                                  mapping = ggplot2::aes(peak_gene.p_raw, color = dist_class, linetype = r_pos_class))
          
        }
        
        density_plots <- c(density_plots, gA2)
        
        # Helper function to retrieve all tables and data aggregation steps for subsequent visualization
        tbl.l = .createTables_peakGeneQC(peakGeneCorrelations.all[indexCur,], networkType_details, colors_vec, range)
        
        # Store in the GRN object for independent analysis
        GRN@stats$peak_genes[[paste0(geneTypesSelected, collapse = "+")]] = tbl.l
        
        xlabel = paste0("Correlation raw\np-value (binned)")
        
        
        xlabels = levels(tbl.l$d_merged$peak_gene.p.raw.class)
        xlabels[setdiff(seq_len(length(xlabels)), c(1, floor(length(xlabels)/2), length(xlabels)))] <- ""
        
        if (res == 0.1){
          gB3 = ggplot2::ggplot(tbl.l$d_merged, ggplot2::aes(.data$peak_gene.p.raw.class, .data$ratio, fill = .data$classAll)) + 
            ggplot2::geom_bar(stat = "identity", position = "dodge", na.rm = TRUE, width = 0.5) + 
            ggplot2::geom_hline(yintercept = 1, linetype = "dotted") + 
            ggplot2::xlab(xlabel) + ggplot2::ylab("Ratio") +
            ggplot2::scale_fill_manual("Class", values = c(dist_class, r_pos_class), 
                                       labels = c("real", "background", "r+ (r>0)", "r- (r<=0)"), 
            ) + # labels vector can be kind of manually specified here because the levels were previosly defined in a certain order
            ggplot2::scale_x_discrete(labels = xlabels_peakGene_praw.class) +
            ggplot2::theme_bw() +  
            #ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8), strip.background = ggplot2::element_blank(), strip.placement = "outside", axis.title.y = ggplot2::element_blank()) +
            # ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8) , axis.title.y = ggplot2::element_blank()) +
            theme_main +
            ggplot2::facet_wrap(~ factor(set), nrow = 2, scales = "free_y", strip.position = "left") 
          
        } else {
          gB3 = gB3 + ggplot2::geom_bar(tbl.l$d_merged, mapping = ggplot2::aes(.data$peak_gene.p.raw.bins, .data$ratio, fill = .data$classAll), stat = "identity", position = "dodge", na.rm = TRUE, width = 0.5)
        }
        
        ratio_plots <- c(ratio_plots, gB3)
        #  plot two FDR plots as well: (fraction of negative / negative+positive) and (fraction of background / background + real)
        # that could give an indication of whether an FDR based on the background or based on the negative correlations would be more stringent
        
        # R PEAK GENE #
        
        xlabel = paste0("Correlation coefficient r")
        
        sum_real = table(peakGeneCorrelations.all[indexCur,]$class)[names(dist_class)[1]]
        sum_rnd  = table(peakGeneCorrelations.all[indexCur,]$class)[names(dist_class)[2]]
        binData.r = peakGeneCorrelations.all[indexCur,] %>%
          dplyr::group_by(class) %>%
          dplyr::count(.data$peak_gene.r.class) %>%
          dplyr::mutate(nnorm = dplyr::case_when(class == !! (names(dist_class)[1]) ~ .data$n / (sum_real / sum_rnd), 
                                                 TRUE ~ as.double(.data$n)))
        
        binData.r <- binData.r %>% 
          mutate(dist_class = dist_class[class])
        
        xlabel = paste0("Correlation coefficient r (binned)")
        
        
        if (res == 0.1){
          gD = ggplot2::ggplot(binData.r, ggplot2::aes(.data$peak_gene.r.class, .data$nnorm, group = .data$class), fill = binData.r$dist_class)+#.data$class)) + 
            #ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(preserve = "single"), na.rm = FALSE, width = 0.5) +
            ggplot2::geom_line(ggplot2::aes(.data$peak_gene.r.class, .data$nnorm, group = .data$class), color = binData.r$dist_class, stat = "identity") +
            #ggplot2::geom_density(ggplot2::aes(.data$peak_gene.r.class, .data$nnorm, group = .data$class, color = .data$class, alpha = .2), stat = "identity") + 
            #ggplot2::scale_fill_manual("Group", labels = names(dist_class), values = dist_class) +
            #ggplot2::scale_color_manual("Group", labels = names(dist_class), values = dist_class) +
            #ggplot2::scale_x_discrete(labels = xlabels_peakGene_r.class2, drop = FALSE) +
            ggplot2::theme_bw() + ggplot2::theme(legend.position = "none") +
            ggplot2::xlab(xlabel) + ggplot2::ylab("Abundance") +
            theme_main  #+
          #ggplot2::scale_y_continuous(labels = scales::scientific)
        } else{
          gD = gD + #ggplot2::geom_bar(binData.r, mapping = ggplot2::aes(.data$peak_gene.r.class, .data$nnorm, group = .data$class, fill = .data$class), stat = "identity", position = ggplot2::position_dodge(preserve = "single"), na.rm = FALSE, width = 0.5) + 
            ggplot2::geom_line(binData.r, mapping = ggplot2::aes(.data$peak_gene.r.class, .data$nnorm, group = .data$class), color = binData.r$dist_class, stat = "identity")
        }
        
        
        rcorr_plots <- c(rcorr_plots, gD)
        
        mainTitle = paste0("Summary QC (TF: ", TFCur, ", gene type: ", paste0(geneTypesSelected, collapse = "+"), ",\n", .prettyNum(range), " bp promoter range)")
        
        plots_all = ( ((gA2) + 
                         patchwork::plot_layout(widths = c(4))) / ((gD) + 
                                                                     patchwork::plot_layout(widths = c(4))) ) + 
          patchwork::plot_layout(heights = c(2,1), guides = 'collect') +
          patchwork::plot_annotation(title = mainTitle, theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
        
        plot(plots_all)
      }
      pageCounter = pageCounter + 1
      
      
    } # end for each TF
    
  }
  ii <- ii + 1
}




####################
# HELPER FUNCTIONS #
####################

.prettyNum <- function(number, verbose = TRUE) {
  prettyNum(number, big.mark = ",", scientific = FALSE)
}

.customLabeler <- function(tbl_freq) {
  tbl_freq_label = paste0(names(tbl_freq), " (", tbl_freq, ")")
  names(tbl_freq_label) = names(tbl_freq)
  ggplot2::as_labeller(tbl_freq_label)
}


.createTables_peakGeneQC <- function(peakGeneCorrelations.all.cur, networkType_details, colors_vec, range) {
  
  d = peakGeneCorrelations.all.cur %>% 
    dplyr::group_by(.data$r_positive, class, .data$peak_gene.p.raw.class) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::ungroup()
  
  # Some classes might be missing, add them here with explicit zeros
  for (r_pos in c(TRUE, FALSE)) {
    for (classCur in networkType_details) {
      for (pclassCur in levels(peakGeneCorrelations.all.cur$peak_gene.p.raw.class)) {
        
        row = which(d$r_positive == r_pos & d$class == classCur & as.character(d$peak_gene.p.raw.class) == as.character(pclassCur))
        if (length(row) == 0) {
          d = tibble::add_row(d, r_positive = r_pos, class = classCur, peak_gene.p.raw.class = pclassCur, n = 0)
        }
      }
    }
  }
  
  # Restore the "ordered" factor for class
  d$peak_gene.p.raw.class = factor(d$peak_gene.p.raw.class, ordered = TRUE, levels =  levels(peakGeneCorrelations.all.cur$peak_gene.p.raw.class))
  
  
  # Normalization factors
  dsum = d %>%
    dplyr::group_by(.data$r_positive, .data$class) %>%
    dplyr::summarise(sum_n = sum(.data$n))
  
  
  # Summarize per bin
  d2 = d %>%
    dplyr::group_by(class, .data$peak_gene.p.raw.class) %>%
    dplyr::summarise(sum_pos = .data$n[.data$r_positive],
                     sum_neg = .data$n[!.data$r_positive]) %>%
    dplyr::mutate(ratio_pos_raw = .data$sum_pos / .data$sum_neg,
                  fraction_pos = .data$sum_pos / (.data$sum_pos + .data$sum_neg),
                  fraction_neg = 1 - .data$fraction_pos) %>%
    dplyr::ungroup()
  
  # Compare between real and background
  normFactor_real = dplyr::filter(dsum, class ==  !! (networkType_details[1])) %>%  dplyr::pull(.data$sum_n) %>% sum() /
    dplyr::filter(dsum, class ==  !! (networkType_details[2])) %>%  dplyr::pull(.data$sum_n) %>% sum()
  
  # ratio_norm not used currently, no normalization necessary here or not even useful because we dont want to normalize the r_pos and r_neg ratios: These are signal in a way. Only when comparing between real and background, we have to account for sample size for corrections
  d3 = d %>%
    dplyr::group_by(.data$peak_gene.p.raw.class, .data$r_positive) %>%
    dplyr::summarise(n_real     = .data$n[class == !! (names(colors_vec)[1]) ],
                     n_background = .data$n[class == !! (names(colors_vec)[2]) ]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(ratio_real_raw = .data$n_real / .data$n_background,
                  ratio_real_norm = .data$ratio_real_raw / normFactor_real,
                  enrichment_real      = .data$n_real / (.data$n_real + .data$n_background),
                  enrichment_real_norm = (.data$n_real / normFactor_real) / ((.data$n_real / normFactor_real) + .data$n_background)) 
  
  
  stopifnot(identical(levels(d2$peak_gene.p.raw.class), levels(d3$peak_gene.p.raw.class)))
  # 2 enrichment bar plots but combined using facet_wrap
  d2$set = "r+ / r-"; d3$set = "real / background" 
  d_merged <- tibble::tibble(peak_gene.p.raw.class = c(as.character(d2$peak_gene.p.raw.class), 
                                                       as.character(d3$peak_gene.p.raw.class)),
                             ratio = c(d2$ratio_pos_raw, d3$ratio_real_norm),
                             classAll = c(as.character(d2$class), d3$r_positive),
                             set = c(d2$set, d3$set)) %>%
    dplyr::mutate(classAll = factor(.data$classAll, levels = c(paste0("real_",range), paste0("random_",range), "TRUE", "FALSE")),
                  peak_gene.p.raw.class = factor(.data$peak_gene.p.raw.class, levels = levels(d2$peak_gene.p.raw.class)))
  
  d4 = tibble::tibble(peak_gene.p.raw.class = unique(d$peak_gene.p.raw.class), 
                      n_rpos_real = NA_integer_, n_rpos_random = NA_integer_,
                      n_rneg_real = NA_integer_, n_rneg_random = NA_integer_,
                      ratio_background_real_rpos_norm = NA_real_,
                      ratio_background_real_rneg_norm = NA_real_)
  
  for (i in seq_len(nrow(d4))) {
    row_d2 = which(d2$class == networkType_details[1] & d2$peak_gene.p.raw.class == d4$peak_gene.p.raw.class[i])
    stopifnot(length(row_d2) == 1)
    d4[i, "n_rpos_real"] = d2[row_d2, "sum_pos"] %>% unlist()
    d4[i, "n_rneg_real"] = d2[row_d2, "sum_neg"] %>% unlist()
    row_d2 = which(d2$class == paste0("random_",range) & d2$peak_gene.p.raw.class == d4$peak_gene.p.raw.class[i])
    d4[i, "n_rpos_random"] = d2[row_d2, "sum_pos"] %>% unlist()
    d4[i, "n_rneg_random"] = d2[row_d2, "sum_neg"] %>% unlist()
    
    row_d3 = which(d3$r_positive == TRUE & d3$peak_gene.p.raw.class == d4$peak_gene.p.raw.class[i])
    d4[i, "ratio_background_real_rpos_norm"] = 1 - d3[row_d3, "ratio_real_norm"] %>% unlist()
    row_d3 = which(d3$r_positive == FALSE & d3$peak_gene.p.raw.class == d4$peak_gene.p.raw.class[i])
    d4[i, "ratio_background_real_rneg_norm"] = 1 - d3[row_d3, "ratio_real_norm"] %>% unlist()
  }
  
  d4 = d4 %>%
    dplyr::mutate(ratio_rneg_rpos_real = .data$n_rneg_real / (.data$n_rneg_real + .data$n_rpos_real),
                  ratio_rneg_rpos_random = .data$n_rneg_random / (.data$n_rneg_random + .data$n_rpos_random),
                  peak_gene.p.raw.class.bin = as.numeric(.data$peak_gene.p.raw.class)) %>%
    dplyr::arrange(.data$peak_gene.p.raw.class.bin)
  
  d4_melt = reshape2::melt(d4, id  = c("peak_gene.p.raw.class.bin", "peak_gene.p.raw.class")) %>%
    dplyr::filter(grepl("ratio", .data$variable))
  
  
  list(d = d, d2 = d2, d3 = d3, d4 = d4, d4_melt = d4_melt, d_merged = d_merged)
  
}

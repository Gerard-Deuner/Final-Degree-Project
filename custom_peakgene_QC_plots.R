###########################################################
# CUSTOM PEAK-GENE QUALITY CONTROL PLOTS (MULTIPLE EGRNS) #
###########################################################

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

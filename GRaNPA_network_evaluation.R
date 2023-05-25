#############################
# GRaNPA NETWORK EVALUATION #
#############################

# load libraries
library(GRaNIE)
library(GRaNPA)
library(dplyr)
library(qs)
library(readxl)
library(data.table)
library(ggpubr)

# define resolution 
res <- 6 # 8, 6, 10, 0.25, 18

# load eGRN 
GRN = fread(paste0("/g/scb/zaugg/deuner/GRaNIE/outputdata/batch_mode/combined_batch_mode_spearman_nomicro/Batch_Mode_Outputs/output_pseudobulk_clusterRes", res, "_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/connections_TFPeak0.2_peakGene0.1.tsv.gz"))
head(GRN)


##########################################
# Day 4 vs Day 2 Differential Expression #
##########################################

# load DE dataset
# (the DE dataset should always contain these three columns: ‘ENSEMBL’, ‘padj’ and ‘logFC’)
DE <- fread("/g/scb/zaugg/deuner/GRaNPA/inputdata/DE_timecourse_Day4vsDay0.tsv")
head(DE)

# Run GRaNPA
GRaNPA_d4d0 = GRaNPA::GRaNPA_main_function(DE_data = DE, 
                                             GRN_matrix_filtered = GRN,
                                             DE_pvalue_th = 0.05,
                                             logFC_th = 2,
                                             num_run = 20, #3
                                             num_run_CR = 20, #2
                                             num_run_random = 20, #3
                                             cores = 20,
                                             importance = "permutation",
                                             ML_type = "regression",
                                             control = "cv",
                                             train_part = 1)

# Check results
p1 <- plot_GRaNPA_density_mod(GRaNPA.object = GRaNPA_d4d0, 
                                  plot_name = "density.pdf", 
                                  outputFolder = ".", width = 4, height = 4) + 
  labs(x = "R²[10-fold cross validation]", y = "Density")

# See the relation between actual log fold change and predicted log fold change
p2 <- plot_GRaNPA_scatter_mod(GRaNPA.object = GRaNPA_d4d0, 
                                  plot_name = "scatter.pdf", 
                                  output = ".", width = 4, height = 4) + 
  labs(x = "Predicted Log Fold Change", y = "True Log Fold Change")

# Investigate important TFs
p3 <- GRaNPA::plot_GRaNPA_TF_imp(GRaNPA.object = GRaNPA_d4d0, plot_name = "TF_imp.pdf", output = ".", width = 4, height = 4) 

pl <- ggarrange(p1, p2, p3, labels = c("A", "B", "C")) %>% 
  annotate_figure(top = text_grob("GRaNPA Analysis - Day 4 vs. Day 0", 
                                        color = "black", face = "bold", size = 14))
pl
ggsave(paste0("/g/scb/zaugg/deuner/GRaNPA/figures/res", res, "Day4vsDay0.png"), pl, device = "png")

##########################################
# Day 2 vs Day 0 Differential Expression #
##########################################

# load DE dataset
# (the DE dataset should always contain these three columns: ‘ENSEMBL’, ‘padj’ and ‘logFC’)
DE <- fread("/g/scb/zaugg/deuner/GRaNPA/inputdata/DE_timecourse_Day2vsDay0.tsv")
head(DE)

# Run GRaNPA
GRaNPA_d2d0 = GRaNPA::GRaNPA_main_function(DE_data = DE, 
                                           GRN_matrix_filtered = GRN,
                                           DE_pvalue_th = 0.05,
                                           logFC_th = 2,
                                           num_run = 20,
                                           num_run_CR = 20,
                                           num_run_random = 20,
                                           cores = 20,
                                           importance = "permutation",
                                           ML_type = "regression",
                                           control = "cv",
                                           train_part = 1)

# Check results
p1 <- plot_GRaNPA_density_mod(GRaNPA.object = GRaNPA_d2d0, 
                              plot_name = "density.pdf", 
                              outputFolder = ".", width = 4, height = 4) + 
  labs(x = "R²[10-fold cross validation]", y = "Density")

# See the relation between actual log fold change and predicted log fold change
p2 <- plot_GRaNPA_scatter_mod(GRaNPA.object = GRaNPA_d2d0, 
                              plot_name = "scatter.pdf", 
                              output = ".", width = 4, height = 4) + 
  labs(x = "Predicted Log Fold Change", y = "True Log Fold Change")

# Investigate important TFs
p3 <- GRaNPA::plot_GRaNPA_TF_imp(GRaNPA.object = GRaNPA_d2d0, plot_name = "TF_imp.pdf", output = ".", width = 4, height = 4) 

pl <- ggarrange(p1, p2, p3, labels = c("A", "B", "C")) %>% 
  annotate_figure(top = text_grob("GRaNPA Analysis - Day 2 vs. Day 0", 
                                  color = "black", face = "bold", size = 14))
pl
ggsave(paste0("/g/scb/zaugg/deuner/GRaNPA/figures/res", res, "Day2vsDay0.png"), pl, device = "png")

##########################################
# Day 4 vs Day 2 Differential Expression #
##########################################

# load DE dataset
# (the DE dataset should always contain these three columns: ‘ENSEMBL’, ‘padj’ and ‘logFC’)
DE <- fread("/g/scb/zaugg/deuner/GRaNPA/inputdata/DE_timecourse_Day4vsDay2.tsv")
head(DE)

# Run GRaNPA
GRaNPA_d4d2 = GRaNPA::GRaNPA_main_function(DE_data = DE, 
                                           GRN_matrix_filtered = GRN,
                                           DE_pvalue_th = 0.05,
                                           logFC_th = 2,
                                           num_run = 20,
                                           num_run_CR = 20,
                                           num_run_random = 20,
                                           cores = 20,
                                           importance = "permutation",
                                           ML_type = "regression",
                                           control = "cv",
                                           train_part = 1)

# Check results
p1 <- plot_GRaNPA_density_mod(GRaNPA.object = GRaNPA_d4d2, 
                              plot_name = "density.pdf", 
                              outputFolder = ".", width = 4, height = 4) + 
  labs(x = "R²[10-fold cross validation]", y = "Density")

# See the relation between actual log fold change and predicted log fold change
p2 <- plot_GRaNPA_scatter_mod(GRaNPA.object = GRaNPA_d4d2, 
                              plot_name = "scatter.pdf", 
                              output = ".", width = 4, height = 4) + 
  labs(x = "Predicted Log Fold Change", y = "True Log Fold Change")

# Investigate important TFs
p3 <- GRaNPA::plot_GRaNPA_TF_imp(GRaNPA.object = GRaNPA_d4d2, plot_name = "TF_imp.pdf", output = ".", width = 4, height = 4) 

pl <- ggarrange(p1, p2, p3, labels = c("A", "B", "C")) %>% 
  annotate_figure(top = text_grob("GRaNPA Analysis - Day 4 vs. Day 2", 
                                  color = "black", face = "bold", size = 14))
pl
ggsave(paste0("/g/scb/zaugg/deuner/GRaNPA/figures/res", res, "Day4vsDay2.png"), pl, device = "png")

####################
# eGRNs Comparison #
####################

# Plot R square
cpl <- GRaNPA::plot_GRaNPA_boxplot(GRaNPA.object_list = list(GRaNPA_d4d0, GRaNPA_d2d0, GRaNPA_d4d2),  
                            name_list = c("Day4 vs Day0", "Day2 vs Day0", "Day4 vs Day2"), 
                            plot_name = "combined_spearman_res8_GRaNPA.pdf", 
                            output = "." , 
                            width = 8,
                            height = 4) + 
  geom_hline(yintercept = 0.05, col = "red", linetype = "dashed")

ggsave(paste0("/g/scb/zaugg/deuner/GRaNPA/figures/res", res, "eGRNs_comparison"), cpl, device = "png")


#########################
# MODIFY SOME FUNCTIONS #
#########################

plot_GRaNPA_scatter_mod <- function(
    GRaNPA.object, ## this should be a object from prediction power pipeline.
    outputFolder,
    plot_name,
    width = 5,
    height = 5
) {
  
  
  data = GRaNPA.object$normal_data
  pred = GRaNPA.object$normal_models[[1]]$finalModel$predictions
  df = data.frame(
    actual = data$y,
    predicted =  pred
  )
  
  
  range.min = min(df)
  range.max = max(df)
  
  
  fig_scat = df %>% dplyr::mutate( col = ifelse(sign(actual*predicted)==1 , "A","B")) %>%
    ggplot2::ggplot( aes ( x = predicted , y = actual )) +
    ggplot2::annotate(geom = "rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
                      fill = "#EA302E", colour = "darkgray", alpha = 0.25) +
    ggplot2::annotate(geom = "rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
                      fill = "#98D3F1", colour = "darkgray", alpha = 0.25) +
    ggplot2::annotate(geom = "rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
                      fill = "#98D3F1", colour = "darkgray", alpha = 0.25) +
    ggplot2::annotate(geom = "rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
                      fill = "#EA302E", colour = "darkgray", alpha = 0.25) +
    ggplot2::geom_point(aes(color = col)) +
    ggplot2::geom_smooth(method = "lm" , se = FALSE, color = "dodgerblue") +
    ggplot2::coord_fixed() + ggplot2::theme_bw() +
    ggplot2::ylim(range.min,range.max) + ggplot2::guides(color = FALSE) + ggplot2::xlim(range.min,range.max) + 
    ggplot2::scale_color_manual(values = c("black","#EA302E")) 
    ggplot2::theme(axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 15, color = "gray25"))
    
    
  
  ggplot2::ggsave(filename = plot_name , plot = fig_scat , device = "pdf" , path = outputFolder,
                  width = width , height = height  )
  return(fig_scat)
}

plot_GRaNPA_density_mod <- function(
    GRaNPA.object, ## this should be a object from prediction power pipeline.
    outputFolder,
    plot_name,
    width = 5,
    height = 5,
    range = NULL
) {
  
  data = data.frame(
    r2 = c(GRaNPA.object$normal_dist$N_ranger_model_rsq, GRaNPA.object$random_dist$RG_ranger_model_rsq, GRaNPA.object$CR_dist$RC_ranger_model_rsq),
    type = rep(c("Actual GRN" , "Randomized GRN" , "QC GRN") , c(nrow(GRaNPA.object$normal_dist), nrow(GRaNPA.object$random_dist), nrow(GRaNPA.object$CR_dist)))
  )
  
  grn.labs <- c("Randomized GRN",  "QC GRN" , "Actual GRN")
  names(grn.labs) <- c("Randomized GRN", "QC GRN", "Actual GRN")
  
  data$type = factor(data$type , grn.labs)
  
  p_fig = data %>% ggplot2::ggplot(aes(r2, fill = type)) + ggplot2::geom_density( alpha = 0.7) +
    ggplot2::theme_bw() + 
    ggplot2::scale_fill_manual( values = wesanderson::wes_palette("Darjeeling1", 3, type = "discrete") , labels = grn.labs , name = "GRN type") +
    ggplot2::xlab(expression(paste("R"^"2"))) + ggplot2::ylab("Density") + ggplot2::ggtitle("")# + ggplot2::guides(fill = FALSE)
    
  if (!is.null(range)){
    p_fig = p_fig + ggplot2::xlim(range)
  }
  
  ggplot2::ggsave(plot = p_fig, filename = plot_name, device = "pdf", path = outputFolder, width = width, height = height)
  return(p_fig)
}




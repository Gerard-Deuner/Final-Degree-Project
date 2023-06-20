########################
# PEARSON VS. SPEARMAN #
########################

library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(qs)
library(GRaNIE)

# Cell Type Pseudobulking - WNN - res: 0.5 - 16 clusters

# PEARSON
# load eGRN
GRN.p <- qread("/g/scb/zaugg/deuner/GRaNIE/outputdata/timecourse_pearson_celltype_wnn/output_pseudobulk_celltype_wnn_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/GRN.qs")
GRN.p

# build eGRN
GRN.p = build_eGRN_graph(GRN.p, forceRerun = TRUE)
GRN.p = visualizeGRN(GRN.p, plotAsPDF = FALSE)

# SPEARMAN
GRN.s <- qread("/g/scb/zaugg/deuner/GRaNIE/outputdata/timecourse_spearman_celltype_wnn/output_pseudobulk_celltype_wnn_RNA_limma_quantile_ATAC_DESeq2_sizeFactors/GRN.qs")
GRN.s

# build eGRN
GRN.s = build_eGRN_graph(GRN.s, forceRerun = TRUE)
GRN.s = visualizeGRN(GRN.s, plotAsPDF = FALSE)

GRN.p@connections$all.filtered #more connections
GRN.s@connections$all.filtered 

# Set up the main directory
path = "/g/scb/zaugg/deuner/GRaNIE"
# Load the seurat object
seuratFile = "timecourse.pp.seuratObject.qs"
timecourse.s <- qread(paste0(path,"/tmp/",seuratFile))

# ## PEARSON
# # peak-gene connections here
# GRN.p@connections$all.filtered[["0"]]
# connections.p <- as.data.frame(GRN.p@connections$all.filtered[["0"]])
# 
# # SCT counts here
# GRN.p@data[["RNA"]][["counts"]]
# counts.p <- as.data.frame(GRN.p@data[["RNA"]][["counts"]])
# counts.p$gene.ENSEMBL <- rownames(counts.p)
# counts.p
# 
# # norm ATAC counts here
# GRN.p@data[["peaks"]][["counts"]]
# peaks.p <- as.data.frame(GRN.p@data[["peaks"]][["counts"]])
# peaks.p$peak.ID <- rownames(peaks.p)
# peaks.p
# 
# # connections + counts
# joined.p <- connections.p %>% inner_join(counts.p, by = "gene.ENSEMBL") %>% inner_join(peaks.p, by = "peak.ID")
# 
# 
# ## SPEARMAN
# # peak-gene connections here
# GRN.s@connections$all.filtered[["0"]]
# connections.s <- as.data.frame(GRN.s@connections$all.filtered[["0"]])

order <- c("diff___tiny", "hiPSC", "hiPSC___start.diff", "diff___hiPSC_like", "diff___NPC_like", "diff", "diff___immature.neuron", "neuron___excitatory", 
          "neuron___development_1", "neuron___development_2", "mature.neuron___adhesion", "neuron___mature.neuron", "mature.neuron___small_1",
          "mature.neuron___small_2", "microglia_1", "microglia_2")

# SCT counts here
GRN.s@data[["RNA"]][["counts"]]
counts.s <- as.data.frame(GRN.s@data[["RNA"]][["counts"]])
counts.s$gene.ENSEMBL <- rownames(counts.s)
counts.s <- counts.s[,c(order, "gene.ENSEMBL")]

# norm ATAC counts here
GRN.s@data[["peaks"]][["counts"]]
peaks.s <- as.data.frame(GRN.s@data[["peaks"]][["counts"]])
peaks.s$peak.ID <- rownames(peaks.s)
peaks.s <- peaks.s[,c(order, "peak.ID")]

# # connections + counts
# joined.s <- connections.s %>% inner_join(counts.s, by = "gene.ENSEMBL") %>% inner_join(peaks.s, by = "peak.ID")

# # check shared TF-peak-links
# tf.peak.gene.p <- GRN.p@connections$all.filtered[["0"]]
# tf.peak.gene.s <- GRN.s@connections$all.filtered[["0"]]
# colnames(tf.peak.gene.p)
# 
# tf.shared.p.s <- inner_join(tf.peak.gene.p, tf.peak.gene.s, by = c("TF.ID")) 
# gene.shared.p.s <- inner_join(tf.peak.gene.p, tf.peak.gene.s, by = c("gene.name"))
# peak.shared.p.s <- inner_join(tf.peak.gene.p, tf.peak.gene.s, by = c("peak.ID"))
# 
# View(tf.peak.gene.p)
# tf.peak.gene.s
# 
# tf.shared.p.s
# gene.shared.p.s
# peak.shared.p.s
# 
# 
# # CHECK IF ASSUMPTIONS ARE HELD!
# # Pearson:
# 
# # select top5 most significant peak-gene links
# top5.p <- joined.p %>% top_n(5, peak_gene.r)
# top5.p.1 <- top5.p[,1:32]
# top5.p.2 <- top5.p[,c(1:16, 33:48)]
# # long format
# top5.p.long.1 <- reshape2::melt(top5.p.1 , id.vars = colnames(top5.p)[1:16], variable.name = "cluster.RNA", value.name = "value.RNA")
# top5.p.long.2 <- reshape2::melt(top5.p.2, id.vars = colnames(top5.p)[1:16], variable.name = "cluster.ATAC", value.name = "value.ATAC")
# top5.p.both <- cbind(top5.p.long.1, top5.p.long.2[,17:ncol(top5.p.long.2)])
# 
# # assign informative colors to the cell type clusters
# cols <- c("#FF6600", "#FF0000", "#CC0033", "#99FFFF", "#0099FF", "#33FF33", "#FF00FF", "#990000", "#FF3300", "#FFCC66", "#660066", "#0000FF", "#FFFF33", "#FF6666", "#FF9933", "#66FF99")
# 
# ggplot(top5.p.both, aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
#   geom_point(alpha = 0.75, size = 2) + 
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Peak-gene link with highest correlation coefficient") + 
#   theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   geom_label_repel(data=top5.p.both[top5.p.both$TF.name == "YY1",], aes(label=as.factor(TF.name), alpha=0.7), size=5, force=1.3, max.overlaps = 30) + 
#   scale_color_manual(values = cols)
# 
# table(joined.p$TF.ID)
# 
# 
# ggplot(top5.p.both, aes(x = value.RNA, y = value.ATAC, col = cluster.RNA, label=TF.name)) + 
#   geom_point(alpha = 0.75, size = 2) + 
#   labs(x = "Gene Expression (SCT) log 10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 10 most correlated peak-genes") + 
#   theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   geom_text(hjust=0, vjust=0)
# 
# # set new ids
# new.ids <- as.character(c(1:5))
# top5.p.both$id <- new.ids
# names(new.ids) <- levels(top5.p.both$id)
# top5.p.both
# 
# ggplot(top5.p.both, aes(x = log10(value.RNA), y = value.ATAC, col = id)) + 
#   geom_point(alpha = 0.75, size = 2) + 
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Id", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-genes") + 
#   theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   geom_label_repel(data=top5.p.both[1:5,], aes(label=as.factor(round(peak_gene.r, 4))), size=5)
# 
# 
# 
# p1 <- ggplot(top5.p.both %>% dplyr::filter(id == "1"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
#   geom_point(alpha = 0.75, size = 2) + 
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "TTop 1 most correlated peak-gene link") + 
#   theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   scale_color_manual(values = cols)
# p2 <-ggplot(top5.p.both %>% dplyr::filter(id == "2"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
#   geom_point(alpha = 0.75, size = 2) + 
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 2 most correlated peak-gene link") + 
#   theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   scale_color_manual(values = cols)
# p3 <-  ggplot(top5.p.both %>% dplyr::filter(id == "3"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
#   geom_point(alpha = 0.75, size = 2) + 
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 3 most correlated peak-gene link") + 
#   theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   scale_color_manual(values = cols)
# p4 <-ggplot(top5.p.both %>% dplyr::filter(id == "4"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
#   geom_point(alpha = 0.75, size = 2) + 
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 4 most correlated peak-gene link") + 
#   theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   scale_color_manual(values = cols)
# p5 <-ggplot(top5.p.both %>% dplyr::filter(id == "5"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
#   geom_point(alpha = 0.75, size = 2) + 
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-gene link") + 
#   theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   scale_color_manual(values = cols)
# 
# ggarrange(p1,p2,p3,p4,p5)
# p1
# p2
# 
# 
# 
# # Spearman:
# 
# # select top5 most significant peak-gene links
# ## add also the pearson coefficient as a columns
# top5.s <- joined.s %>% top_n(5, peak_gene.r)
# top5.s.1 <- top5.s[,1:32]
# top5.s.2 <- top5.s[,c(1:16, 33:48)]
# # long format
# top5.s.long.1 <- reshape2::melt(top5.s.1 , id.vars = colnames(top5.s)[1:16], variable.name = "cluster.RNA", value.name = "value.RNA")
# top5.s.long.2 <- reshape2::melt(top5.s.2, id.vars = colnames(top5.s)[1:16], variable.name = "cluster.ATAC", value.name = "value.ATAC")
# top5.s.both <- cbind(top5.s.long.1, top5.s.long.2[,17:ncol(top5.s.long.2)])
# 
# # set new ids
# new.ids <- as.character(c(1:5))
# top5.s.both$id <- new.ids
# names(new.ids) <- levels(top5.s.both$id)
# new.ids
# 
# ggplot(top5.s.both, aes(x = log10(value.RNA), y = value.ATAC, col = id)) + 
#   geom_point(alpha = 0.75, size = 2) + 
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Id", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-genes") + 
#   theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   geom_label_repel(data=top5.s.both[1:5,], aes(label=as.factor(round(peak_gene.r, 4))), size=5)
# 
# 
# ggplot(top5.s.both, aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
#   geom_point(alpha = 0.75, size = 2) + 
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Peak-gene link with highest correlation coefficient") + 
#   theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   scale_color_manual(values = cols)
# 
# p1 <- ggplot(top5.s.both %>% dplyr::filter(id == "1"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
#   geom_point(alpha = 0.75, size = 2) + 
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 1 most correlated peak-gene link") + 
#   theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   scale_color_manual(values = cols)
# p2 <-ggplot(top5.s.both %>% dplyr::filter(id == "2"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
#   geom_point(alpha = 0.75, size = 2) + 
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 2 most correlated peak-gene link") + 
#   theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   scale_color_manual(values = cols)
# p3 <-  ggplot(top5.s.both %>% dplyr::filter(id == "3"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
#   geom_point(alpha = 0.75, size = 2) + 
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 3 most correlated peak-gene link") + 
#   theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   scale_color_manual(values = cols)
# p4 <-ggplot(top5.s.both %>% dplyr::filter(id == "4"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
#   geom_point(alpha = 0.75, size = 2) + 
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 4 most correlated peak-gene link") + 
#   theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   scale_color_manual(values = cols)
# p5 <-ggplot(top5.s.both %>% dplyr::filter(id == "5"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
#   geom_point(alpha = 0.75, size = 2) + 
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-gene link") + 
#   theme_classic() + 
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   scale_color_manual(values = cols)
# 
# ggarrange(p1,p2,p3,p4,p5)
# p1
# p4
# p5

#########################
# ALL PEAK - GENE LINKS #
#########################

# redoit but with all peak-gene links and putting pearson and spearman correlation coefficients in the same dataframe
# extract peak-genes from GRN ran with Spearman correlation
#pg.links <- GRN.s@connections$peak_genes$`0` # IMPORTANT LINE
col.names <- c("peak.ID", "gene.ENSEMBL", "peak_gene.distance", "peak_gene.r", "peak_gene.p_raw") 
pg.links <- GRN.s@connections$all.filtered$`0`
pg.links <- pg.links[col.names]
colnames(pg.links)[4] <- "spearman.r"
# add pearson correlation and rawp values
pg.links$pearson.r <- rep(0, nrow(pg.links))#GRN.p@connections$peak_genes$`0`$peak_gene.r
pg.links$pearson.p_raw <- rep(0, nrow(pg.links))#GRN.p@connections$peak_genes$`0`$peak_gene.p_raw

# add id column
pg.links$id <- as.character(c(1:nrow(pg.links)))

# add gene expresion and peak accessibility counts for each connection
pg.joined.s <- pg.links %>% inner_join(counts.s, by = "gene.ENSEMBL") %>% inner_join(peaks.s, by = "peak.ID")
pg.joined.s$id <- NULL

############
# SPEARMAN #
############

# select top5 most significant peak-gene links
## add also the pearson coefficient as a columns
pg.top5.s <- pg.joined.s %>% top_n(5, spearman.r)
pg.top5.s.1 <- pg.top5.s[,1:23]
pg.top5.s.2 <- pg.top5.s[,c(1:7, 24:39)]

# long format
pg.top5.s.long.1 <- reshape2::melt(pg.top5.s.1 , id.vars = colnames(pg.top5.s)[1:7], variable.name = "cluster.RNA", value.name = "value.RNA")
pg.top5.s.long.2 <- reshape2::melt(pg.top5.s.2, id.vars = colnames(pg.top5.s)[1:7], variable.name = "cluster.ATAC", value.name = "value.ATAC")
pg.top5.s.both <- cbind(pg.top5.s.long.1, pg.top5.s.long.2[,8:ncol(pg.top5.s.long.2)])
pg.top5.s.both

# set new ids
new.ids <- rep(as.character(c(1:5)), nrow(pg.top5.s.both)/5)
pg.top5.s.both$ids <- new.ids


corMethod <- "spearman"

cols <- c("gray", "#FFFF00", "#FFCC33", "#FF9966", "#FF6633", "red", "#990033", "#996699", "#9933FF", "#6600CC", "#330099", "#3333FF","#0000FF", "#3366FF" , "#00C000", "green")
vals <- c("Diff - small", "hiPSC", "Start Diff", "Diff - hiPSC-like", "Diff - NPC-like", "Diff", "Diff - Immature Neuron", "Excitatory Neuron", 
          "Neuron Development 1", "Neuron Development 2", "Mature Neuron - adhesion", "Mature Neuron", "Mature Neuron - small 1",
          "Mature Neuron - small 2", "Microglia 1", "Microglia 2")
digits_round <- 3
pg.titleCur = paste0("Peak=", pg.top5.s.both %>% dplyr::filter(ids == "1") %>% pull(peak.ID) %>% dplyr::first(), ", Gene=", pg.top5.s.both %>% dplyr::filter(ids == "1") %>% pull(gene.ENSEMBL) %>% dplyr::first(), 
                     ":\nCor = ", round(pg.top5.s.both %>% dplyr::filter(ids == "1") %>% pull(spearman.r) %>% dplyr::first(), digits_round), " (", corMethod, "), raw p-value = ", round(pg.top5.s.both %>% dplyr::filter(ids == "1") %>% pull(peak_gene.p_raw) %>% dplyr::first(), digits_round),
                     "\n", "# eGRN link: highest correlation")

x_min = min(-0.01, min(pg.top5.s.both$value.RNA, na.rm = TRUE))
y_min = min(-0.01, min(pg.top5.s.both$value.ATAC, na.rm = TRUE))

pgs1 <- ggplot(pg.top5.s.both %>% dplyr::filter(ids == "1"), aes(x = value.RNA, y = value.ATAC, col = cluster.RNA)) + 
  geom_smooth(method = "lm", formula = "y ~ x", col = "black", na.rm = TRUE) + 
  #geom_point(size = rel(0.5)) + 
  geom_point(alpha = 0.75, size = 2) + 
  xlim(x_min, NA) + ylim(y_min, NA) +
  labs(x = "(Normalized) gene expression", y = "(Normalized) peak accessibility", col = "Cluster / Cell Type", title = pg.titleCur) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 10)) + 
  scale_color_manual(values = cols,
                     labels = vals)

ggsave("/g/scb/zaugg/deuner/figs/spearman_corr_plot_peakgene_top1link.png", pgs1, device = "png")

# ggplot(pg.top5.s.both, aes(x = log10(value.RNA), y = value.ATAC, col = id)) +
#   geom_point(alpha = 0.75, size = 2) +
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Id", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-genes") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# p1 <- ggplot(pg.top5.s.both %>% dplyr::filter(id == "316796"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) +
#   geom_point(alpha = 0.75, size = 2) +
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 1 most correlated peak-gene link. Id: 316796") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = cols)
# p2 <-ggplot(pg.top5.s.both %>% dplyr::filter(id == "380465"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) +
#   geom_point(alpha = 0.75, size = 2) +
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 2 most correlated peak-gene link. Id: 380465") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = cols)
# p3 <-  ggplot(pg.top5.s.both %>% dplyr::filter(id == "484736"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) +
#   geom_point(alpha = 0.75, size = 2) +
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 3 most correlated peak-gene link. Id: 484736") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = cols)
# p4 <-ggplot(pg.top5.s.both %>% dplyr::filter(id == "714900"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) +
#   geom_point(alpha = 0.75, size = 2) +
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 4 most correlated peak-gene link. Id: 714900") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = cols)
# p5 <-ggplot(pg.top5.s.both %>% dplyr::filter(id == "714906"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) +
#   geom_point(alpha = 0.75, size = 2) +
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-gene link. Id: 714906") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = cols)
# 
# ggarrange(p1,p2,p3,p4,p5)
# p1

###########
# PEARSON #
###########

# redoit but with all peak-gene links and putting pearson and spearman correlation coefficients in the same dataframe
# extract peak-genes from GRN ran with Spearman correlation
#pg.links <- GRN.s@connections$peak_genes$`0` # IMPORTANT LINE
col.names <- c("peak.ID", "gene.ENSEMBL", "peak_gene.distance", "peak_gene.r", "peak_gene.p_raw") 
pg.links <- GRN.p@connections$all.filtered$`0`
pg.links <- pg.links[col.names]
colnames(pg.links)[4] <- "pearson.r"
# add pearson correlation and rawp values
pg.links$spearman.r <- rep(0, nrow(pg.links)) #GRN.p@connections$peak_genes$`0`$peak_gene.r
pg.links$spearman.p_raw <- rep(0, nrow(pg.links)) #GRN.p@connections$peak_genes$`0`$peak_gene.p_raw

# add id column
pg.links$id <- as.character(c(1:nrow(pg.links)))

# add gene expresion and peak accessibility counts for each connection
pg.joined.s <- pg.links %>% inner_join(counts.s, by = "gene.ENSEMBL") %>% inner_join(peaks.s, by = "peak.ID")
pg.joined.s$id <- NULL

# select top5 most significant peak-gene links
## add also the pearson coefficient as a columns
pg.top5.s <- pg.joined.s %>% top_n(5, pearson.r)
pg.top5.s.1 <- pg.top5.s[,1:23]
pg.top5.s.2 <- pg.top5.s[,c(1:7, 24:39)]
# long format
pg.top5.s.long.1 <- reshape2::melt(pg.top5.s.1 , id.vars = colnames(pg.top5.s)[1:7], variable.name = "cluster.RNA", value.name = "value.RNA")
pg.top5.s.long.2 <- reshape2::melt(pg.top5.s.2, id.vars = colnames(pg.top5.s)[1:7], variable.name = "cluster.ATAC", value.name = "value.ATAC")
pg.top5.s.both <- cbind(pg.top5.s.long.1, pg.top5.s.long.2[,8:ncol(pg.top5.s.long.2)])
pg.top5.s.both

# set new ids
new.ids <- rep(as.character(c(1:5)), nrow(pg.top5.s.both)/5)
pg.top5.s.both$idp <- new.ids

# peak-gene links #

# correlation: pearson

corMethod <- "pearson"

digits_round <- 3
pg.titleCur = paste0("Peak=", pg.top5.s.both %>% dplyr::filter(idp == "1") %>% pull(peak.ID) %>% dplyr::first(), ", Gene=", pg.top5.s.both %>% dplyr::filter(idp == "1") %>% pull(gene.ENSEMBL) %>% dplyr::first(), 
                     ":\nCor = ", round(pg.top5.s.both %>% dplyr::filter(idp == "1") %>% pull(pearson.r) %>% dplyr::first(), digits_round), " (", corMethod, "), raw p-value = ", round(pg.top5.s.both %>% dplyr::filter(idp == "1") %>% pull(peak_gene.p_raw) %>% dplyr::first(), digits_round),
                     "\n", "# eGRN link: highest correlation")

x_min = min(-0.01, min(pg.top5.s.both$value.RNA, na.rm = TRUE))
y_min = min(-0.01, min(pg.top5.s.both$value.ATAC, na.rm = TRUE))

pgp1 <- ggplot(pg.top5.s.both %>% dplyr::filter(idp == "1"), aes(x = value.RNA, y = value.ATAC, col = cluster.RNA)) + 
  geom_smooth(method = "lm", formula = "y ~ x", col = "black", na.rm = TRUE) + 
  #geom_point(size = rel(0.5)) + 
  geom_point(alpha = 0.75, size = 2) + 
  xlim(x_min, NA) + ylim(y_min, NA) +
  labs(x = "", y = "(Normalized) peak accessibility", col = "Cluster / Cell Type", title = pg.titleCur) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 10)) + 
  scale_color_manual(values = cols,
                     labels = vals)

ggsave("/g/scb/zaugg/deuner/figs/pearson_corr_plot_peakgene_top1link.png", pgp1, device = "png")

# ggplot(pg.top5.s.both, aes(x = log10(value.RNA), y = value.ATAC, col = id)) +
#   geom_point(alpha = 0.75, size = 2) +
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Id", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-genes") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# p1 <- ggplot(pg.top5.s.both %>% dplyr::filter(id == "146686"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) +
#   geom_point(alpha = 0.75, size = 2) +
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "TTop 1 most correlated peak-gene link. Id: 146686") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = cols)
# p2 <-ggplot(pg.top5.s.both %>% dplyr::filter(id == "18914"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) +
#   geom_point(alpha = 0.75, size = 2) +
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 2 most correlated peak-gene link. Id: 18914") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = cols)
# p3 <-  ggplot(pg.top5.s.both %>% dplyr::filter(id == "441623"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) +
#   geom_point(alpha = 0.75, size = 2) +
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 3 most correlated peak-gene link. Id: 441623") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = cols)
# p4 <-ggplot(pg.top5.s.both %>% dplyr::filter(id == "634342"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) +
#   geom_point(alpha = 0.75, size = 2) +
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 4 most correlated peak-gene link. Id: 634342") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = cols)
# p5 <-ggplot(pg.top5.s.both %>% dplyr::filter(id == "757121"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) +
#   geom_point(alpha = 0.75, size = 2) +
#   labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-gene link. Id: 757121") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = cols)
# 
# ggarrange(p1,p2,p3,p4,p5)
# p1

#################
# TF-PEAK LINKS #
#################

############
# SPEARMAN #
############

# redoit but with all TF-peak links and putting pearson and spearman correlation coefficients in the same dataframe
# extract TF-peak from GRN ran with Spearman correlation
#tfp.links <- GRN.s@connections$TF_peaks$`0`$main # IMPORTANT LINE
col.names <- c("TF.ID",                  "TF_peak.r_bin",          "TF_peak.r",             
               "TF_peak.fdr",            "TF.name",       "peak.ID",               
               "TF_peak.fdr_direction",  "TF_peak.connectionType")
tfp.links <- GRN.s@connections$all.filtered$`0`
tfp.links <- tfp.links[,col.names]
colnames(tfp.links)[3] <- "spearman.r"
TF.ENSEMBL = GRN.s@annotation$TFs %>% dplyr::select(TF.ID, TF.ENSEMBL)
tfp.links <- inner_join(tfp.links, TF.ENSEMBL, by = "TF.ID")
names(tfp.links)[9] <- "gene.ENSEMBL"

# add gene expresion and peak accessibility counts for each connection
tfp.joined.s <- tfp.links %>% inner_join(counts.s, by = "gene.ENSEMBL") %>% inner_join(peaks.s, by = "peak.ID")

# select top5 most significant peak-gene links
## add also the pearson coefficient as a columns
tfp.top5.s <- tfp.joined.s %>% top_n(5, spearman.r)
tfp.top5.s <- tfp.top5.s[1:5,]
tfp.top5.s.1 <- tfp.top5.s[,1:25]
tfp.top5.s.2 <- tfp.top5.s[,c(1:9, 26:41)]
# long format
tfp.top5.s.long.1 <- reshape2::melt(tfp.top5.s.1 , id.vars = colnames(tfp.top5.s)[1:9], variable.name = "cluster.RNA", value.name = "value.RNA")
tfp.top5.s.long.2 <- reshape2::melt(tfp.top5.s.2, id.vars = colnames(tfp.top5.s)[1:9], variable.name = "cluster.ATAC", value.name = "value.ATAC")
tfp.top5.s.both <- cbind(tfp.top5.s.long.1, tfp.top5.s.long.2[,10:ncol(tfp.top5.s.long.2)])
tfp.top5.s.both %>% names

# set new ids
new.ids <- rep(as.character(c(1:5)), nrow(tfp.top5.s.both)/5)
tfp.top5.s.both$ids <- new.ids

tfp.top5.s.both$TF.ID


corMethod <- "spearman"

digits_round <- 3
tfp.titleCur = paste0("Peak=", tfp.top5.s.both %>% dplyr::filter(ids == "1") %>% pull(peak.ID) %>% dplyr::first(), ", TF=", tfp.top5.s.both %>% dplyr::filter(ids == "1") %>% pull(TF.ID) %>% dplyr::first(), 
                     ":\nCor = ", round(tfp.top5.s.both %>% dplyr::filter(ids == "1") %>% pull(spearman.r) %>% dplyr::first(), digits_round), " (", corMethod, "), raw p-value = ", round(tfp.top5.s.both %>% dplyr::filter(ids == "1") %>% pull(TF_peak.fdr) %>% dplyr::first(), digits_round),
                     "\n", "# eGRN link: highest correlation")

x_min = min(-0.01, min(tfp.top5.s.both$value.RNA, na.rm = TRUE))
y_min = min(-0.01, min(tfp.top5.s.both$value.ATAC, na.rm = TRUE))

tfps1 <- ggplot(tfp.top5.s.both %>% dplyr::filter(ids == "1"), aes(x = value.RNA, y = value.ATAC, col = cluster.RNA)) + 
  geom_smooth(method = "lm", formula = "y ~ x", col = "black", na.rm = TRUE) + 
  #geom_point(size = rel(0.5)) + 
  geom_point(alpha = 0.75, size = 2) + 
  xlim(x_min, NA) + ylim(y_min, NA) +
  labs(x = "(Normalized) TF expression", y = "", col = "Cluster / Cell Type", title = tfp.titleCur) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  scale_color_manual(values = cols,
                     labels = vals)

ggsave("/g/scb/zaugg/deuner/figs/spearman_corr_plot_TFpeak_top1link.png", tfps1, device = "png")



###########
# PEARSON #
###########

# redoit but with all TF-peak links and putting pearson and spearman correlation coefficients in the same dataframe
# extract TF-peak from GRN ran with Spearman correlation
#tfp.links <- GRN.p@connections$TF_peaks$`0`$main # IMPORTANT LINE
col.names <- c("TF.ID",                  "TF_peak.r_bin",          "TF_peak.r",             
               "TF_peak.fdr",            "TF.name",       "peak.ID",               
               "TF_peak.fdr_direction",  "TF_peak.connectionType")
tfp.links <- GRN.p@connections$all.filtered$`0`
tfp.links <- tfp.links[,col.names]
colnames(tfp.links)[3] <- "pearson.r"
TF.ENSEMBL = GRN.s@annotation$TFs %>% dplyr::select(TF.ID, TF.ENSEMBL)
tfp.links <- inner_join(tfp.links, TF.ENSEMBL, by = "TF.ID")
names(tfp.links)[9] <- "gene.ENSEMBL"

# add gene expresion and peak accessibility counts for each connection
tfp.joined.s <- tfp.links %>% inner_join(counts.s, by = "gene.ENSEMBL") %>% inner_join(peaks.s, by = "peak.ID")

# select top5 most significant peak-gene links
## add also the pearson coefficient as a columns
tfp.top5.s <- tfp.joined.s %>% top_n(5, pearson.r)
tfp.top5.s <- tfp.top5.s[1:5, ]
tfp.top5.s.1 <- tfp.top5.s[,1:25]
tfp.top5.s.2 <- tfp.top5.s[,c(1:9, 26:41)]
# long format
tfp.top5.s.long.1 <- reshape2::melt(tfp.top5.s.1 , id.vars = colnames(tfp.top5.s)[1:9], variable.name = "cluster.RNA", value.name = "value.RNA")
tfp.top5.s.long.2 <- reshape2::melt(tfp.top5.s.2, id.vars = colnames(tfp.top5.s)[1:9], variable.name = "cluster.ATAC", value.name = "value.ATAC")
tfp.top5.s.both <- cbind(tfp.top5.s.long.1, tfp.top5.s.long.2[,10:ncol(tfp.top5.s.long.2)])
tfp.top5.s.both %>% names

# set new ids
new.ids <- rep(as.character(c(1:5)), nrow(tfp.top5.s.both)/5)
tfp.top5.s.both$ids <- new.ids


corMethod <- "pearson"

digits_round <- 3
tfp.titleCur = paste0("Peak=", tfp.top5.s.both %>% dplyr::filter(ids == "1") %>% pull(peak.ID) %>% dplyr::first(), ", TF=", tfp.top5.s.both %>% dplyr::filter(ids == "1") %>% pull(TF.ID) %>% dplyr::first(), 
                      ":\nCor = ", round(tfp.top5.s.both %>% dplyr::filter(ids == "1") %>% pull(pearson.r) %>% dplyr::first(), digits_round), " (", corMethod, "), raw p-value = ", round(tfp.top5.s.both %>% dplyr::filter(ids == "1") %>% pull(TF_peak.fdr) %>% dplyr::first(), digits_round),
                      "\n", "# eGRN link: highest correlation")

x_min = min(-0.01, min(tfp.top5.s.both$value.RNA, na.rm = TRUE))
y_min = min(-0.01, min(tfp.top5.s.both$value.ATAC, na.rm = TRUE))

tfpp1 <- ggplot(tfp.top5.s.both %>% dplyr::filter(ids == "1"), aes(x = value.RNA, y = value.ATAC, col = cluster.RNA)) + 
  geom_smooth(method = "lm", formula = "y ~ x", col = "black", na.rm = TRUE) + 
  #geom_point(size = rel(0.5)) + 
  geom_point(alpha = 0.75, size = 2) + 
  xlim(x_min, NA) + ylim(y_min, NA) +
  labs(x = "", y = "", col = "Cluster / Cell Type", title = tfp.titleCur) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  scale_color_manual(values = cols,
                     labels = vals)

ggsave("/g/scb/zaugg/deuner/figs/pearson_corr_plot_TFpeak_top1link.png", tfps1, device = "png")

###############
# FINAL PLOTS #
###############

fp <- ggarrange(pgp1, pgs1, tfpp1, tfps1, ncol = 2, nrow = 2, common.legend = TRUE, labels = c("A", "B", "C", "D"), legend = "right")
fp
ggsave("/g/scb/zaugg/deuner/figs/correlationplots.png", fp, device = "png")


counts <- counts.s
peaks <- peaks.s
head(counts)

counts.long <- tidyr::pivot_longer(counts, names(counts)[1:16], names_to = "cluster", values_to = "expression")
peaks.long <- tidyr::pivot_longer(peaks, names(peaks)[1:16], names_to = "cluster", values_to = "accessibility")


dist <- ggplot(counts.long, aes(expression)) + 
  geom_density(mapping = aes(fill = "red"), alpha = .2, col = "#333333", show.legend = TRUE) + 
  geom_density(peaks.long, mapping = aes(accessibility, fill = "#00CCFF"), alpha = .2, col = "#333333", show.legend = TRUE) + 
  xlim(c(0,25)) + 
  labs(y = "Density", x = "Activity") + 
  theme_classic() +
  #scale_fill_identity(name = 'Activity', guide = 'legend',labels = c("Gene Expression", "Peak Accessibility")) +
  scale_fill_manual(name = 'Activity', 
                      values =c('red'='red','#00CCFF'='#00CCFF'), labels = c("Gene Expression", "Peak Accessibility"))

timecourse.s@meta.data$celltype_wnn %>% levels()
# vals <- c("Diff - Immature Neuron", "Mature Neuron - small 1", "Excitatory Neuron", "Diff - hiPSC-like", "Start Diff", "Diff",
#           "Microglia 1", "Mature Neuron", "Mature Neuron - adhesion", "Neuron Development 1", "Microglia 2",
#           "hiPSC", "Diff - NPC-like", "Mature Neuron - small 2", "Neuron Development 2", "Diff - small")
nms <- c("hiPSC",                    "diff___NPC_like",          "diff",                     "diff___immature.neuron",
         "neuron___mature.neuron",   "neuron___excitatory",      "diff___hiPSC_like",        "mature.neuron___adhesion",
         "hiPSC___start.diff",       "microglia_1",              "microglia_2",              "neuron___development_1",
         "neuron___development_2",   "mature.neuron___small_1",  "mature.neuron___small_2",  "diff___tiny")
vals <- c("Diff - small", "hiPSC", "Start Diff", "Diff - hiPSC-like", "Diff - NPC-like", "Diff", "Diff - Immature Neuron", "Excitatory Neuron", 
                  "Neuron Development 1", "Neuron Development 2", "Mature Neuron - adhesion", "Mature Neuron", "Mature Neuron - small 1",
                  "Mature Neuron - small 2", "Microglia 1", "Microglia 2")
order <- c("diff___tiny", "hiPSC", "hiPSC___start.diff", "diff___hiPSC_like", "diff___NPC_like", "diff", "diff___immature.neuron", "neuron___excitatory", 
           "neuron___development_1", "neuron___development_2", "mature.neuron___adhesion", "neuron___mature.neuron", "mature.neuron___small_1",
           "mature.neuron___small_2", "microglia_1", "microglia_2")
names(cols) <- order
names(vals) <- order

ct <- DimPlot(timecourse.s, reduction = "wnn.umap", group.by = "celltype_wnn", label = F, label.size = 2.5, repel = TRUE) + 
  labs(x = "UMAP 1", y = "UMAP 2", col = "") + 
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) + 
  scale_color_manual(values = unname(cols[nms]), labels = unname(vals[nms]))

legends <- ggarrange(get_legend(dist), get_legend(pgp1), nrow=2, heights = c(0.2, 0.8), widths = c(0.5, 0.5))

# Combining plots
rm_legend <- function(p){p + theme(legend.position = "none")}
plots <- ggarrange(rm_legend(ct), rm_legend(dist), rm_legend(pgp1),
                   rm_legend(tfpp1), rm_legend(pgs1), rm_legend(tfps1), nrow = 3, ncol = 2, labels = c("A", "B", "C", "D", "E", "F"))

# plots + merged legends
afp <- ggarrange(plots, legends, widths = c(0.75, 0.25))

tiff(paste0("/g/scb/zaugg/deuner/figs/allcorrelationplots", ".tiff"), units="in", width=13, height=9, res=300, type = "cairo")
afp
dev.off()


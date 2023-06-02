########################
# PEARSON VS. SPEARMAN #
########################

library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(ggrepel)

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

## PEARSON
# TF-peak-gene connections here
GRN.p@connections$all.filtered[["0"]]
connections.p <- as.data.frame(GRN.p@connections$all.filtered[["0"]])

# SCT counts here
GRN.p@data[["RNA"]][["counts"]]
counts.p <- as.data.frame(GRN.p@data[["RNA"]][["counts"]])
counts.p$gene.ENSEMBL <- rownames(counts.p)
counts.p

# norm ATAC counts here
GRN.p@data[["peaks"]][["counts"]]
peaks.p <- as.data.frame(GRN.p@data[["peaks"]][["counts"]])
peaks.p$peak.ID <- rownames(peaks.p)
peaks.p

# connections + counts
joined.p <- connections.p %>% inner_join(counts.p, by = "gene.ENSEMBL") %>% inner_join(peaks.p, by = "peak.ID")


## SPEARMAN
# peak-gene connections here
GRN.s@connections$all.filtered[["0"]]
connections.s <- as.data.frame(GRN.s@connections$all.filtered[["0"]])

# SCT counts here
GRN.s@data[["RNA"]][["counts"]]
counts.s <- as.data.frame(GRN.s@data[["RNA"]][["counts"]])
counts.s$gene.ENSEMBL <- rownames(counts.s)
counts.s

# norm ATAC counts here
GRN.s@data[["peaks"]][["counts"]]
peaks.s <- as.data.frame(GRN.s@data[["peaks"]][["counts"]])
peaks.s$peak.ID <- rownames(peaks.s)
peaks.s

# connections + counts
joined.s <- connections.s %>% inner_join(counts.s, by = "gene.ENSEMBL") %>% inner_join(peaks.s, by = "peak.ID")

# check shared TF-peak-links
tf.peak.gene.p <- GRN.p@connections$all.filtered[["0"]]
tf.peak.gene.s <- GRN.s@connections$all.filtered[["0"]]
colnames(tf.peak.gene.p)

tf.shared.p.s <- inner_join(tf.peak.gene.p, tf.peak.gene.s, by = c("TF.ID")) 
gene.shared.p.s <- inner_join(tf.peak.gene.p, tf.peak.gene.s, by = c("gene.name"))
peak.shared.p.s <- inner_join(tf.peak.gene.p, tf.peak.gene.s, by = c("peak.ID"))

View(tf.peak.gene.p)
tf.peak.gene.s

tf.shared.p.s
gene.shared.p.s
peak.shared.p.s


# CHECK IF ASSUMPTIONS ARE HELD!
# Pearson:

# select top5 most significant peak-gene links
top5.p <- joined.p %>% top_n(5, peak_gene.r)
top5.p.1 <- top5.p[,1:32]
top5.p.2 <- top5.p[,c(1:16, 33:48)]
# long format
top5.p.long.1 <- reshape2::melt(top5.p.1 , id.vars = colnames(top5.p)[1:16], variable.name = "cluster.RNA", value.name = "value.RNA")
top5.p.long.2 <- reshape2::melt(top5.p.2, id.vars = colnames(top5.p)[1:16], variable.name = "cluster.ATAC", value.name = "value.ATAC")
top5.p.both <- cbind(top5.p.long.1, top5.p.long.2[,17:ncol(top5.p.long.2)])

# assign informative colors to the cell type clusters
cols <- c("#FF6600", "#FF0000", "#CC0033", "#99FFFF", "#0099FF", "#33FF33", "#FF00FF", "#990000", "#FF3300", "#FFCC66", "#660066", "#0000FF", "#FFFF33", "#FF6666", "#FF9933", "#66FF99")

ggplot(top5.p.both, aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Peak-gene link with highest correlation coefficient") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_label_repel(data=top5.p.both[top5.p.both$TF.name == "YY1",], aes(label=as.factor(TF.name), alpha=0.7), size=5, force=1.3, max.overlaps = 30) + 
  scale_color_manual(values = cols)

table(joined.p$TF.ID)


ggplot(top5.p.both, aes(x = value.RNA, y = value.ATAC, col = cluster.RNA, label=TF.name)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log 10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 10 most correlated peak-genes") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_text(hjust=0, vjust=0)

# set new ids
new.ids <- as.character(c(1:5))
top5.p.both$id <- new.ids
names(new.ids) <- levels(top5.p.both$id)
top5.p.both

ggplot(top5.p.both, aes(x = log10(value.RNA), y = value.ATAC, col = id)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Id", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-genes") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_label_repel(data=top5.p.both[1:5,], aes(label=as.factor(round(peak_gene.r, 4))), size=5)



p1 <- ggplot(top5.p.both %>% dplyr::filter(id == "1"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "TTop 1 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p2 <-ggplot(top5.p.both %>% dplyr::filter(id == "2"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 2 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p3 <-  ggplot(top5.p.both %>% dplyr::filter(id == "3"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 3 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p4 <-ggplot(top5.p.both %>% dplyr::filter(id == "4"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 4 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p5 <-ggplot(top5.p.both %>% dplyr::filter(id == "5"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)

ggarrange(p1,p2,p3,p4,p5)
p1
p2

# Spearman:

# select top5 most significant peak-gene links
## add also the pearson coefficient as a columns
top5.s <- joined.s %>% top_n(5, peak_gene.r)
top5.s.1 <- top5.s[,1:32]
top5.s.2 <- top5.s[,c(1:16, 33:48)]
# long format
top5.s.long.1 <- reshape2::melt(top5.s.1 , id.vars = colnames(top5.s)[1:16], variable.name = "cluster.RNA", value.name = "value.RNA")
top5.s.long.2 <- reshape2::melt(top5.s.2, id.vars = colnames(top5.s)[1:16], variable.name = "cluster.ATAC", value.name = "value.ATAC")
top5.s.both <- cbind(top5.s.long.1, top5.s.long.2[,17:ncol(top5.s.long.2)])

# set new ids
new.ids <- as.character(c(1:5))
top5.s.both$id <- new.ids
names(new.ids) <- levels(top5.s.both$id)
new.ids

ggplot(top5.s.both, aes(x = log10(value.RNA), y = value.ATAC, col = id)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Id", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-genes") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_label_repel(data=top5.s.both[1:5,], aes(label=as.factor(round(peak_gene.r, 4))), size=5)


ggplot(top5.s.both, aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Peak-gene link with highest correlation coefficient") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)

p1 <- ggplot(top5.s.both %>% dplyr::filter(id == "1"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 1 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p2 <-ggplot(top5.s.both %>% dplyr::filter(id == "2"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 2 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p3 <-  ggplot(top5.s.both %>% dplyr::filter(id == "3"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 3 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p4 <-ggplot(top5.s.both %>% dplyr::filter(id == "4"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 4 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p5 <-ggplot(top5.s.both %>% dplyr::filter(id == "5"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-gene link") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)

ggarrange(p1,p2,p3,p4,p5)
p1
p4
p5

# check links overlaps

# merge links

# redoit but with all peak-gene links and putting pearson and spearman correlation coefficients in the same dataframe
# extract peak-genes from GRN ran with Spearman correlation 
pg.links <- GRN.s@connections$peak_genes$`0` # IMPORTANT LINE
colnames(pg.links)[4] <- "spearman.r" 
# add pearson correlation values
pg.links$pearson.r <- GRN.p@connections$peak_genes$`0`$peak_gene.r

# add id column
pg.links$id <- as.character(c(1:nrow(pg.links)))

# add gene expresion and peak accessibility counts for each connection
pg.joined.s <- pg.links %>% inner_join(counts.s, by = "gene.ENSEMBL") %>% inner_join(peaks.p, by = "peak.ID")

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
new.ids <- as.character(c(1:5))
names(new.ids) <- levels(pg.top5.s.both$id)
new.ids

ggplot(pg.top5.s.both, aes(x = log10(value.RNA), y = value.ATAC, col = id)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Id", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-genes") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5))


p1 <- ggplot(pg.top5.s.both %>% dplyr::filter(id == "316796"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 1 most correlated peak-gene link. Id: 316796") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p2 <-ggplot(pg.top5.s.both %>% dplyr::filter(id == "380465"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 2 most correlated peak-gene link. Id: 380465") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p3 <-  ggplot(pg.top5.s.both %>% dplyr::filter(id == "484736"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 3 most correlated peak-gene link. Id: 484736") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p4 <-ggplot(pg.top5.s.both %>% dplyr::filter(id == "714900"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 4 most correlated peak-gene link. Id: 714900") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p5 <-ggplot(pg.top5.s.both %>% dplyr::filter(id == "714906"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "SPEARMAN CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-gene link. Id: 714906") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)

ggarrange(p1,p2,p3,p4,p5)
p1

###########
# PEARSON #
###########

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
new.ids <- as.character(c(1:5))
names(new.ids) <- levels(pg.top5.s.both$id)
new.ids

ggplot(pg.top5.s.both, aes(x = log10(value.RNA), y = value.ATAC, col = id)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Id", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-genes") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5))


p1 <- ggplot(pg.top5.s.both %>% dplyr::filter(id == "146686"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "TTop 1 most correlated peak-gene link. Id: 146686") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p2 <-ggplot(pg.top5.s.both %>% dplyr::filter(id == "18914"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 2 most correlated peak-gene link. Id: 18914") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p3 <-  ggplot(pg.top5.s.both %>% dplyr::filter(id == "441623"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 3 most correlated peak-gene link. Id: 441623") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p4 <-ggplot(pg.top5.s.both %>% dplyr::filter(id == "634342"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 4 most correlated peak-gene link. Id: 634342") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)
p5 <-ggplot(pg.top5.s.both %>% dplyr::filter(id == "757121"), aes(x = log10(value.RNA), y = value.ATAC, col = cluster.RNA)) + 
  geom_point(alpha = 0.75, size = 2) + 
  labs(x = "Gene Expression (SCT) log10", y = "Peak Accessibility", col = "Cluster / Cell Type", title = "PEARSON CORRELATION: Peak-Gene links", caption = "Top 5 most correlated peak-gene link. Id: 757121") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = cols)

ggarrange(p1,p2,p3,p4,p5)
p1

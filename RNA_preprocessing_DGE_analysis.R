###############################################################
# BULK RNA-SEQ TIMECOURSE DATA PREPROCESSING AND DGE ANALYSIS #
###############################################################

## Source of the data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118307
## Paper: https://www.sciencedirect.com/science/article/pii/S2405471218303582?via%3Dihub

# Load packages
library(edgeR)  # 'edgeR' also loads 'limma' 
library(EnsDb.Hsapiens.v86)
library(Glimma)
library(SummarizedExperiment)
library(factoextra)
library(pheatmap)
library(readxl)

# Set data directory
dir <- "/g/scb/zaugg/deuner/GRaNPA/inputdata/"

# Load raw bulk RNA counts
RNA = read_xlsx(paste0(dir, "GSE118307_RNA-Seq_timecourse_Supplement.xlsx"), col_names = TRUE, skip = 1)

#Data Cleaning
RNA <- as.data.frame(RNA[!duplicated(RNA$Gene.ID), ])  # keep only the first (arbitrary)
rownames(RNA) <- RNA$Gene.ID
RNA[,1] <- NULL
head(RNA)

#Subset the Wild Type cohorts (remove KO)
RNA <- RNA %>%
  dplyr::select(starts_with("WT"))

# Obtain gene annotations
genes <- genes(EnsDb.Hsapiens.v86)

# Build SummarizedExperiment
se <- SummarizedExperiment(assay = list("counts" = RNA),
                           colData = colnames(RNA))
names(colData(se)) <- "sample"
se$day <- as.factor(rep(0:4, each = 7))

# Remove lowly expressed genes
keep <- filterByExpr(se, group = se)
se <- se[keep, ]

# Save the se in a .RData file
save(se, file = paste0(dir, "se.RData"))
load(paste0(dir, "se.RData")) #load it for further use

# Compute TMM normalization factors
dgl <- calcNormFactors(se, method = "TMM")

# Calculate CPM (Counts Per Million Reads)
assays(se)$CPM <- cpm(se)

# Calculate CPM based on effective library size (using TMM normalization factors)
assays(se)$TMM <- cpm(dgl, normalized.lib.sizes = TRUE)

PCA <- prcomp(log2(t(assays(se)$TMM)+1))
fviz_pca_ind(PCA, addEllipses = F,
             col.ind = se$sample,
             pointsize = 3) 
fviz_pca_ind(PCA, addEllipses = F,
             col.ind = se$day,
             pointsize = 3)

# heatmap
pheatmap(log2(assays(se)$TMM + 1),
         show_rownames = FALSE, 
         annotation_col = as.data.frame(colData(se)))

###############################################################
# BULK RNA-SEQ TIMECOURSE DATA PREPROCESSING AND DGE ANALYSIS #
###############################################################

## Source of the data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118307
## Paper: https://www.sciencedirect.com/science/article/pii/S2405471218303582?via%3Dihub

# Done it for: Day4 vs. Day0, Day4 vs. Day2, Day2 vs. Day0.

# Load packages
library(edgeR)  # 'edgeR' also loads 'limma' 
library(EnsDb.Hsapiens.v86)
library(Glimma)
library(limma)
library(SummarizedExperiment)
library(factoextra)
library(pheatmap)
library(readxl)
library(dplyr)

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

# # Remove lowly expressed genes (I want DGE data for all genes)
# keep <- filterByExpr(se, group = se)
# se <- se[keep, ]

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

# # heatmap
# pheatmap(log2(assays(se)$TMM + 1),
#          show_rownames = FALSE, 
#          annotation_col = as.data.frame(colData(se)))

# Create design matrix and contrast matrix of iPSCs (day 0) vs mature neurons (day 4)
d4d0_se <- se[,se$day %in% c(4, 0)]
day <- as.factor(d4d0_se$day)
design <- model.matrix(~0 + day)
#colnames(design)[1:2] <- levels(day)
head(design)

contrast <- makeContrasts(
  day4 - day0,
  levels = colnames(design)
)
contrast

# Remove mean-variance relationship from count data
y <- voom(dgl[, dgl$samples$day %in% c(4,0)], design, plot = TRUE)


# Fit linear models for all pairwise comparisons
fit <- lmFit(y, design)
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)

# Examine number of Differentially Expressed genes
results <- topTable(fit, sort.by = "logFC", n = Inf) # Get all results, ordered by P-value
head(results)

results.5fdr <-  subset(results, adj.P.Val < 0.05) 
nrow(results.5fdr)

# Compare the number of differentially (up and down) expressed genes
dt <- decideTests(fit, adjust.method = "fdr", p.value = 0.05)
summary(dt)

# glMD Plot
#glMDPlot(fit, status = dt, side.main="ENSEMBL", counts = log2(assays(se)$TMM+1), groups = se$day)

# Volcano Plot
dgl$samples$group <- dgl$samples$day  # Allows coloring by cohort
volcanoplot(fit, dge = dgl, style = "p-value", highlight = 0, names = fit$genes$ID, hl.col="blue",
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35) 

# load gene ENSEMBL annotations
require('biomaRt')
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'hgnc_symbol',
    'ensembl_gene_id',
    'gene_biotype'),
  uniqueRows = TRUE)
colnames(annotLookup)[1:2] <- c("gene", "gene.ENSEMBL") 

# Generate DE dataframe to give as input to GRaNPA
results$gene <- vapply(strsplit(rownames(results),"[.]"), `[`, 1, FUN.VALUE=character(1))
DE <- results %>%
  dplyr::rename(padj = adj.P.Val) %>%
  dplyr::inner_join(annotLookup, by = "gene") %>%
  dplyr::select(ENSEMBL = gene.ENSEMBL, padj, logFC)

fwrite(DE, paste0(dir, "DE_timecourse_Day4vsDay0.tsv"))

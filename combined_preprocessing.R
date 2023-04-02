#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# COMBINED DATASET PREPROCESSING #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load libraries
library(Seurat)
library(qs)
library(sctransform)
library(ggplot2)
library(BiocParallel)
library(scDblFinder)
library(dplyr)

# Set up working directory
work.dir <- "/g/scb/zaugg/deuner/GRaNIE/"

# Load merged Seurat Object
combined.path <- "/g/scb/zaugg/deuner/SCENIC+/inputdata/" 
combined.s <- qread(paste0(combined.path, "multiome.combined.seuratObject.qs"))
combined.s

# Create new seurat Object starting with the RNA counts
comb.s <- CreateSeuratObject(combined.s@assays$RNA@counts, assay = "RNA")
rm(combined.s)

##  QC and selecting cells for further analysis

# Compute % of mitochondrial genes for each cell
comb.s[["percent.mt"]] <- PercentageFeatureSet(comb.s, pattern = "^MT-")
comb.s[["percent.ribo"]] <- PercentageFeatureSet(comb.s, pattern = "^RP[SL]")

# Visualize QC metrics as a violin plot
VlnPlot(comb.s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Visualize some more metrics
plot1 <- FeatureScatter(comb.s, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(comb.s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_hline(yintercept=11500, linetype="dashed", color = "darkred")
plot1 + plot2

# Filter out low quality cells
comb.s <- subset(comb.s, subset = nFeature_RNA > 500 & nFeature_RNA < 11500 & percent.mt < 15)

# Run doublet finder
# https://github.com/plger/scDblFinder
comb.s <- scDblFinder(as.SingleCellExperiment(comb.s))
comb.s <- as.Seurat(comb.s)

# filter out doublets
comb.s <- comb.s %>% subset(scDblFinder.class == "singlet")

# Normalization
comb.s <- SCTransform(object=comb.s, method = "glmGamPoi", return.only.var.genes = FALSE, conserve.memory = TRUE, vars.to.regress = c("percent.mt", "percent.ribo"))

# Find highly variable genes
comb.s <- FindVariableFeatures(comb.s, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(comb.s), 10)

# plot variable features
plot <- VariableFeaturePlot(comb.s)
plot <- LabelPoints(plot = plot, points = top10, repel = TRUE)
plot

# Linear dimensionality reduction (PCA)
comb.s <- RunPCA(comb.s, features = VariableFeatures(object = comb.s))

# Examine and visualize PCA results a few different ways
print(comb.s[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(comb.s, dims = 1:2, reduction = "pca")

# Plot PCA
Idents(comb.s) <- comb.s@meta.data$ident
DimPlot(comb.s, reduction = "pca")

# Visualize PCs
DimHeatmap(comb.s, dims = 1:15, cells = 500, balanced = TRUE)

# Elbow plot
ElbowPlot(comb.s)

# Cluster cells
comb.s <- FindNeighbors(comb.s, dims = 1:15)
comb.s <- FindClusters(comb.s, resolution = 0.5)

# UMAP
comb.s <- RunUMAP(comb.s, dims = 1:15, umap.method = "umap-learn")
DimPlot(comb.s, reduction = "umap")
DimPlot(comb.s, reduction = "umap", group.by = "orig.ident")
uwot <- RunUMAP(comb.s, dims = 1:15, umap.method = "uwot")
DimPlot(object = uwot, reduction = 'umap') + NoLegend()
DimPlot(object = uwot, reduction = 'umap', group.by = "orig.ident") + NoLegend()

# Save object
qsave(comb.s, paste0(work.dir, "tmp/combined.pp.seuratObject.qs"))

# Read object
comb.s <- qread(paste0(work.dir, "tmp/combined.pp.seuratObject.qs"))

# visualize technical variables
FeaturePlot(comb.s, reduction = "umap", features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"))

# Find marker genes (differentially expressed)
# find markers for every cluster compared to all remaining cells, report only the positive ones
comb.markers <- FindAllMarkers(comb.s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
comb.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

comb.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(comb.s, features = top10$gene) + NoLegend()

# Cell type Annotation
comb.s@meta.data$num_nCount_RNA <- as.numeric(comb.s@meta.data$nCount_RNA)
DimPlot(comb.s, reduction = "umap", group.by = "num_nCount_RNA")

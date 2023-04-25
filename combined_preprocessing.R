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
library(monocle3)
library(harmony)
library(Signac)

#######################################
# NOTES ON THE STANDARD PREPROCESSING #
#######################################

# ATAC (/g/scb/zaugg/marttine/dariaMultiome/ArchR/dariaMultiome.R)
# 1. Mapping/fragment file generation with chromap against Hg38
# 2. Bin-cell matrix generation with ArchR
# 3. Filter for > 1000 fragments, > 0.12 Fraction fragments in promoter, TSS score > 3 
# 4. TF-IDF normalisation
# 5. Iterative LSI for dimension reduction
# 
# RNA (/g/scb/zaugg/marttine/dariaMultiome/RNA/dariaMultiomeRNA.R)
# 1. Mapping with STARsolo against Hg38
# 2. Filtering with Seurat (nFeatures > 500, mito %, ribo %)
# 3. Normalisation with SCT
# 4. Doublet calling/removal with DoubletFinder
# 5. Save seurat objects “/g/scb/zaugg/marttine/dariaMultiome/RNA/{sample}-seuratObject.qs”
# 
# Find cell barcodes present in both RNA and ATAC modalities and save “/g/scb/zaugg/marttine/dariaMultiome/cellBarcodes.rna.atac.{sample}.txt“
# 
# ATAC continued (/g/scb/zaugg/marttine/dariaMultiome/ArchR/dariaMultiome.R)
# 1. Re-run iterative LSI
# 2. Louvain clustering
# 3. Peak calling on called clusters (q < 0.01)
# 4. Save peak matrices “/g/scb/zaugg/marttine/dariaMultiome/peakMatrix.{sample}.mtx” (unnormalised and non-binary)
# 5. Save peak ranges “/g/scb/zaugg/marttine/dariaMultiome/peakRanges.{sample}.qs”
# 6. Save ArchR output/objects “/g/scb/zaugg/marttine/dariaMultiome/ArchR/{sample}/
#   
#   Combined (/g/scb/zaugg/marttine/dariaMultiome/combined/multiome.R)
# 1. Export peak matrix and read into Seurat object with RNA
# 2. Redo PCA (RNA)
# 3. Redo TF-IDF and SVD (ATAC)
# 4. Joint embedding with WNN
# 5. Save Seurat objects /g/scb/zaugg/marttine/dariaMultiome/combined/{sample}.seuratObject.qs

# Set up working directory
work.dir <- "/g/scb/zaugg/deuner/GRaNIE/"

# Load merged Seurat Object
combined.path <- "/g/scb/zaugg/deuner/SCENIC+/inputdata/" 
comb.s <- qread(paste0(combined.path, "multiome.combined.seuratObject.qs"))
comb.s

#######################
# RNA REPREPROCESSING #
#######################

# Set RNA data as default assay
DefaultAssay(comb.s) <- "RNA"

##  QC and selecting cells for further analysis

# Compute % of mitochondrial genes for each cell
comb.s[["percent.mt"]] <- PercentageFeatureSet(comb.s, pattern = "^MT-")
comb.s[["percent.ribo"]] <- PercentageFeatureSet(comb.s, pattern = "^RP[SL]")

# Visualize QC metrics as a violin plot
Idents(comb.s) <- comb.s@meta.data$orig.ident
VlnPlot(comb.s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2)

# Visualize some more metrics
plot1 <- FeatureScatter(comb.s, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_hline(yintercept=15, linetype="dashed", color = "darkred")
plot2 <- FeatureScatter(comb.s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_hline(yintercept=11500, linetype="dashed", color = "darkred")
plot3 <- FeatureScatter(comb.s, feature1 = "nCount_RNA", feature2 = "percent.ribo") + geom_hline(yintercept=15, linetype="dashed", color = "darkred")
plot1 + plot2 + plot3

# Filter out low quality cells
comb.s <- subset(comb.s, subset = nFeature_RNA > 500 & nFeature_RNA < 11500 & percent.mt < 15)

# Run doublet finder
# https://github.com/plger/scDblFinder
sc <- as.SingleCellExperiment(comb.s)
sc@colData$samples <- sc@colData$orig.ident
sc@colData
comb.s <- scDblFinder(as.SingleCellExperiment(comb.s), samples = "orig.ident")
comb.s <- as.Seurat(comb.s)

DefaultAssay(comb.s) <- "RNA"

# filter out doublets
comb.s <- comb.s %>% subset(scDblFinder.class == "singlet")

# Normalization - SCTransform regressing out mt and ribo variables
mt.genes <-  grep(pattern = "^MT-", x = rownames(x = comb.s), value = TRUE)
ribo.genes <-  grep(pattern = "^RP[SL]", x = rownames(x = comb.s), value = TRUE)
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

# Visualize possible batch effects
FeaturePlot(comb.s, reduction = "pca", features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"))

# Visualize PCs
DimHeatmap(comb.s, dims = 1:15, cells = 500, balanced = TRUE)

# Elbow plot
ElbowPlot(comb.s, ndims = 40)

# Decide optimal number of PCs that explain most of the variation in the data
## We can calculate where the principal components start to elbow by taking the larger value of:
## The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
## The point where the percent change in variation between the consecutive PCs is less than 0.1%.
## We will start by calculating the first metric:

# Determine percent of variation associated with each PC
pct <- comb.s[["pca"]]@stdev / sum(comb.s[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 #41

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 #18
# Usually, we would choose the minimum of these two metrics as the PCs covering the majority of the variation in the data.

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs #18

# Cluster cells
comb.s <- FindNeighbors(comb.s, dims = 1:18)
comb.s <- FindClusters(comb.s, resolution = 0.5)
comb.s[["res.0.5"]] <- Idents(comb.s)

# UMAP
#comb.s <- RunUMAP(comb.s, dims = 1:18, umap.method = "umap-learn")
comb.s<- RunUMAP(comb.s, umap.method = "umap-learn", dims = 1:18, min.dist = 0.5, spread = 1, negative.sample.rate = 5, 
                 n.epochs = NULL, n.components = 2, learning.rate = 1) 
DimPlot(comb.s, reduction = "umap")
DimPlot(comb.s, reduction = "umap", group.by = "orig.ident")

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

##########################################
# Cell type Annotation based on RNA data #
##########################################

## From cell type to cluster annotation
# hiPSC markers
hiPSC.markers <- unique(c("5T4", "ABCG2", "B18R", "CBX2", "CD9", "CD30", "TNFRSF8", "CD117", "CDX2", 
                          "CHD1", "DNMT3B", "DPPA2", "DPPA4", "DPPA5", "ESG1", "EPCAM", "TROP1", "NR3B2",
                          "ERVMER34", "ESGP", "FBXO15", "FGF4", "FGF5", "FOXD3", "GCNF", "NR6A1", "GDF3",
                          "POU5F1", "OCT4", "409B2", "B7", "HOIK1", "KUCG2", "WIBJ2", "WTC", "KLF4, NANOG"))
FeaturePlot(comb.s, reduction = "umap", features = hiPSC.markers) + labs(caption = "hiPSC markers")

# NPC markers
NPC.markers <- unique(c("HES5", "GFAP", "NES", "VIM", "SOX2", "BLBP", "PAX6", "Nestin, ABCG2", "FGFR1"))
FeaturePlot(comb.s, reduction = "umap", features = NPC.markers) + labs(caption = "NPC markers")

# mature neurons markers
mat.neuron.markers <- unique(c("NeuN", "MAP2", "PSD95", "NEFM", "NEFH", "SYP", "SLC17A6", "VGLUT2", "RBFOX3", "SAP90", "DLG4"))
FeaturePlot(comb.s, reduction = "umap", features = mat.neuron.markers) + labs(caption = "mature neurons markers")

# immature neurons markers
immat.neuron.markers <- unique(c("DCX", "TBR1", "PSANCAM", "Ascl1", "ASCL1", "DCX", "TUBB3", "NeuroD1", "NEUROD1",
                                 "TBR1", "STMN1","NGN2"))
FeaturePlot(comb.s, reduction = "umap", features = immat.neuron.markers) + labs(caption = "immature neurons markers")

# general neuronal markers
neuron.markers <- unique(c("NEUN", "MAP2", "PSD95", "vGLUT2", "NMDAR2B", "GAP43", "GAT1", "GAD65", "GAD67", "GIRK2", 
                           "NURR1", "LMX1B", "FOXA2")) 
FeaturePlot(comb.s, reduction = "umap", features = neuron.markers) + labs(caption = "neuronal markers")

# microglia markers
microglia.markers <- unique(c("TMEM119", "CD11B", "CD45", "IBA1", "F80", "CD68", "CD40",  "CD11B", "P2RY12", "CD14", "CD80", 
                              "CD115", "CX3CR1", "F4", "F40", "FCER1G", "TMEM119", "CD11b", "CD45", "Iba1", "CX3CR1", "F4", 
                              "F80", "CD68", "CD40"))
FeaturePlot(comb.s, reduction = "umap", features = microglia.markers) + labs(caption = "microglia markers")

# differentiating cells markers
diff.markers <- unique(c("AADC", "DAT", "LMX1B", "MAP2", "SOX2", "SOX1", "MAP2A", "MAP2B", "TUJ1", "GFAP", "GABA", "GAD65", 
                         "TH", "CHAT"))
FeaturePlot(comb.s, reduction = "umap", features = diff.markers) + labs(caption = "differentiating cells markers")

DimPlot(comb.s, reduction = "umap", label = TRUE) + DimPlot(comb.s, reduction = "umap", group.by = "orig.ident")

# Define basic cell type clusters
basic.cluster.ids <- c("NPC-1", #0
                       "NPC-2", #1
                       "neuron-1", #2
                       "NPC-3", #3
                       "NPC-5", #4 
                       "neuron-2", #5
                       "neuron-3", #6
                       "diff-NPC", #7
                       "hiPSC-1", #8
                       "diff", #9
                       "neuron-NPC-like", #10
                       "microglia", #11
                       "diff-neuron", #12
                       "hiPSC-2", #13
                       "neuron-4", #14
                       "neuron-5" #15
)
Idents(comb.s)

# Add new cell annotation to the Seurat Object
names(basic.cluster.ids) <- levels(comb.s)
comb.s <- RenameIdents(comb.s, basic.cluster.ids)
comb.s <- comb.s %>% AddMetaData(metadata = Idents(comb.s), col.name = "celltype")
DimPlot(comb.s, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "celltype")
DimPlot(comb.s, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "orig.ident")
DimPlot(comb.s, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "res.0.5")

# Save object
qsave(comb.s, paste0(work.dir, "tmp/combined.pp.seuratObject.qs"))

# Read object
comb.s <- qread(paste0(work.dir, "tmp/combined.pp.seuratObject.qs"))

####################################
# TRAJECTORY - PSEUDOTIME ANALYSIS #
####################################

# Use monocle3
## https://github.com/satijalab/seurat/issues/2833

# Extract data from seurat object
gene_annotation <- as.data.frame(rownames(comb.s@reductions[["pca"]]@feature.loadings),
                                 row.names = rownames(comb.s@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"


cell_metadata <- as.data.frame(comb.s@assays[["SCT"]]@counts@Dimnames[[2]],
                               row.names = comb.s@assays[["SCT"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"

New_matrix <- comb.s@assays[["SCT"]]@counts
New_matrix <- New_matrix[rownames(comb.s@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix

# create cell data object
cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)


# Add cluster information
recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

list_cluster <- comb.s@active.ident
names(list_cluster) <- comb.s@assays[["SCT"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <-comb.s@reductions[["umap"]]@cell.embeddings

#cds_from_seurat@preprocess_aux$gene_loadings <- timecourse.s@reductions[["pca"]]@feature.loadings

# create trajectory inference graph
cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = F)

# visualize tranjectory analysis
plot_cells(cds_from_seurat, 
           color_cells_by = "cluster",
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=4)+ 
  ggtitle("Combined - Trajectory Analysis") + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5))

# compute pseudotime
cds_from_seurat <- order_cells(cds_from_seurat, reduction_method = "UMAP", root_cells = colnames(cds_from_seurat[,clusters(cds_from_seurat) == "hiPSC-1"]))

# visualize pseudotime
plot_cells(cds_from_seurat,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph=FALSE,
           graph_label_size=1.5) + 
  ggtitle("Combined - Pseudotime") + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5))

cds_from_seurat$monocle3_pseudotime <- pseudotime(cds_from_seurat) 
cds_from_seurat$seurat_clusters <- comb.s@meta.data$seurat_clusters
data.pseudo <-as.data.frame(colData(cds_from_seurat))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = seurat_clusters)) + 
  geom_boxplot() + 
  ylab("ordered clusters by pseudotime")

# find genes that change as function of pseudotime
genes <- graph_test(cds_from_seurat, neighbor_graph = "principal_graph", cores = 4)

genes %>% 
  arrange(q_value) %>%
  filter(status == "OK") %>%
  head()

gene.names <- genes %>% 
  arrange(q_value) %>%
  filter(status == "OK") %>%
  head() %>%
  rownames()

FeaturePlot(comb.s, features = gene.names)

# visualize pseudotime in seurat
comb.s$pseudotime <- pseudotime(cds_from_seurat)
FeaturePlot(comb.s, features = "pseudotime", label = TRUE)

# save the seurat object
Idents(comb.s) <- comb.s[["celltype"]]
qsave(comb.s, paste0(work.dir, "tmp/combined.pp.seuratObject.qs"))

# Read object
comb.s <- qread(paste0(work.dir, "tmp/combined.pp.seuratObject.qs"))


########################
# ATAC REPREPROCESSING #
########################

DefaultAssay(comb.s) <- "ATAC"

comb.s@assays$ATAC@data

## Quality Control (QC)

# Visualize QC metrics as a violin plot
Idents(comb.s) <- comb.s@meta.data$orig.ident
VlnPlot(comb.s, features = c("nFeature_ATAC", "nCount_ATAC", "percent.mt", "percent.ribo"), ncol = 2)
# Visualize some more metrics
plot1 <- FeatureScatter(comb.s, feature1 = "nCount_ATAC", feature2 = "percent.mt") + geom_hline(yintercept=15, linetype="dashed", color = "darkred")
plot2 <- FeatureScatter(comb.s, feature1 = "nCount_ATAC", feature2 = "nFeature_ATAC") + geom_hline(yintercept=11500, linetype="dashed", color = "darkred")
plot1 + plot2

## Normalisation and Dimensionality Reduction
comb.s@assays$ATAC@data <- comb.s@assays$ATAC@counts # remove prenormalized data
comb.s <- RunTFIDF(comb.s)
comb.s <- FindTopFeatures(comb.s, min.cutoff = 10)
comb.s <- RunSVD(comb.s)
comb.s

## Non-linear dimensionality reduction (UMAP)
# We exclude the first dimension as this is typically correlated with sequencing depth
DepthCor(comb.s)
# Visualize ATAC space
comb.s <- RunUMAP(comb.s, reduction = "lsi", dims = 2:20, reduction.name = "umap.atac", reduction.key = "atacUMAP_", umap.method = "umap-learn", assay = "ATAC")
p1 <- DimPlot(comb.s, reduction = "umap.atac", group.by = "celltype", label = TRUE) +  ggtitle("ATAC")
p2 <- DimPlot(comb.s, reduction = "umap.atac", group.by = "orig.ident", label = FALSE)  + ggtitle("ATAC")
p1 + p2

p3 <- DimPlot(comb.s, reduction = "umap", group.by = "orig.ident", label = TRUE) +  ggtitle("RNA")
p4 <- DimPlot(comb.s, reduction = "umap.atac", group.by = "orig.ident", label = FALSE)  + ggtitle("ATAC")
p3 + p4

# Check whether any technical variable drives the UMAP representation
FeaturePlot(comb.s, reduction = "umap.atac", features = c("nCount_ATAC", "nFeature_ATAC", "percent.mt", "percent.ribo"))

# Make UMAP on a subset of regions that are at least covered in all the datasets (as clustering is based on dataset)
#comb.s.2 <- comb.s
#comb.s.s <- subset(comb.s.2, rownames(comb.s.2) %in% )

comb.s.2 <- RunTFIDF(comb.s.2)
#comb.s.2 <- FindTopFeatures(comb.s.2, min.cutoff = "q0", assay = "ATAC")
comb.s.2 <- subset(comb.s.2, rownames(comb.s.2))
comb.s.2 <- RunSVD(comb.s.2)
comb.s.2 <- RunUMAP(comb.s.2, reduction = "lsi", dims = 2:20, reduction.name = "umap.atac", reduction.key = "atacUMAP_", umap.method = "umap-learn")
# Correct for batch effect using RunHarmony()
# comb.s <- RunHarmony(object = comb.s, group.by.vars = 'orig.ident', reduction = 'lsi', assay.use = 'peak', project.dim = FALSE, dims.use = 2:20)
# comb.s <- RunUMAP(comb.s, reduction = "harmony", dims = 2:20, reduction.name = "umap.atac.harmony", reduction.key = "atacUMAP_", umap.method = "umap-learn")
# p1 <- DimPlot(comb.s, reduction = "umap.atac.harmony", group.by = "celltype", label = TRUE) +  ggtitle("ATAC")
# p2 <- DimPlot(comb.s, reduction = "umap.atac", group.by = "orig.ident", label = FALSE)  + ggtitle("ATAC")
# p1 + p2
# 
# p3 <- DimPlot(comb.s, reduction = "umap", group.by = "celltype", label = TRUE) +  ggtitle("RNA")
# p4<- DimPlot(comb.s, reduction = "umap.atac", group.by = "celltype", label = FALSE)  + ggtitle("ATAC")
# p3 + p4

############################
# JOINT EMBEDDING WITH WNN #
############################

comb.s <- FindMultiModalNeighbors(comb.s, reduction.list = list("pca", "lsi"), dims.list = list(1:18, 2:20))
comb.s <- RunUMAP(comb.s, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
comb.s <- FindClusters(comb.s, graph.name = "wsnn", algorithm = 3, verbose = FALSE) #SLM algorithm

# Show UMAPs for RNA, ATAC and WNN
p1 <- DimPlot(comb.s, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(comb.s, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(comb.s, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

DimPlot(comb.s, reduction = "wnn.umap", group.by = "orig.ident", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN") + 
  DimPlot(comb.s, reduction = "wnn.umap", group.by = "wsnn_res.0.8", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

p1 <- DimPlot(comb.s, reduction = "umap", group.by = "orig.ident", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(comb.s, reduction = "umap.atac", group.by = "orig.ident", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(comb.s, reduction = "wnn.umap", group.by = "orig.ident", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))


# Clustering with WNN 
comb.s <- FindClusters(comb.s, graph.name = "wsnn", algorithm = 3, verbose = FALSE, res = c(0.1, seq(0.25, 1, 0.25), seq(2,10,1), seq(12,20,2))) #SLM algorithm
p1 <- DimPlot(comb.s, reduction = "wnn.umap", group.by = "wsnn_res.0.5", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p2 <- DimPlot(comb.s, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN - RNA Cell Annotation")
p1 | p2

# Adapting Cell annotation performed with RNA clusters to WNN clusters
ann.equival <- c("NPC-1", #0
                 "NPC-2", #1
                 "diff-NPC", #2
                 "neuron-2", #3
                 "NPC-3", #4
                 "neuron-3", #5
                 "NPC-4", #6
                 "diff", #7
                 "neuron-NPC-like", #8
                 "hiPSC-1", #9
                 "NPC-5", #10
                 "neuron-1", #11
                 "NPC-6", #12
                 "neuron-5", #13
                 "neuron-6", #14
                 "hiPSC-2",  #15
                 "diff-neuron",  #16
                 "microglia-1",  #17
                 "NPC-7",  #18
                 "microglia-2",  #19
                 "neuron-7",  #20
                 "neuron-4",  #21
                 "neuron-8",  #22
                 "neuron-9",  #23
                 "undefined",  #24
                 "NPC-8"  #25
)

Idents(comb.s)
Idents(comb.s) <- "wsnn_res.0.5"

names(ann.equival) <- levels(comb.s)
comb.s <- RenameIdents(comb.s, ann.equival)
comb.s[["celltype_wnn"]] <- Idents(comb.s)

p1 <- DimPlot(comb.s, reduction = "wnn.umap", group.by = "wsnn_res.0.5", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p2 <- DimPlot(comb.s, reduction = "wnn.umap", group.by = "celltype_wnn", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN - RNA Cell Annotation")
p1 | p2

# Save seurat object
qsave(comb.s, paste0(work.dir, "tmp/combined.pp.seuratObject.qs"))

# Read object
comb.s <- qread(paste0(work.dir, "tmp/combined.pp.seuratObject.qs"))
FeaturePlot(comb.s, reduction = "umap", features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"))

###########################
# REMOVE MICROGLIAL CELLS #
###########################

# Read the seurat
comb.s <- qread(paste0(work.dir, "tmp/combined.pp.seuratObject.qs"))

DimPlot(comb.s, reduction = "wnn.umap", group.by = "wsnn_res.0.5", label = TRUE, label.size = 2.5, repel = TRUE) + theme(legend.position = "none")
DimPlot(comb.s, reduction = "wnn.umap", group.by = "celltype_wnn", label = TRUE, label.size = 2.5, repel = TRUE) + theme(legend.position = "none")

# remove microglia and very small clusters
comb.s2 <- comb.s %>%
  subset(subset = wsnn_res.0.5 != 24) %>%
  subset(subset = wsnn_res.0.5 != 19) %>%
  subset(subset = wsnn_res.0.5 != 17) %>%
  subset(subset = celltype_wnn != "NPC-8")
  

# check if I removed the cells correctly
table(comb.s2$celltype_wnn)


# take a look at the UMAP without those cell types
DimPlot(comb.s2, reduction = "wnn.umap", group.by = "celltype_wnn", label = TRUE, label.size = 2.5, repel = TRUE) + theme(legend.position = "none")

# Save seurat object
qsave(comb.s2, paste0(work.dir, "tmp/combined.pp.nomicro.seuratObject.qs"))

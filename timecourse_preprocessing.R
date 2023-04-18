##############################
# TIMECOURSE REPREPROCESSING #
##############################

# Load libraries
library(GRaNIE)
library(Seurat)
library(qs)
library(dplyr)
library(reticulate)
library(ggplot2)
library(monocle3)

# Set working directory path
work.dir = "/g/scb/zaugg/deuner/GRaNIE/"

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


# load "raw" seurat object
timecourse.s <- qread(paste0(work.dir, "inputdata/timecourse.seuratObject.qs"))

# visualize the preprocessed data 
# Visualize QC metrics as a violin plot
Idents(timecourse.s) <- timecourse.s@meta.data$orig.ident
VlnPlot(timecourse.s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Define QC thresholds
nfeatures.threshold <- 11500
percent.mt.threshold <- 15

# Visualize some more metrics
plot1 <- FeatureScatter(timecourse.s, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  geom_hline(yintercept=percent.mt.threshold, linetype="dashed", color = "darkred")
plot2 <- FeatureScatter(timecourse.s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept=nfeatures.threshold, linetype="dashed", color = "darkred")
plot3 <- FeatureScatter(timecourse.s, feature1 = "nCount_RNA", feature2 = "percent.ribo") + 
  geom_hline(yintercept=15, linetype="dashed", color = "darkred")
plot1 + plot2 + plot3

# Remove some outliers (just 3)
timecourse.s <- subset(timecourse.s, subset = nFeature_RNA < nfeatures.threshold)
timecourse.s <- subset(timecourse.s, subset = percent.mt < percent.mt.threshold)

# Non random distribution of percent.mt, normalise with SCTtransform excluding this variable
timecourse.s <- SCTransform(object=timecourse.s, assay = "RNA", method = "glmGamPoi", return.only.var.genes = FALSE, conserve.memory = TRUE, vars.to.regress = "percent.mt")

# Feature selection
timecourse.s <- FindVariableFeatures(timecourse.s, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(timecourse.s), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(timecourse.s)
plot1 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1

# Save seurat object 
qsave(timecourse.s, paste0(work.dir, "tmp/timecourse.norm.seuratObject.qs"))

# Read seurat object
timecourse.s <- qread(paste0(work.dir, "tmp/timecourse.norm.seuratObject.qs"))

# PCA
DimPlot(timecourse.s, reduction = "pca")

# Visualize if data clusters by some technical variable
FeaturePlot(timecourse.s, reduction = "pca", features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"))

# Analyze PCs
VizDimLoadings(timecourse.s, dims = 1:4, reduction = "pca")
ElbowPlot(timecourse.s)

# decide optimal clustering algorithm for modularity optimization and play with UMAPs
## 1 = original Louvain algorithm; 
## 2 = Louvain algorithm with multilevel refinement; 
## 3 = SLM algorithm; 
## 4 = Leiden algorithm

timecourse.s <- FindNeighbors(timecourse.s, k.param = 10, dims = 1:15, annoy.metric = "euclidean")

# louvain.umap <-FindClusters(timecourse.s, res = 0.5, algorithm = 1, graph.name = "SCT_snn") %>% RunUMAP(dims = 1:15, umap.method = "umap-learn") %>% DimPlot(reduction = "umap")
# louvain.mlr.umap <- FindClusters(timecourse.s, res = 0.5, algorithm = 2, graph.name = "SCT_nn") %>% RunUMAP(dims = 1:15, umap.method = "umap-learn") %>% DimPlot(reduction = "umap")
# slm.umap <- FindClusters(timecourse.s, res = 0.5, algorithm = 3, graph.name = "SCT_nn") %>% RunUMAP(dims = 1:15, umap.method = "umap-learn") %>% DimPlot(reduction = "umap")
# leiden.umap <- FindClusters(timecourse.s, res = 0.5, algorithm = 4, graph.name = "SCT_nn") %>% RunUMAP(dims = 1:15, umap.method = "umap-learn") %>% DimPlot(reduction = "umap")


# louvain.umap
# louvain.mlr.umap
# slm.umap
# leiden.umap

# decide optimal number of dimension by doing the proper PCA analysis and decide optimal cluster resolution
timecourse.s <-FindClusters(timecourse.s, res = 0.5, algorithm = 1, graph.name = "SCT_nn")
timecourse.s <- RunUMAP(timecourse.s, umap.method = "umap-learn", dims = 1:15, min.dist = 0.5, spread = 1, negative.sample.rate = 5, 
                        n.epochs = NULL, n.components = 2, learning.rate = 1) 
DimPlot(timecourse.s, reduction = "umap", label = TRUE)

# Save seurat object
qsave(timecourse.s, paste0(work.dir, "tmp/timecourse.norm.seuratObject.qs"))

# Read seurat object
timecourse.s <- qread(paste0(work.dir, "tmp/timecourse.norm.seuratObject.qs"))

# Visualize if data clusters by some technical variable
FeaturePlot(timecourse.s, reduction = "umap", features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"))

## Cell Type annotation

##################################################
# from cluster markers to to celltype annotation #
##################################################

markers <- FindAllMarkers(timecourse.s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) -> top3

markers %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = avg_log2FC) -> top1

markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) -> top2

DoHeatmap(timecourse.s, features = top10$gene) + NoLegend()

print(top3, n = 50)

DimPlot(timecourse.s, reduction = "umap", label = TRUE)

FeaturePlot(timecourse.s, reduction = "umap", features = top1$gene)

DimPlot(timecourse.s, reduction = "umap", label = TRUE)

FeaturePlot(timecourse.s, reduction = "umap", features = "POU4F1")

names(new.cluster.ids.manual.1) <- levels(timecourse.s)
timecourse.s <- RenameIdents(timecourse.s, new.cluster.ids)
DimPlot(timecourse.s, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

###############################################
# from celltype markers to cluster annotation #
###############################################

# BDNF
FeaturePlot(timecourse.s, reduction = "umap", features = "POU4F1") + labs(title = "BDNF")

# immature neuron
FeaturePlot(timecourse.s, reduction = "umap", features = c("DCX", "TBR1"))

# mature neuron
FeaturePlot(timecourse.s, reduction = "umap", features = c("RBFOX3", "MAP2", "SYP"))

# gluneuron
FeaturePlot(timecourse.s, reduction = "umap", features = c("SLC17A7", "SLC17A6"))

# NPC
FeaturePlot(timecourse.s, reduction = "umap", features = c("GFAP", "NES", "VIM", "SOX2", "BLBP", "PAX6"))

# immature neuron
FeaturePlot(timecourse.s, reduction = "umap", features = c("DCX", "PSANCAM"))

# mature neuron
FeaturePlot(timecourse.s, reduction = "umap", features = c("NeuN", "TUBB3", "MAP2"))

# iPSC
FeaturePlot(timecourse.s, reduction = "umap", features = c("NANOG", "KLF4"))

# hiPSC
FeaturePlot(timecourse.s, reduction = "umap", features = c("POU5F1", "OCT4"))  
# NPC
FeaturePlot(timecourse.s, reduction = "umap", features = c("HES5"))  
# neural maturation
FeaturePlot(timecourse.s, reduction = "umap", features = c("SLC17A6", "VGLUT2"))

# hiPSC and immature neurons
FeaturePlot(timecourse.s, reduction = "umap", features = c("TBR2", "MASH1", "Ascl1", "ASCL1", "DCX", "TUBB3", "NeuroD1", "NEUROD1",
                                                           "TBR1", "STMN1","NGN2"))
# mature neurons
FeaturePlot(timecourse.s, reduction = "umap", features = c("NeunN", "MAP2", "PSD95", "NEFM", "NEFH", "SYP")) 
# glutamatergic neurons
FeaturePlot(timecourse.s, reduction = "umap", features = c("VGLUT1", "VGLUT2", "NMDAR1", "NMDAR2B", "GLS", "GLUL")) 
# GABAnergic neurons
FeaturePlot(timecourse.s, reduction = "umap", features = c("GAT1", "SLC6A1", "GABRA1", "GABRA2", "GAD65", "GAD67"))
# microglia
FeaturePlot(timecourse.s, reduction = "umap", features = c("TMEM119", "CD11b", "CD45", "Iba1", "CX3CR1", "F4", "F80", "CD68", "CD40")) 

# ipsc
FeaturePlot(timecourse.s, reduction = "umap", features = c("409B2", "B7", "HOIK1", "KUCG2", "WIBJ2", "WTC"))

## ALL
# microglia
FeaturePlot(timecourse.s, reduction = "umap", features = c("TMEM119", "CD11B", "CD45", "IBA1", "F80", "CD68", "CD40",  "CD11B", "P2RY12", "CD14", "CD80", "CD115", "CX3CR1", "F4", "F40", "FCER1G"))
# neuron
FeaturePlot(timecourse.s, reduction = "umap", features = c("NEUN", "MAP2", "PSD95", "vGLUT2", "NMDAR2B", "GAP43", "GAT1", "GAD65", "GAD67", "GIRK2", "NURR1", "LMX1B", "FOXA2"))
# hiPSC
FeaturePlot(timecourse.s, reduction = "umap", features = c("5T4", "ABCG2", "B18R", "CBX2", "CD9", "CD30", "TNFRSF8", "CD117", "CDX2", 
                                                           "CHD1", "DNMT3B", "DPPA2", "DPPA4", "DPPA5", "ESG1", "EPCAM", "TROP1", "NR3B2",
                                                           "ERVMER34", "ESGP", "FBXO15", "FGF4", "FGF5", "FOXD3", "GCNF", "NR6A1", "GDF3"))
# diff
FeaturePlot(timecourse.s, reduction = "umap", features = c("AADC", "DAT", "LMX1B", "MAP2", "SOX2", "SOX1", "MAP2A", "MAP2B", "TUJ1", "GFAP", "GABA", "GAD65", "TH", "CHAT"))

#----------------------------------------------------------

# hiPSC markers
hiPSC.markers <- unique(c("5T4", "ABCG2", "B18R", "CBX2", "CD9", "CD30", "TNFRSF8", "CD117", "CDX2", 
                          "CHD1", "DNMT3B", "DPPA2", "DPPA4", "DPPA5", "ESG1", "EPCAM", "TROP1", "NR3B2",
                          "ERVMER34", "ESGP", "FBXO15", "FGF4", "FGF5", "FOXD3", "GCNF", "NR6A1", "GDF3",
                          "POU5F1", "OCT4", "409B2", "B7", "HOIK1", "KUCG2", "WIBJ2", "WTC", "KLF4, NANOG"))
FeaturePlot(timecourse.s, reduction = "umap", features = hiPSC.markers) + labs(caption = "hiPSC markers")

# NPC markers
NPC.markers <- unique(c("HES5", "GFAP", "NES", "VIM", "SOX2", "BLBP", "PAX6", "Nestin, ABCG2", "FGFR1"))
FeaturePlot(timecourse.s, reduction = "umap", features = NPC.markers) + labs(caption = "NPC markers")

# mature neurons markers
mat.neuron.markers <- unique(c("NeuN", "MAP2", "PSD95", "NEFM", "NEFH", "SYP", "SLC17A6", "VGLUT2", "RBFOX3", "SAP90", "DLG4"))
FeaturePlot(timecourse.s, reduction = "umap", features = mat.neuron.markers) + labs(caption = "mature neurons markers")

# immature neurons markers
immat.neuron.markers <- unique(c("DCX", "TBR1", "PSANCAM", "Ascl1", "ASCL1", "DCX", "TUBB3", "NeuroD1", "NEUROD1",
                                 "TBR1", "STMN1","NGN2"))
FeaturePlot(timecourse.s, reduction = "umap", features = immat.neuron.markers) + labs(caption = "immature neurons markers")

# general neuronal markers
neuron.markers <- unique(c("NEUN", "MAP2", "PSD95", "vGLUT2", "NMDAR2B", "GAP43", "GAT1", "GAD65", "GAD67", "GIRK2", 
                           "NURR1", "LMX1B", "FOXA2")) 
FeaturePlot(timecourse.s, reduction = "umap", features = neuron.markers) + labs(caption = "neuronal markers")

# microglia markers
microglia.markers <- unique(c("TMEM119", "CD11B", "CD45", "IBA1", "F80", "CD68", "CD40",  "CD11B", "P2RY12", "CD14", "CD80", 
                              "CD115", "CX3CR1", "F4", "F40", "FCER1G", "TMEM119", "CD11b", "CD45", "Iba1", "CX3CR1", "F4", 
                              "F80", "CD68", "CD40"))
FeaturePlot(timecourse.s, reduction = "umap", features = microglia.markers) + labs(caption = "microglia markers")

# differentiating cells markers
diff.markers <- unique(c("AADC", "DAT", "LMX1B", "MAP2", "SOX2", "SOX1", "MAP2A", "MAP2B", "TUJ1", "GFAP", "GABA", "GAD65", 
                         "TH", "CHAT"))
FeaturePlot(timecourse.s, reduction = "umap", features = diff.markers) + labs(caption = "differentiating cells markers")

# cells treated with BDNF
FeaturePlot(timecourse.s, reduction = "umap", features = "POU4F1") + labs(caption = "Cells treated with BDNF")


##############################################
# establish a cell type annotation consensus #
##############################################

DimPlot(timecourse.s, reduction = "umap", label = TRUE)

new.cluster.ids <- c("diff - NPC-like", #0
                     "hiPSC", #1
                     "diff", #2
                     "diff - immature.neuron", #3
                     "neuron - mature.neuron", #4 #Gap43
                     "microglia", #5
                     "diff - hiPSC-like", #6
                     "hiPSC - start.diff", #7
                     "neuron - development", #8
                     "neuron - excitatory", #9
                     "mature.neuron - adhesion", #10 #Map2
                     "neuron - mature.neuron", #11
                     "diff - RUNX1", #12
                     "diff - tiny" #13
)

basic.cluster.ids <- c("diff", #0
                       "hiPSC", #1
                       "diff", #2
                       "diff", #3
                       "neuron", #4 #Gap43
                       "microglia", #5
                       "diff", #6
                       "hiPSC", #7
                       "neuron", #8
                       "neuron", #9
                       "neuron", #10 #Map2
                       "neuron", #11
                       "diff", #12
                       "diff" #13
)

# Add new cell annotation to the Seurat Object
names(new.cluster.ids) <- levels(timecourse.s)
timecourse.s <- RenameIdents(timecourse.s, new.cluster.ids)
DimPlot(timecourse.s, reduction = "umap", label = T, pt.size = 0.5)

timecourse.s <- timecourse.s %>% AddMetaData(metadata = Idents(timecourse.s), col.name = "celltype")
names(basic.cluster.ids) <- new.cluster.ids
timecourse.s <- RenameIdents(timecourse.s, basic.cluster.ids)
timecourse.s <- timecourse.s %>% AddMetaData(metadata = Idents(timecourse.s), col.name = "basic_celltype")

Idents(timecourse.s) <- timecourse.s@meta.data$celltype

p1 <- DimPlot(timecourse.s, reduction = "umap", group.by = "basic_celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("Basic Cell Type Annotation") +
  theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(timecourse.s, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("Accurate Cell Type Annotation") +
  theme(plot.title = element_text(hjust = 0.5))
p1 | p2

# Save seurat object
qsave(timecourse.s, paste0(work.dir, "tmp/timecourse.norm.seuratObject.qs"))

# Read seurat object
timecourse.s <- qread(paste0(work.dir, "tmp/timecourse.norm.seuratObject.qs"))

####################################
# TRAJECTORY - PSEUDOTIME ANALYSIS #
####################################

# Use monocle3
## https://github.com/satijalab/seurat/issues/2833

# Extract data from seurat object
gene_annotation <- as.data.frame(rownames(timecourse.s@reductions[["pca"]]@feature.loadings),
                                 row.names = rownames(timecourse.s@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"


cell_metadata <- as.data.frame(timecourse.s@assays[["SCT"]]@counts@Dimnames[[2]],
                               row.names = timecourse.s@assays[["SCT"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"


New_matrix <- timecourse.s@assays[["SCT"]]@counts
New_matrix <- New_matrix[rownames(timecourse.s@reductions[["pca"]]@feature.loadings), ]
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



list_cluster <- timecourse.s@active.ident
names(list_cluster) <- timecourse.s@assays[["SCT"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster



cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <-timecourse.s@reductions[["umap"]]@cell.embeddings


#cds_from_seurat@preprocess_aux$gene_loadings <- timecourse.s@reductions[["pca"]]@feature.loadings

# create trajectory inference graph
cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = F)

# visualize trajectory analysis
plot_cells(cds_from_seurat, 
           color_cells_by = "cluster",
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=4)+ 
  ggtitle("Timecourse - Trajectory Analyses") + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5))

# compute pseudotime
cds_from_seurat <- order_cells(cds_from_seurat, reduction_method = "UMAP", root_cells = colnames(cds_from_seurat[,clusters(cds_from_seurat) == "hiPSC"]))

# visualize pseudotime
plot_cells(cds_from_seurat,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph=FALSE,
           graph_label_size=1.5) +
  ggtitle("Timecourse - Pseudotime") + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5))

cds_from_seurat$monocle3_pseudotime <- pseudotime(cds_from_seurat) 
cds_from_seurat$seurat_clusters <- timecourse.s@meta.data$seurat_clusters
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

FeaturePlot(timecourse.s, features = gene.names)

# visualize pseudotime in seurat
timecourse.s$pseudotime <- pseudotime(cds_from_seurat)
Idents(timecourse.s) <- timecourse.s$celltype
FeaturePlot(timecourse.s, features = "pseudotime", label = TRUE)

Idents(timecourse.s) <- timecourse.s$basic_celltype
FeaturePlot(timecourse.s, features = "pseudotime", label = TRUE)

Idents(timecourse.s) <- timecourse.s$celltype

# cluster cells according to pseudotime
timecourse.s$pseudotime_clusters_n7 <- cut(timecourse.s$pseudotime, c(0, 10, 20, 30, 40, 50, 60, Inf), c(1:7), include.lowest = TRUE)
timecourse.s$pseudotime_clusters_n14 <- cut(timecourse.s$pseudotime, seq(0, 70, by = 5), c(1:14), include.lowest = TRUE)

# visualize them
DimPlot(timecourse.s, reduction = "umap", group.by = "pseudotime_clusters_n7", label = TRUE, label.size = 2.5, repel = TRUE) +
  theme(plot.title = element_text(hjust = 0.5)) +  ggtitle("Pseudotime Clustering")
DimPlot(timecourse.s, reduction = "umap", group.by = "pseudotime_clusters_n14", label = TRUE, label.size = 2.5, repel = TRUE) +
  theme(plot.title = element_text(hjust = 0.5)) +  ggtitle("Pseudotime Clustering")

# Save the seurat object
qsave(timecourse.s, paste0(work.dir, "tmp/timecourse.pp.seuratObject.qs"))

# Read the seurat object
timecourse.s <- qread(paste0(work.dir, "tmp/timecourse.pp.seuratObject.qs"))


########################
# ATAC REPREPROCESSING # (needed because we have filtered out some cells)
########################

# set the ATAC assay as the default assay
DefaultAssay(timecourse.s) <- "ATAC"

## Quality Control (QC)

# Visualize QC metrics as a violin plot
Idents(timecourse.s) <- timecourse.s[["orig.ident"]]
VlnPlot(timecourse.s, features = c("nFeature_ATAC", "nCount_ATAC", "percent.mt", "percent.ribo"), ncol = 2)

## Normalisation and Dimensionality Reduction
Idents(timecourse.s) <- timecourse.s[["celltype"]]
timecourse.s <- RunTFIDF(timecourse.s)
timecourse.s <- FindTopFeatures(timecourse.s, min.cutoff = "q0")
timecourse.s <- RunSVD(timecourse.s)

## Non-linear dimensionality reduction (UMAP)
# We exclude the first dimension as this is typically correlated with sequencing depth
DepthCor(timecourse.s)
# Visualize ATAC space
timecourse.s <- RunUMAP(timecourse.s, reduction = "lsi", dims = 2:20, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
p1 <- DimPlot(timecourse.s, reduction = "umap.atac", group.by = "celltype", label = TRUE) +  ggtitle("ATAC")
p2 <- DimPlot(timecourse.s, reduction = "umap.atac", group.by = "orig.ident", label = FALSE)  + ggtitle("ATAC")
p1 + p2

p3 <- DimPlot(timecourse.s, reduction = "umap", group.by = "celltype", label = FALSE) +  ggtitle("RNA")
p4 <- DimPlot(timecourse.s, reduction = "umap.atac", group.by = "celltype", label = FALSE)  + ggtitle("ATAC")
p3 + p4

# Check whether any technical variable drives the UMAP representation
FeaturePlot(timecourse.s, reduction = "umap.atac", features = c("nCount_ATAC", "nFeature_ATAC", "percent.mt", "percent.ribo"))


############################
# JOINT EMBEDDING WITH WNN # (with new RNA and ATAC preprocessed data)
############################

timecourse.s <- FindMultiModalNeighbors(timecourse.s, reduction.list = list("pca", "lsi"), dims.list = list(1:20, 2:20))
timecourse.s <- RunUMAP(timecourse.s, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
timecourse.s <- FindClusters(timecourse.s, graph.name = "wsnn", algorithm = 3, verbose = FALSE, res = 0.8) #SLM algorithm

# Show UMAPs for RNA, ATAC and WNN
p1 <- DimPlot(timecourse.s, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(timecourse.s, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(timecourse.s, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

DimPlot(timecourse.s, reduction = "wnn.umap", group.by = "orig.ident", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN") + 
  DimPlot(timecourse.s, reduction = "wnn.umap", group.by = "wsnn_res.0.8", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

# Save seurat object
qsave(timecourse.s, paste0(work.dir, "tmp/timecourse.pp.seuratObject.qs"))

# Read the seurat object
timecourse.s <- qread(paste0(work.dir, "tmp/timecourse.pp.seuratObject.qs"))

# Clustering with WNN 
timecourse.s <- FindClusters(timecourse.s, graph.name = "wsnn", algorithm = 3, verbose = FALSE, res = 0.5) #SLM algorithm
p1 <- DimPlot(timecourse.s, reduction = "wnn.umap", group.by = "wsnn_res.0.5", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p2 <- DimPlot(timecourse.s, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN - RNA Cell Annotation")
p1 | p2

# Adapting Cell annotation performed with RNA clusters to WNN clusters
ann.equival <- c("hiPSC", #0
                 "diff - NPC-like", #1
                 "diff", #2
                 "diff - immature.neuron", #3
                 "neuron - mature.neuron", #4
                 "neuron - excitatory", #5
                 "diff - hiPSC-like", #6
                 "mature.neuron - adhesion", #7
                 "hiPSC - start.diff", #8
                 "microglia-1", #9
                 "microglia-2", #10
                 "neuron - development-1", #11
                 "neuron - development-2", #12
                 "mature.neuron - small-1", #13
                 "mature.neuron - small-2", #14
                 "diff - tiny"  #15
)

names(ann.equival) <- levels(timecourse.s)
timecourse.s <- RenameIdents(timecourse.s, ann.equival)
timecourse.s[["celltype_wnn"]] <- Idents(timecourse.s)

# visualize it
p1 <- DimPlot(timecourse.s, reduction = "wnn.umap", group.by = "wsnn_res.0.5", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p2 <- DimPlot(timecourse.s, reduction = "wnn.umap", group.by = "celltype_wnn", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN Cell Annotation")
p1 | p2
p2

# Add sampleID column 
timecourse.s$sampleID <- timecourse.s$orig.ident

# Set SCT to active assay
DefaultAssay(timecourse.s) <- "SCT"

# Save seurat object
qsave(timecourse.s, paste0(work.dir, "tmp/timecourse.pp.seuratObject.qs"))

# Read the seurat
timecourse.s <- qread(paste0(work.dir, "tmp/timecourse.pp.seuratObject.qs"))

###########################
# REMOVE MICROGLIAL CELLS #
###########################

# Read the seurat
timecourse.s <- qread(paste0(work.dir, "tmp/timecourse.pp.seuratObject.qs"))

DimPlot(timecourse.s, reduction = "wnn.umap", group.by = "wsnn_res.5", label = TRUE, label.size = 2.5, repel = TRUE) + theme(legend.position = "none")

# remove microglia and very small clusters
timecourse.s2 <- timecourse.s %>%
  subset(subset = wsnn_res.5 != 12) %>%
  subset(subset = wsnn_res.5 != 22) %>%
  subset(subset = wsnn_res.5 != 29) %>%
  subset(subset = wsnn_res.5 != 37) %>%
  subset(subset = wsnn_res.5 != 35) %>%
  subset(subset = wsnn_res.5 != 40)

# check if I removed the cells correctly
table(timecourse.s2$celltype_wnn)

# remove robust remaining cells
timecourse.s2 <- timecourse.s2 %>%
  subset(subset = celltype_wnn != "mature.neuron - small-1") %>%
  subset(subset = celltype_wnn != "diff - tiny")

# check again if I removed the cells correctly
table(timecourse.s2$celltype_wnn)

# take a look at the UMAP without those cell types
DimPlot(timecourse.s2, reduction = "wnn.umap", group.by = "celltype_wnn", label = TRUE, label.size = 2.5, repel = TRUE) + theme(legend.position = "none")

# visualize pseudotime clustering
DimPlot(timecourse.s2, reduction = "umap", group.by = "pseudotime_clusters_n7", label = TRUE, label.size = 2.5, repel = TRUE) +
  theme(plot.title = element_text(hjust = 0.5)) +  ggtitle("Pseudotime Clustering")
DimPlot(timecourse.s2, reduction = "umap", group.by = "pseudotime_clusters_n14", label = TRUE, label.size = 2.5, repel = TRUE) +
  theme(plot.title = element_text(hjust = 0.5)) +  ggtitle("Pseudotime Clustering")

# get markers for each pseudotime cluster
## n = 7
Idents(timecourse.s2) <- "pseudotime_clusters_n7"
markers7 <- FindAllMarkers(timecourse.s2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers7 %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC) -> top1_7
top1_7
DotPlot(timecourse.s2, features = top1_7$gene) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
DoHeatmap(timecourse.s2, features = unique(top1_7$gene))

markers7 %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC) -> top3_7
top3_7
DotPlot(timecourse.s2, features = unique(top3_7$gene)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
DoHeatmap(timecourse.s2, features = unique(top3_7$gene))

## n = 14
Idents(timecourse.s2) <- "pseudotime_clusters_n14"
markers14 <- FindAllMarkers(timecourse.s2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers14 %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC) -> top1_14
top1_14
DotPlot(timecourse.s2, features = top1_14$gene)
DoHeatmap(timecourse.s2, features = unique(top1_14$gene))

# set Idents to celltype again
Idents(timecourse.s2) <- "celltype_wnn"

# Save seurat object
qsave(timecourse.s2, paste0(work.dir, "tmp/timecourse.pp.nomicro.seuratObject.qs"))

# Read the seurat
timecourse.s2 <- qread(paste0(work.dir, "tmp/timecourse.pp.nomicro.seuratObject.qs"))

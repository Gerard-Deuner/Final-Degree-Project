#########
# PANDO #
#########

# https://quadbiolab.github.io/Pando/

# Load Packages
library(Pando)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(qs)
library(Signac)
library(dplyr)
library(tidyr)
library(purrr)
library(EnsDb.Hsapiens.v86)

# Define dataset  
dataset <- "timecourse" # timecourse | combined

# Set data directory
data.dir <- paste0("/g/scb/zaugg/deuner/GRaNIE/tmp/", dataset, ".pp.nomicro.seuratObject.qs")

# Load seurat object
seurat_object <- qread(data.dir)

# Get chromatin accessibility data
fragpath <- '/g/scb/zaugg/deuner/SCENIC+/inputdata/timecourse_fragments.tsv'
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- 'UCSC'
 
#seurat_object <- CreateSeuratObject(counts = rna.s.obj@assays$RNA@counts ,assay = 'RNA')
# atac = CreateChromatinAssay(
#   counts = seurat_object@assays$ATAC@counts,
#   sep = c(':', '-'),
#   fragments = fragpath,
#   annotation = annotation
# )
# seurat_object[['peaks']] = atac

# Get motif data
data(motifs)

# Add genome annotation for ChromatinAssay
DefaultAssay(seurat_object) <- "ATAC"
genome(seurat_object) <- "hg38"

# Add annotation
## get peak ranges (timecourse)
seqlevelsStyle(annotation) <- 'UCSC'
peakRanges <- qread("/g/scb/zaugg/marttine/dariaMultiome/peakRanges.timecourse.qs")
Annotation(seurat_object[["ATAC"]]) <- peakRanges


# Initiate GRN object and select candidate regions
seurat_object <- initiate_grn(seurat_object, 
                              peak_assay = "ATAC",
                              rna_assay = "RNA",
                              exclude_exons = T)

# Scan candidate regions for TF binding motifs
seurat_object <- find_motifs(
  seurat_object,
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# Code that fixes an error ecountered in infer_grn(). It insures a one-to-one mapping between cand_regions and peak_ranges.
regions <- NetworkRegions(seurat_object)
cand_ranges <- regions@ranges #%>% as.data.frame() %>% dplyr::rename("seqnames" = "names") %>% makeGRangesFromDataFrame()
peak_ranges <- StringToGRanges(rownames(GetAssay(seurat_object, assay='ATAC')))
peak_overlaps <- findOverlaps(cand_ranges, peak_ranges)

peak_overlaps_df <- data.frame(peak_overlaps)

peak_overlaps_df <- peak_overlaps_df[!duplicated(peak_overlaps_df[['queryHits']]), ]
peak_overlaps_df <- peak_overlaps_df[!duplicated(peak_overlaps_df[['subjectHits']]), ]

new_motif_data <- regions@motifs@data[peak_overlaps_df[['queryHits']], ]
new_regions_peaks <- peak_overlaps_df[['subjectHits']]
new_ranges_data <- regions@ranges[peak_overlaps_df[['queryHits']]]

seurat_object@grn@regions@motifs@data <- new_motif_data
seurat_object@grn@regions@peaks <- new_regions_peaks
seurat_object@grn@regions@ranges <- new_ranges_data

# Select variable features
seurat_object <- Seurat::FindVariableFeatures(seurat_object, assay='RNA', nfeatures = 100)

# Infer gene regulatory network
seurat_object <- infer_grn(seurat_object,
                           parallel = F,
                           verbose = 2,
                           genes = VariableFeatures(seurat_object))

# Print inferred coefficients
coef(seurat_object)

# Find gene and regulatory modules 
seurat_object <- find_modules(seurat_object)

# Print modules
modules <- NetworkModules(seurat_object)

# Save Network output
write.csv(modules@meta, paste0("/g/scb/zaugg/deuner/Pando/outputdata/", dataset, ".GRN.tsv"))

# counts <- Read10X_h5('../data/pbmc/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5')
# counts$Peaks=counts$Peaks[startsWith(rownames(counts$Peaks),'chr'),]
# 
# fragpath <- '../data/pbmc/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz'
# annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(annotation) <- 'UCSC'
# 
# pbmc <- CreateSeuratObject(counts = counts$`Gene Expression`,assay = 'RNA')
# atac = CreateChromatinAssay(
#   counts = counts$Peaks,
#   sep = c(':', '-'),
#   fragments = fragpath,
#   annotation = annotation
# )
# pbmc[['peaks']] = atac
# 
# pbmc_grn = initiate_grn(pbmc)
# pbmc_grn = find_motifs(pbmc_grn, pfm=motifs, genome=BSgenome.Hsapiens.UCSC.hg38)
# 
# pbmc_grn <- FindVariableFeatures(pbmc_grn, nfeatures=100)
# 
# pbmc_grn = inferd_grn(pbmc_grn, genes=VariableFeatures(pbmc_grn), parallel=F, verbose = 2)
# pbmc_grn = find_modules(pbmc_grn)

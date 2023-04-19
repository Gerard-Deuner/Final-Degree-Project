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

# Define dataset  
dataset <- "timecourse" # timecourse | combined

# Set data directory
data.dir <- paste0("/g/scb/zaugg/deuner/GRaNIE/tmp/", dataset, ".pp.seuratObject.qs")

# Load seurat object
seurat_object <- qread(data.dir)

# Get motif data
data(motifs)

# Add genome annotation for ChromatinAssay
DefaultAssay(seurat_object) <- "ATAC"
genome(seurat_object) <- "hg38"

# Add annotation
## get peak ranges (timecourse)
peakRanges <- qread("/g/scb/zaugg/marttine/dariaMultiome/peakRanges.timecourse.qs")
Annotation(seurat_object[["ATAC"]]) <- peakRanges

# Select variable features
seurat_object <- Seurat::FindVariableFeatures(seurat_object, assay='RNA')

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

# Code that fixes an error ecountered in infer_grn(). It inseured a one-to-one mapping between cand_regions and peak_ranges.
regions <- NetworkRegions(seurat_object)
cand_ranges <- peakRanges@ranges %>% as.data.frame() %>% dplyr::rename("seqnames" = "names") %>% makeGRangesFromDataFrame()
peak_ranges <- StringToGRanges(rownames(GetAssay(seurat_object, assay='ATAC')))
peak_overlaps <- findOverlaps(cand_ranges, peak_ranges)

peak_overlaps_df <- data.frame(peak_overlaps)

peak_overlaps_df <- peak_overlaps_df[!duplicated(peak_overlaps_df[['queryHits']]), ]
peak_overlaps_df <- peak_overlaps_df[!duplicated(peak_overlaps_df[['subjectHits']]), ]

new_motif_data <- peakRanges@motifs@data[peak_overlaps_df[['queryHits']], ]
new_regions_peaks <- peak_overlaps_df[['subjectHits']]
new_ranges_data <- peakRanges@ranges[peak_overlaps_df[['queryHits']]]

seurat_object@grn@regions@motifs@data <- new_motif_data
seurat_object@grn@regions@peaks <- new_regions_peaks
seurat_object@grn@regions@ranges <- new_ranges_data

# Infer gene regulatory network
seurat_object <- infer_grn(seurat_object,
                           peak_to_gene_method = 'Signac',
                           method = 'glm',
                           parallel = TRUE)

# Print inferred coefficients
coef(seurat_object)

# Find gene and regulatory modules 
test_srt <- find_modules(test_srt)

# Print modules
modules <- NetworkModules(test_srt)

# Save Network output
write.csv(modules@meta, paste0("/g/scb/zaugg/deuner/Pando/outputdata/", dataset, ".GRN.tsv"))

overlaps <- findOverlaps(seurat_object@assays$ATAC@ranges, seurat_object@assays$ATAC@ranges)
all(data.frame(overlaps)[['queryHits']] == data.frame(overlaps)[['subjectHits']])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Data Extraction for each Dataset Specific Seurat Object separately #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load libraries
library(qs)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(tibble)

# Load Seurat Objects 
timecourse <- qread("/g/scb/zaugg/marttine/dariaMultiome/combined/timecourse.seuratObject.qs")
neuron <- qread("/g/scb/zaugg/marttine/dariaMultiome/combined/Neuron.seuratObject.qs")
NPC <- qread("/g/scb/zaugg/marttine/dariaMultiome/combined/NPC.seuratObject.qs")
cocultured <- qread("/g/scb/zaugg/marttine/dariaMultiome/combined/cocultured28.seuratObject.qs")

Seurat.objects <- c(timecourse, neuron, NPC, cocultured)

samples <- c("timecourse", "neuron", "NPC", "cocultured")

for (i in 1:4){
  
  # Set up sample name and its corresponding Seurat Object
  sample = samples[i]    
  object = Seurat.objects[[i]]

  print(paste("Running on sample", sample))
  
  # Export RNA raw counts matrix
  write.table(as.matrix(object@assays$RNA@counts), 
              paste0('../inputdata/',sample,'_RNA_counts.csv'), 
              sep = ',', row.names = T, col.names = T, quote = F)
  
  print(paste("RNA counts stored for", sample))
  
  # Export ATAC-seq peak counts matrix
  write.table(as.matrix(object@assays$ATAC@counts), 
              paste0('../inputdata/',sample,'_peak_counts.csv'), 
              sep = ',', row.names = T, col.names = T, quote = F)
  
  print(paste("Peak counts stored for", sample))
  
  # Export cell type annotations
  write.table(as.matrix(Cells(object)), 
              paste0('../inputdata/',sample,'_metadata.csv'), 
              sep = ',', row.names = T, col.names = T, quote = F)
  
  print(paste("metadata stored for", sample))
}

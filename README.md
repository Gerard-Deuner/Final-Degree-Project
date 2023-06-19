# Benchmarking & Evaluating Single-cell Enhancer-Gene Regulatory Networks 

### Final Degree Project | Bachelor's Degree in Bioinformatics | ESCI-UPF

## Project Abstract

Enhancer-inclusive gene regulatory networks allow us to characterize and provide functional understanding of regulatory interactions underlying phenotypes in a cell-type specific manner by the reconstruction of tripartite networks which describe transcription factor-enhancer-gene associations. GRaNIE is a method that infers gene regulatory networks including enhancers from bulk transcriptomics (RNA-seq) and chromatin accessibility (ATAC-seq) data. Bulk data presents constraints that mask cell-type- and cell-state-specific gene expression and chromatin accessibility data as the counts represent an average of the activity of a population of cells. Hence, to overcome these limitations the approach requires its adjustment at the single-cell resolution to take into account inter-cell variation and consequently truly capture specific cell-type and cell-state regulatory differences; Therefore, a proper benchmark at the single-cell layer of the approach remains to be accomplished. 

We propose a way of benchmarking GRaNIE at the single-cell resolution through a validation-based approach that assesses the network's ability to predict real regulatory interactions from cell-type-specific promoter capture Hi-C, eQTL and ChIP-Seq data. Furthermore, we report further validation through a network biological comprehensive analysis alongside an evaluation of the eGRNs predictive power by means of GRaNPA, a machine learning framework that assesses the network's capability of  predicting cell-type-specific differential gene expression data. The results show GRaNIE, applied to single-cell hiPSC-neuron timecourse datasets, is capable of inferring neuronal differentiation and development relevant regulatory associations as well as predicting cell-type specific iPSC against neuron DGE at the single-cell level. 

## Scripts Descriptions
### Datasets Preprocessing
* #### :page_facing_up: **timecourse_preprocessing.R**
  Preprocessing performed on the timecourse dataset.
* #### :page_facing_up: **combined_preprocessing.R**
  Preprocessing performed on the combined dataset.
* #### :page_facing_up: **data_extraction.R**
  Extraction of raw data (RNA counts, ATAC counts, metadata) from the standardly preprocessed dataset-specific-seurat objects for further dataset-customized preprocessing.
* #### :page_facing_up: **SeuratToAnndata.R**
  Conversion from seurat objects to AnnData objects to permit the usage of the preprocessed datasets for eGRN inference via SCENIC+ which is only available in Python/Scanpy.
* #### :page_facing_up: **get_fragments_from_archr.R**
  Get the ATAC fragments from the ArchR projects for each dataset.

### Single-cell GRaNIE
* #### :page_facing_up: **GRaNIE_batch_mode.R**
  Run GRaNIE in batch mode (across a wide range of cluster resolutions on the integrated WNN space) generating an eGRN per resolution.
  #### :page_facing_up: **GRaNIE_specific_resolution.R**
  Run GRaNIE in for a specific cluster resolution.
  #### :page_facing_up: **GRaNIE_helper_functions.R**
  Helper functions to execute GRaNIE for single-cell datasets. Functions to preprocess the data, run GRaNIE for a specific metadata column and in batch mode.
  #### :page_facing_up: **GRaNIE_metadata.R**
  Run GRaNIE pseudobulking the cells according a metadata column e.g. celltype. 
  #### :page_facing_up: **timecourse_scGRaNIE_test.R**
  Run GRaNIE on the timecourse dataset for initial size tests, different cell type annotations and batch mode are utilized.
  #### :page_facing_up: **network_enrichment_analysis.R**
  Perform gene enrichment and TF enrichment analyses on a given eGRN.
  #### :page_facing_up: **custom_peakgene_QC_plots.R**
  Customize GRaNIE original code that creates the peak-gene quality control plots in a qay that different eGRNs can be analyzed at the same time.   
  #### :page_facing_up: **MergePeakGeneQCPDFs.sh**
  Create a single pdf document with the final peak-gene quality control plots for different eGRNs (each network-specific QC plot in a different page)


### single-cell GRaNIE Benchmark
#### Correlation Method Analysis - *Pearson and Spearman correlation comparison on GRaNIE at the single-cell level.*
#### Overlap Analysis - *Overlap analysis of GRaNIE across resolutions.*
#### Validations - *Validation analyses on GRaNIE networks using promoter capture Hi-C, eQTL and ChIP-Seq data.*
#### GRaNPA - *Validation analyses on GRaNIE networks using GRaNPA.*
### single cell-eGRN Inference Methods
#### SCENIC+
#### Pando



**The scripts descriptions can be found within each respective directory.** 



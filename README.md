# Benchmarking & Evaluating Single-cell Enhancer-Gene Regulatory Networks 

## Final Degree Project | Bachelor's Degree in Bioinformatics | ESCI-UPF

### Project Abstract

Enhancer-inclusive gene regulatory networks allow us to characterize and provide functional understanding of regulatory interactions underlying phenotypes in a cell-type specific manner by the reconstruction of tripartite networks which describe transcription factor-enhancer-gene associations. GRaNIE is a method that infers gene regulatory networks including enhancers from bulk transcriptomics (RNA-seq) and chromatin accessibility (ATAC-seq) data. Bulk data presents constraints that mask cell-type- and cell-state-specific gene expression and chromatin accessibility data as the counts represent an average of the activity of a population of cells. Hence, to overcome these limitations the approach requires its adjustment at the single-cell resolution to take into account inter-cell variation and consequently truly capture specific cell-type and cell-state regulatory differences; Therefore, a proper benchmark at the single-cell layer of the approach remains to be accomplished. 

We propose a way of benchmarking GRaNIE at the single-cell resolution through a validation-based approach that assesses the network's ability to predict real regulatory interactions from cell-type-specific promoter capture Hi-C, eQTL and ChIP-Seq data. Furthermore, we report further validation through a network biological comprehensive analysis alongside an evaluation of the eGRNs predictive power by means of GRaNPA, a machine learning framework that assesses the network's capability of  predicting cell-type-specific differential gene expression data. The results show GRaNIE, applied to single-cell hiPSC-neuron timecourse datasets, is capable of inferring neuronal differentiation and development relevant regulatory associations as well as predicting cell-type specific iPSC against neuron DGE at the single-cell level. 



### Scripts Descriptions

write here a description for each script

(maybe include figure of GRN)

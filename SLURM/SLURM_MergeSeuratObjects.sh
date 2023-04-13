#!/bin/bash
#SBATCH --job-name MERGE_SEURAT_OBJECTS
#SBATCH --account deuner
#SBATCH --mem 10000
#SBATCH --nodes 1
#SBATCH --ntasks 4
#SBATCH --partition htc
#SBATCH --time 04:00:00
#SBATCH --output MERGING_output.txt
#SBATCH --error MERGING_error.txt
 
echo "I, ${USER}, am running on host ${HOSTNAME}"

echo "STARTING JOB"
 
# Load R module
module load R

 
# Executing an R script
Rscript /g/scb/zaugg/deuner/SCENIC+/code/data_extraction.R
 
sleep 600

echo "FINISHED JOB"

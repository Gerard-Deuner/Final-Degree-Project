#!/bin/bash
#SBATCH --job-name GRANIE_BATCH_MODE_COMBINED_SPEARMAN
#SBATCH -A zaugg
#SBATCH --nodes 1
#SBATCH --ntasks 16
#SBATCH --mem 128GB
#SBATCH --time 24:00:00
#SBATCH --output output_combined_spearman.txt
#SBATCH --error error_combined_spearman.txt
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=gerard.deuner@embl.de 

echo "I, ${USER}, am running on host ${HOSTNAME}"
 
echo "STARTING JOB"
 
# Load R module
module load R
module load R-bundle-Bioconductor 
 
# Executing an R script
Rscript --vanilla /g/scb/zaugg/deuner/GRaNIE/code/GRaNIE_batch_mode.R combined spearman

sleep 600

echo "FINISHED JOB"

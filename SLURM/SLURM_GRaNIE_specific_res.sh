#!/bin/bash
#SBATCH --job-name GRANIE_16
#SBATCH -A zaugg
#SBATCH --nodes 1
#SBATCH --ntasks 16
#SBATCH --mem 128GB
#SBATCH --time 24:00:00
#SBATCH --output output_granie_comb_ps_16
#SBATCH --error error_granie_comb_ps_16
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=gerard.deuner@embl.de 

echo "I, ${USER}, am running on host ${HOSTNAME}"
 
echo "STARTING JOB"
 
# Load R module
module load R
module load R-bundle-Bioconductor/3.15-foss-2021b-R-4.2.0
 
# Executing an R script
Rscript --vanilla /g/scb/zaugg/deuner/GRaNIE/code/GRaNIE_specific_resolution.R combined pearson 16

sleep 600

echo "FINISHED JOB"

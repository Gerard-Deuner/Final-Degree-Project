#!/bin/bash
#SBATCH --job-name PANDO_TIMECOURSE
#SBATCH -A zaugg
#SBATCH --nodes 1
#SBATCH --ntasks 16
#SBATCH --mem 128GB
#SBATCH --time 24:00:00
#SBATCH --output output_pando_timecourse.txt
#SBATCH --error error_pando_timecourse.txt
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=gerard.deuner@embl.de 

echo "I, ${USER}, am running on host ${HOSTNAME}"
 
echo "STARTING JOB"
 
# Load R module
module load R
module load R-bundle-Bioconductor 
 
# Executing an R script
Rscript --vanilla /g/scb/zaugg/deuner/Pando/code/Pando.R

sleep 600

echo "FINISHED JOB"

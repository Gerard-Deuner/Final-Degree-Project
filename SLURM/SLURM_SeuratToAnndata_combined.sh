#!/bin/bash
#SBATCH --job-name SEURATTOANNDATA_COMBINED
#SBATCH -A zaugg
#SBATCH --nodes 1
#SBATCH --ntasks 16
#SBATCH --mem 128GB
#SBATCH --time 16:00:00
#SBATCH --output output_SeuratToAnndata_combined.txt
#SBATCH --error error_SeuratToAnndata_combined.txt
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=gerard.deuner@embl.de 

echo "I, ${USER}, am running on host ${HOSTNAME}"
 
echo "STARTING JOB"
 
# Load R module
module load R
module load R-bundle-Bioconductor 
 
# Executing an R script
Rscript --vanilla /g/scb/zaugg/deuner/SCENIC+/code/SeuratToAnndata.R combined

sleep 600

echo "FINISHED JOB"

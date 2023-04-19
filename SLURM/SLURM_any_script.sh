#!/bin/bash
#SBATCH --job-name MYJOB
#SBATCH -A zaugg
#SBATCH --nodes 1
#SBATCH --ntasks 16
#SBATCH --mem 128GB
#SBATCH --time 24:00:00
#SBATCH --output output_myjob.txt
#SBATCH --error error_myjob.txt
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=gerard.deuner@embl.de 

echo "I, ${USER}, am running on host ${HOSTNAME}"
 
echo "STARTING JOB"
 
# Load R module
module load R
module load R-bundle-Bioconductor/3.15-foss-2021b-R-4.2.0
 
 
# Executing an R script
# Introduce here any script and its arguments
Rscript --vanilla /g/scb/zaugg/deuner/valdata/code/ChiPSeq_setup.R combined

sleep 600

echo "FINISHED JOB"

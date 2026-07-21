#!/bin/bash
#SBATCH --mem=60G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-22

i=${SLURM_ARRAY_TASK_ID}

echo "Running chr${i}"

module purge
module load R/4.4.1

Rscript 06.coloc_eQTL_caQTL.R "CD4T_chromatin" "All_CD4T_cells" ${i}

echo "Done"

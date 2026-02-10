#!/bin/bash
#SBATCH --job-name=005_PreparecovsFile_percell_ncells         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --partition=dev  
#SBATCH --mem=16G
#SBATCH --time=07:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=ALL       # send email when job ends or fails
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH -o logs/005_PreparecovsFile_percell_ncells-%j.out

module load R/4.3.1

#file=/gpfs/commons/home/mmatos/cd4_QTL_analysis/02_Gene_expression/cell_lvl2.txt

#export celltype=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $file | awk '{print $1}' | xargs)

Rscript 005_PreparecovsFile_ncells.R #$celltype


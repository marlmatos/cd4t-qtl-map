#!/bin/bash
#SBATCH --job-name=prep_obj         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --partition=dev  
#SBATCH --mem=200G
#SBATCH --time=07:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=ALL       # send email when job ends or fails
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH -o logs/001_prep_obj-%j.out

module load R/4.3.1

R CMD BATCH 001_preparing_seurat_obj_for_QTL.R logs/001_preparing_seurat_obj_for_QTL.Rout
 
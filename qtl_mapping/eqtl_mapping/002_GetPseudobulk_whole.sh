#!/bin/bash
#SBATCH --job-name=002_pseudobulk_obj         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --partition=bigmem  
#SBATCH --mem=800G
#SBATCH --time=07:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=ALL       # send email when job ends or fails
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH -o logs/002_pseudobulk_obj-%j.out

module load R/4.3.1

R CMD BATCH 002_GetPseudobulk_whole.R logs/002_GetPseudobulk_whole.Rout

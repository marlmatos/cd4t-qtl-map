#!/bin/bash
#SBATCH --job-name=liftover        # create a short name for your job
#SBATCH -p cpu 
#SBATCH -t 02:00:00
#SBATCH --cpus-per-task=16                   
#SBATCH --mem=200G                     # Job memory reques
#SBATCH --mail-type=END,FAIL             # send email when job begins
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH --output=liftover-%j.out

module load R
Rscript chrombonet_var_liftover_hg38_to_hg19.R  
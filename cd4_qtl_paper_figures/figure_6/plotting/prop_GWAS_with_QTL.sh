#!/bin/bash
#SBATCH --job-name=coloc_sum        # create a short name for your job
#SBATCH -p cpu 
#SBATCH -t 1-00:00:00
#SBATCH --cpus-per-task=4                   
#SBATCH --mem=60G                     # Job memory reques
#SBATCH --mail-type=END,FAIL             # send email when job begins
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH --output=coloc_sum-%j.out

module load R
Rscript prop_GWAS_with_QTL.R
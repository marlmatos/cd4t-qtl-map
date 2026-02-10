#!/bin/bash
#SBATCH --job-name=finemapping_summary        # create a short name for your job
#SBATCH -p cpu 
#SBATCH -t 1-00:00:00
#SBATCH --cpus-per-task=10                   
#SBATCH --mem=100G                        # Job memory reques
#SBATCH --mail-type=END,FAIL             # send email when job begins
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH --output=finemapping_summary-%j.out

module load R
Rscript creating_summary_gwas_chrombpnet_coloc.R
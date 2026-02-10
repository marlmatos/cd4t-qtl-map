#!/bin/bash
#SBATCH --job-name=bpnet_gwas        # create a short name for your job
#SBATCH -p cpu 
#SBATCH -t 02:00:00
#SBATCH --cpus-per-task=4                   
#SBATCH --mem=48G                     # Job memory reques
#SBATCH --mail-type=END,FAIL             # send email when job begins
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH --output=bpnet_gwas-%j.out

module load R
Rscript prop_cbpnet_gwas.R  
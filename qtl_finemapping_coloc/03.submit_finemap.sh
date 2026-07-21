#!/bin/bash
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --job-name=finemap
#SBATCH -o ./Reports/finemap-%j.out
#SBATCH --time=7-00:00:00     # 7 days

module purge
module load R/4.4.1

# Finemap regions
Rscript 03.susie_finemap_eQTLs.R $1 $2 $3 $4 $5

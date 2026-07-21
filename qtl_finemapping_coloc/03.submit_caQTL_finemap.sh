#!/bin/bash
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=finemap
#SBATCH -o ./Reports/finemap-%j.out

module purge
module load R/4.4.1

# Finemap regions
Rscript 03.susie_finemap_caQTLs.R $1 $2 $3 $4 $5

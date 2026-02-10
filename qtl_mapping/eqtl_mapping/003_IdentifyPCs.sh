#!/bin/bash
#SBATCH --job-name=003_IdentifyPCs         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --partition=dev  
#SBATCH --mem=100G
#SBATCH --time=07:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=ALL       # send email when job ends or fails
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH -o logs/003_IdentifyPCs-%j.out

module load R/4.3.1

R CMD BATCH 003_IdentifyPCs.R logs/003_IdentifyPCs.Rout

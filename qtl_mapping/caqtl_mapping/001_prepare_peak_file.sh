#!/bin/bash
#SBATCH --job-name=001_prepare_peak_file         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --partition=dev  
#SBATCH --mem=30G
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=ALL       # send email when job ends or fails
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH -o logs/001_prepare_peak_file-%j.out

module load R/4.3.1

R CMD BATCH 001_prepare_peak_file.R logs/001_prepare_peak_file.Rout

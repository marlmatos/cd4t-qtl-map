#!/bin/bash
#SBATCH --job-name=04_ProcessExpression         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --partition=pe2  
#SBATCH --mem=30G
#SBATCH --time=07:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=ALL       # send email when job ends or fails
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH -o logs/04_ProcessExpression-%j.out

module load R/4.3.1


Rscript 004_ProcessExpression.R 


#!/bin/bash
#SBATCH --job-name=clust        # create a short name for your job
#SBATCH -p cpu 
#SBATCH -t 1-00:00:00
#SBATCH --cpus-per-task=4                   
#SBATCH --mem=30G                     # Job memory reques
#SBATCH --mail-type=END,FAIL             # send email when job begins
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH --output=clust-%j.out

module load R
Rscript 01_run_peak_coaccessibility_analysiS.R
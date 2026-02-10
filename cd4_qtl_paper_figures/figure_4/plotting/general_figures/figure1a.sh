#!/bin/bash
#SBATCH --job-name=plot        # create a short name for your job
#SBATCH -p cpu 
#SBATCH -t 1-00:00:00
#SBATCH --cpus-per-task=4                   
#SBATCH --mem=100G                     # Job memory reques
#SBATCH --mail-type=END,FAIL             # send email when job begins
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH --output=plot-%j.out

module load R
Rscript figure_upset_plot.R
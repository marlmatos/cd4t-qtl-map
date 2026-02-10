#!/bin/bash
#SBATCH --job-name=plotting        # create a short name for your job
#SBATCH -p bigmem 
#SBATCH -t 00:30:00
#SBATCH --cpus-per-task=4                   
#SBATCH --mem=400G                     # Job memory reques
#SBATCH --mail-type=END,FAIL             # send email when job begins
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH --output=plotting-%j.out

module load R
Rscript figure1_main_plot_rel_annotations.R
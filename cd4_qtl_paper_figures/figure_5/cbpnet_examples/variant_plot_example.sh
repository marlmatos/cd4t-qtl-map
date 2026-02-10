#!/bin/bash
#SBATCH --job-name=plot_bpnet        # create a short name for your job
#SBATCH -p cpu 
#SBATCH -t 00:15:00
#SBATCH --cpus-per-task=4                   
#SBATCH --mem=100G                     # Job memory reques
#SBATCH --mail-type=END,FAIL             # send email when job begins
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH --output=logs/plot_bpnet-%j.out

# keep Python from seeing ~/.local site-packages
# before conda activate
export ADDR2LINE="${ADDR2LINE:-addr2line}"
# 
# 
# module purge
# source ~/.bashrc
# conda activate renviron_ne1
module load R
export RETICULATE_PYTHON="${CONDA_PREFIX}/bin/python"
export PYTHONNOUSERSITE=1


Rscript qtl_variant_plot_example_TRBV28.R 
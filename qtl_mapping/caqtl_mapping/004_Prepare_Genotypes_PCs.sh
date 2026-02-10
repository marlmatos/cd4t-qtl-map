#!/bin/bash
#SBATCH --job-name=run_nb        # Job name
#SBATCH -p pe2
#SBATCH -t 7-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G                     # Job memory reques
#SBATCH --mail-type=END,FAIL             # send email when job begins
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH --time=24:00:00             # Wall time limit (days-hrs:min:sec)
#SBATCH --output=logs/004_Prepare_Genotypes_PCs.out    # Path to the standard output and error files relative to the working

#This script is to run jupyter noteboks in back ground on hpc
module load miniconda3/3.22.0
source activate genomic_tools #preinstall jupyter in the conda environment in which your job will be conducted


NOTEBOOK_TO_RUN=004_Prepare_Genotypes_PCs.ipynb
OUTPUT="~/cd4_caQTL_analysis/variant_to_peak_QTL/run_101424_cpm_tmm_maf5_fdr5_50kb/results/004_genotypes/004_Prepare_Genotypes_PCs-ran.ipynb" 

jupyter nbconvert --to notebook --execute ${NOTEBOOK_TO_RUN} --output $OUTPUT

#you can output in different format like html, pdf ,etc.
#See https://nbconvert.readthedocs.io/en/latest/usage.html

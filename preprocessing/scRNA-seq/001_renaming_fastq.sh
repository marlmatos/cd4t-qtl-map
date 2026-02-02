#!/bin/bash -l
#SBATCH --job-name=bcftools     # Job name
#SBATCH --mail-type=ALL	          # Mail events (NONE, BEGIN, END, FAIL, ALL
#SBATCH --mail-user=mmmatos@nygenome.org   
#SBATCH -p pe2 
#SBATCH -t 7-00:00:00
#SBATCH --cpus-per-task=4                   
#SBATCH --mem=32G                     # Job memory request
#SBATCH --output=001_renaming_fastqs-%j.out

source ~./bashrc
python 001_renaming_fastqs.py

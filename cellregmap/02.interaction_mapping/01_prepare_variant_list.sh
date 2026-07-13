#!/bin/bash

#SBATCH --job-name=CRM5                  # Job name
#SBATCH --nodes=1                # node count               
#SBATCH --partition=pe2                # Partition Name
#SBATCH --mem=48G
#SBATCH --time=24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --output=logs/mofa_%A.log               # Standard output and error log
#SBATCH --error=logs/mofa_%A.err
#SBATCH --mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mmatos@nygenome.org      # Where to send mail 

module load miniconda3/3.22.0

source activate renviron

python 01_prepare_variant_list.py 

#!/bin/bash

#SBATCH --job-name=caQTL_tensor                 # Job name
#SBATCH --partition=gpu                      # Partition Name
#SBATCH --mem=8G
#SBATCH --gres=gpu
#SBATCH --output=logs/1mb/001_caQTL_narrowpeaks_tensor-1mb_%A_%a.log               # Standard output and error log
#SBATCH --error=logs/1mb/001_caQTL_narrowpeaks_tensor-1mb_%A_%a.err
#SBATCH --mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --array=1-22
#SBATCH --mail-user=mmatos@nygenome.org      # Where to send mail 

module purge
module load TensorFlow R

#check for packages
which R
which python
which tensorqtl


# Define the chromosome number based on SLURM_ARRAY_TASK_ID
chromosome=$SLURM_ARRAY_TASK_ID


python 001_caQTL_narrowpeaks_tensor-1mb.py "$chromosome"

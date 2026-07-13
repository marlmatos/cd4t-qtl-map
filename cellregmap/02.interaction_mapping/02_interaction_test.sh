#!/bin/bash
#SBATCH --job-name=cellregmap                  # Job name
#SBATCH --nodes=1                # node count               
#SBATCH --partition=cpu                # Partition Name
#SBATCH --mem=120G
#SBATCH --cpus-per-task=4  # Request 4 CPU cores for the task
#SBATCH --time=2-00:00:00          # total run time limit (HH:MM:SS)
#SBATCH --output=logs/cellregmap_%A_%a.log               # Standard output and error log
#SBATCH --array=0-6142
#SBATCH --error=logs/cellregmap_%A_%a.err

source activate cellregmap


## SLURM ARRAY ##
echo "$SLURM_ARRAY_TASK_ID"

python 02_interaction_test_102024_one_var_pergene.py ${SLURM_ARRAY_TASK_ID}

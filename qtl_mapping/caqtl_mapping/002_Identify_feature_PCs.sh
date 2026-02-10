#!/bin/bash
#SBATCH --job-name=002_Identify_feature_PCs         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --partition=dev  
#SBATCH --mem=100G
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=ALL       # send email when job ends or fails
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH -o logs/002_Identify_feature_PCs-%j.out

# List of filtered tables
#filtered_tables=("filtered_qsmooth_norm_cpm.txt")

# Select the file based on the SLURM_ARRAY_TASK_ID
#selected_table=${filtered_tables[$SLURM_ARRAY_TASK_ID]}


selected_table="filtered_qsmooth_norm_cpm.txt"
# Load necessary modules
module load R/4.3.1

# Run the R script with the array index passed as an argument
Rscript 002_Identify_feature_PCs.R $selected_table
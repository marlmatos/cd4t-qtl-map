#!/bin/bash
#SBATCH --job-name=003_Process_input_peaks         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --partition=dev  
#SBATCH --mem=100G
#SBATCH --time=00:30:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=ALL       # send email when job ends or fails
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH -o logs/003_Process_input_peaks-%j.out


##-------for array only--------
# List of filtered tables
#filtered_tables=("filtered_raw.txt" "filtered_tmm_cpm.txt" "filtered_vst.txt" "stringent_filtered_tmm_cpm.txt" "stringent_filtered_vst.txt")

# Select the file based on the SLURM_ARRAY_TASK_ID
#selected_table=${filtered_tables[$SLURM_ARRAY_TASK_ID]}
###-----------------------------



selected_table="filtered_qsmooth_norm_cpm.txt"

# Load necessary modules
module load R/4.3.1


Rscript 003_Process_input_peaks.R $selected_table

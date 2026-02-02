#!/bin/bash
#SBATCH --job-name=averg_contr_bw           # Job name
#SBATCH --nodes=1                          # node count
#SBATCH --partition=cpu                   # Partition Name
#SBATCH --mem=100G                          # Memory per node (increased for merge)
#SBATCH --cpus-per-task=8                  # CPU cores for the task
#SBATCH --time=02:00:00                    # Time limit (HH:MM:SS)
#SBATCH --output=logs/averg_contr_bw_%j.log # Standard output log
#SBATCH --error=logs/averg_contr_bw_%j.err  # Standard error log
#SBATCH --mail-type=END,FAIL              # Mail events
#SBATCH --mail-user=mmatos@nygenome.org   # Where to send mail

# Load necessary modules
source activate renviron_ne1

python /gchm/cd4_chrombpnet/scripts/9.average_contribution_scores_h5.py
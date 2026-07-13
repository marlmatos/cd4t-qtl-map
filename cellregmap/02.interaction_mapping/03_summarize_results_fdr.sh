#!/bin/bash
#SBATCH --job-name=summary_cellregmap                  # Job name
#SBATCH --nodes=1                # node count               
#SBATCH --partition=cpu                # Partition Name
#SBATCH --mem=30G
#SBATCH --cpus-per-task=4  # Request 4 CPU cores for the task
#SBATCH --time=24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --output=logs/summary_%A.log               # Standard output and error log
#SBATCH --error=logs/summary_%A.err

#SBATCH --mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)

#SBATCH --mail-user=mmatos@nygenome.org      # Where to send mail 


source activate renviron_ne1


python 03_summarize_results_fdr.py \
  --res 34 \
  --results_path ~/cd4_CellRegMap/002_interaction_analysis/results/results_042426/many_lead_snps_per_gene_res34_5factors/p_values \
  --out_results ~/cd4_CellRegMap/002_interaction_analysis/results/results_042426/many_lead_snps_per_gene_res34_5factors/fdr_corrected \
  --outfile 2026_04_26_many_lead_snps_per_gene_res34_5factors_summary_bonferroni_fdr01.csv \
  --correction bonferroni_fdr \
  --alpha 0.01


  python 03_summarize_results_fdr.py \
  --res 34 \
  --results_path ~/cd4_CellRegMap/002_interaction_analysis/results/results_042426/one_lead_snp_per_gene_res34_5factors/p_values \
  --out_results ~/cd4_CellRegMap/002_interaction_analysis/results/results_042426/one_lead_snp_per_gene_res34_5factors/fdr_corrected \
  --outfile 2026_04_26_one_lead_snp_per_gene_res34_5factors_summary_fdr01.csv \
  --correction fdr \
  --alpha 0.01

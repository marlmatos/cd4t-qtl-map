#!/bin/bash
#SBATCH -p cpu
#SBATCH -J concat
#SBATCH -c 1
#SBATCH --mem=64G
#SBATCH -o logs/concat.out
#SBATCH -e logs/concat.err

set -euo pipefail

in_dir="/gchm/cd4_qtl_paper_figures/figure_3/data/mediation_res_filter_cis_only"

# make sure directories exist
mkdir -p logs
mkdir -p "${in_dir}/concatenated"

# ---------- caQTL ----------
prefix="caqtl_CD4T_combined_sig_shared_inpeak_findr_results_mediation_raw_annotated_"
out="${in_dir}/concatenated/caqtl_mediation_cis_annotated_allchr.tsv"

head -n 1 "${in_dir}/${prefix}chr1.tsv" > "$out"
for c in $(seq 1 22); do
  tail -n +2 "${in_dir}/${prefix}chr${c}.tsv" >> "$out"
done

# ---------- eQTL ----------
prefix="eqtl_all_CD4T_cells_gene_expression_mediation_raw_annotated_"
out="${in_dir}/concatenated/eqtl_mediation_cis_annotated_allchr.tsv"

head -n 1 "${in_dir}/${prefix}chr1.tsv" > "$out"
for c in $(seq 1 22); do
  tail -n +2 "${in_dir}/${prefix}chr${c}.tsv" >> "$out"
done
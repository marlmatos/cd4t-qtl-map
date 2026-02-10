#!/bin/bash
#SBATCH -p cpu
#SBATCH -J cis_filter
#SBATCH -c 1
#SBATCH --mem=64G
#SBATCH -o logs/cis_filter.out
#SBATCH -e logs/cis_filter.err

python /gchm/cd4_qtl_paper_figures/figure_3/analysis/mediation_dec25/med_res_filtering/mediation_results_filter.py --atype 'gene'
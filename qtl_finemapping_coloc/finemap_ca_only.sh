#!/bin/bash

numbers=($(seq 1 22))

for i in "${numbers[@]}"; do

echo "Finemapping chr${i}"

NAME="CD4T_chromatin"
DIR=$(pwd)
EQTL="/gchm/cd4_caQTL_analysis/variant_to_peak_QTL/run_012625_qc_aware_qsmooth_CPM_MAF5_FDR5_1MB/results/006_caQTLs/filtered_qsmooth_norm_cpm_1mb/cd4_qsmooth_cpm_chromatin_narrowpeaks.cis_qtl_pairs.chr${i}.csv"
INDEP="/gchm/cd4_caQTL_analysis/variant_to_peak_QTL/run_012625_qc_aware_qsmooth_CPM_MAF5_FDR5_1MB/results/006_caQTLs/filtered_qsmooth_norm_cpm_1mb/cd4_qsmooth_cpm_chromatin_narrowpeaks.independent_cis_QTL_pairs.chr${i}.csv"

# Finemap regions
sbatch 03.submit_caQTL_finemap.sh \
        "${EQTL}" \
        "regions/${NAME}.chr${i}.regions.txt" \
        "Ashkenazi_LD_matrices/" \
        "${NAME}" \
        367

done

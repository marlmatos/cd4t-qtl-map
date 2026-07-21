#!/bin/bash

numbers=($(seq 1 22))

for i in "${numbers[@]}"; do

echo "Finemapping chr${i}"


NAME="All_CD4T_cells"
EQTL="/gchm/cd4_QTL_analysis/03_Run_cisQTL_perchr/analysis/cis_eQTLs/all_CD4T_cells_MAF5/allcells.cis_qtl_pairs.chr${i}.parquet"
DIR=$(pwd)

sbatch 03.submit_finemap.sh \
    ${EQTL} \
    "regions/${NAME}.chr${i}.regions.txt" \
    "Ashkenazi_LD_matrices/" \
    ${NAME} \
    367

done

echo "Finemap process complete!"

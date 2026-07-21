#!/bin/bash

module purge
module load R/4.4.1

EQTL_NAME="CD4T_chromatin"

# obtain list of GWAS file names
find /gchm_lab/collab_gwas/finemapping_autoimmune/output/nathan_completed/ \
  -type d -name "*preprocessed" -exec basename {} \; > preprocessed_folder_names.txt

while read line; do
  echo "Submitting job for $line"

  # Create a job script for each line
  sbatch --cpus-per-task=8 --mem=40G --nodes=1 -o ./Reports/coloc-$line-%j.out --wrap="Rscript 05.Run_coloc_immune_caQTL_GWAS.R $line ${EQTL_NAME}"

done < preprocessed_folder_names.txt

echo "All jobs submitted!"

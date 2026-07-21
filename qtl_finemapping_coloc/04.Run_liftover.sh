#!/bin/bash
#SBATCH --mem=60G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --array=3

i=${SLURM_ARRAY_TASK_ID}

module purge
module load bedtools
module load R/4.4.1

mkdir -p bedfiles
mkdir -p SuSiE_finemap_lbf_GRC37

DIR="/gchm_lab/collab/marlis_pj/coloc"
QTLS=("All_CD4T_cells" "CD4T_chromatin")

for QTL in "${QTLS[@]}"; do

  echo "lifting over ${QTL}"

  mkdir -p "SuSiE_finemap_lbf_GRC37/${QTL}"

  # Use R to combine files because bash was acting weird
  Rscript 04.liftover_finemapped_QTLs.R ${i} ${QTL}

  echo "Done!"

done

echo "Done!"

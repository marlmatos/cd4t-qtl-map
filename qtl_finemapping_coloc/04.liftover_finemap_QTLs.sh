#!/bin/bash
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-22

i=${SLURM_ARRAY_TASK_ID}

module purge
module load bedtools
module load R/4.4.1

mkdir -p bedfiles
mkdir -p SuSiE_finemap_lbf_GRC37

DIR="/gchm_lab/collab/marlis_pj/coloc"
QTLS=("All_CD4T_cells")

for QTL in "${QTLS[@]}"; do

  echo ${QTL}

  mkdir -p SuSiE_finemap_lbf_GRC37/${QTL}

  # Format lbf file as a bed file
  awk -F'\t' '
  BEGIN {
    OFS = "\t"
    # Define required columns
    required_cols = "chr pos PIP variant_id.x lbf_1 lbf_2 lbf_3 lbf_4 lbf_5 lbf_6 lbf_7 lbf_8 lbf_9 lbf_10 beta slope_se pval_nominal af a0 a1 variant_id.y start_distance gene"
    split(required_cols, req_arr, " ")
  }
  
  NR == 1 {
    # Map column names to field indices
    for (i = 1; i <= NF; i++) {
      header[$i] = i
    }
  
    # Check all required columns exist
    for (j in req_arr) {
      col = req_arr[j]
      if (!(col in header)) {
        print "ERROR: Missing required column \"" col "\" in header." > "/dev/stderr"
        exit 1
      }
    }
  
    # Print output header
    print "chr", "start", "end", "PIP", "variant_id.x", "lbf_1", "lbf_2", "lbf_3", "lbf_4", "lbf_5",
          "lbf_6", "lbf_7", "lbf_8", "lbf_9", "lbf_10", "beta", "slope_se", "pval_nominal", "af",
          "a0", "a1", "variant_id.y", "start_distance", "gene"
    next
  }
  
  NF >= 26 {
    chr = "chr"$(header["chr"])
    pos = $(header["pos"])
    start = pos - 1
    end = pos
  
    print chr, start, end, $(header["PIP"]), $(header["variant_id.x"]), $(header["lbf_1"]), $(header["lbf_2"]), $(header["lbf_3"]),
          $(header["lbf_4"]), $(header["lbf_5"]), $(header["lbf_6"]), $(header["lbf_7"]), $(header["lbf_8"]), $(header["lbf_9"]),
          $(header["lbf_10"]), $(header["beta"]), $(header["slope_se"]), $(header["pval_nominal"]), $(header["af"]),
          $(header["a0"]), $(header["a1"]), $(header["variant_id.y"]), $(header["start_distance"]), $(header["gene"])
  }
   ' SuSiE_finemap_lbf/${QTL}/${QTL}_chr${i}_lbfs.txt > SuSiE_finemap_lbf/${QTL}/${QTL}_chr${i}_lbfs.bed  

  # First step is to obtain rsid from position ranges
  # Run bedtools intersect on the data (excluding header) and append result
  tail -n +2 SuSiE_finemap_lbf/${QTL}/${QTL}_chr${i}_lbfs.bed | \
  bedtools intersect -wa -wb -a stdin -b $home/genome/dbsnp151_GRCh38p7/bed_chr_${i}.bed.gz >> bedfiles/grc38_rsids_chr${i}.txt

  # Use R to combine files because bash was acting weird
  Rscript 04.liftover_finemapped_QTLs.R ${i}

done

echo "Done!"

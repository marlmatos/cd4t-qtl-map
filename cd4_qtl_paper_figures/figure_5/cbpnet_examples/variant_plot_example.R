.libPaths(c("/gpfs/commons/home/mmatos/R/x86_64-pc-linux-gnu-library/4.4","/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"))
setwd("~/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/cbpnet_examples")
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(rtracklayer)
  library(BSgenome)
  library(data.table)
  library(BSgenome.Hsapiens.UCSC.hg38)
})
source("~/cd4_qtl_paper_figures/utils/track_plots_helpers.R")
source("~/cd4_qtl_paper_figures/utils/cbpnet_extract_allelic_contrib_plus_atrr_bigwigs.R")
args = commandArgs(trailingOnly=TRUE)

# ----------------------- RESOURCES -------------------------------------------
atrr         <- "/gpfs/commons/home/mmatos/cd4_chrombpnet/chrombpnet_model_b7/variant_contribution_scores/averaged_cd4_tcells_AJ_common_variants.shap.counts.h5"
data_h5      <- "/gpfs/commons/home/mmatos/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_variant_prediction_scores.h5"
tsv          <- "/gpfs/commons/home/mmatos/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_scores/average_cd4_tcells_AJ_common_variants.mean.variant_scores.tsv"
chrom_sizes  <- "/gpfs/commons/home/mmatos/resources/genome/hg38.chrom.sizes"
motifs_path  <- "~/cd4_chrombpnet/chrombpnet_model_b7/tfmodisco_motifs_hocomoco_jaspar_cisbp/motif_instances/motifs_with_tf_annotated_filt_c90_40seq_unique_variant_hits_granges.rds"
bw_path <- "/gpfs/commons/groups/lappalainen_lab/mmatos/cd4_aging_project/data/ATAC-seq/atac_preprocessing/nextflow_tmp/a4/4c07c68cd77cef99a8b320e8e66238/cd4_atac.bigWig"


# # ----------------------- CLI ARGUMENT -----------------------------------------
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) < 1) {
#   stop('Usage: Rscript plot_variant_tracks.R "<variant_id>"\n',
#        'Allowed formats: "chr6_113306741_G_A" or "6:113306741[b38]G,A"\n', call. = FALSE)
# }

# Ensure exporter initialized once
if (!exists(".init_variant_bw_exporter")) {
  stop(".init_variant_bw_exporter() not found. Did you source the Python helper?")
}
.init_variant_bw_exporter()

# ----------------------- Variant PARSER ----------------------------------------
# Only:
#   - "chr6_113306741_G_A"
#   - "6:113306741[b38]G,A"
parse_variant_id <- function(variant_id) {
  x <- trimws(variant_id)
  mA <- stringr::str_match(x, "^chr([0-9XYMT]+)_([0-9]+)_[ACGT]_[ACGT]$")
  if (!is.na(mA[1,1])) return(list(chr = paste0("chr", mA[1,2]), pos = as.integer(mA[1,3])))
  mB <- stringr::str_match(x, "^([0-9XYMT]+):([0-9]+)\\[b\\d+\\][ACGT],[ACGT]$")
  if (!is.na(mB[1,1])) return(list(chr = paste0("chr", mB[1,2]), pos = as.integer(mB[1,3])))
  stop('Unrecognized variant_id format. Allowed only: "chr6_113306741_G_A" or "6:113306741[b38]G,A".')
}

# #get variants to plot
# vars_motifs<-read_tsv("~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap/motif_variant_overlap_all.tsv")
# 
# vars_motifs_top16 <- vars_motifs %>%
#   filter(match_confidence == "strong") %>%
#   group_by(TF_annotation) %>%
#   slice_max(hit_coefficient_global, n = 16, with_ties = FALSE) %>%
#   ungroup()
# 
# variants <- vars_motifs_top16$variant_id
# writeLines(variants, "/gpfs/commons/home/mmatos/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/cbpnet_examples/variant_motifHits_plot.txt")

#var_file=fread("/gpfs/commons/home/mmatos/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/cbpnet_examples/variant_motifHits_plot.txt", header = FALSE)
# variant <-"chr5_50489690_A_G"
# variants <- as.character(var_file$V1)


 
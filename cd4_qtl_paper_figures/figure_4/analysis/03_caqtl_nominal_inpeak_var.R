################################
## Quantifying concoordance between caQTLs and chromBPNet predictions for variants within peaks
## Author: Marliette Matos
## Date: 10/13/2025
################################

##Read and filter caQTL nominal caQTL for variants wihtin peaks

.libPaths(c("/gchm/R/x86_64-pc-linux-gnu-library/4.4", "/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"))
setwd("~/cd4_qtl_paper_figures/figure_4")

library(dplyr)
options(bitmapType="cairo")

#####################
#read chromBPNet results
#####################

cBPNet_var<-read.delim2("/gchm/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_scores/average_cd4_tcells_AJ_common_variants.mean.variant_scores.tsv")
cBPNet_variants <-cBPNet_var$variant_id

#####################
#read caQTL Tensor nominal results
#####################

dir <- "/gchm/cd4_caQTL_analysis/variant_to_peak_QTL/run_012625_qc_aware_qsmooth_CPM_MAF5_FDR5_1MB/results/006_caQTLs/filtered_qsmooth_norm_cpm_1mb"

# Grab every file that matches “…chr{N}_credible_sets.txt”
nom_files <- list.files(
  dir,
  pattern = "^cd4_qsmooth_cpm_chromatin_narrowpeaks\\.cis_qtl_pairs\\.chr(\\d|1\\d|2[0-2])\\.csv$",
  full.names = TRUE
)

if (length(nom_files) == 0)
  stop("No files found in ", finemapping_dir)

message("Reading ", length(nom_files), " nominal association files…")

# Read them all, preserving column classes, then row-bind
res <- do.call(
  rbind,
  lapply(nom_files, function(f) {
    read.delim2(f, stringsAsFactors = FALSE)
  })
)

res <- res %>%
  mutate(var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])", 
                       "chr\\1_\\2_\\3_\\4", 
                       variant_id)) %>% filter(var_id %in% cBPNet_variants)

write.csv(res, "data/CD4T_caQTL_nomical_assoc_inpeak_vars.csv", row.names = FALSE)

#####################
# filter only for variant associations where the variant in within the peak for which it's associated to
#####################

#res <- read.delim2("data/CD4T_caQTL_nomical_assoc_inpeak_vars.csv", sep = ",")
res <- res %>%
  mutate(
    var_chr = str_split_fixed(var_id, "_", 4)[,1],
    var_pos = as.numeric(str_split_fixed(var_id, "_", 4)[,2]),

  )
# Make GRanges for variants from `res`
variant_gr <- GRanges(
  seqnames = res$var_chr,
  ranges = IRanges(start = res$var_pos, end = res$var_pos),
  peakID = res$phenotype_id  # for matching to the associated peak
)

#get peak coordinates
peak_path<-"/gchm/ATAC-seq_analysis/diff_accesibility_ana/results/peak_counts/RAW_cd4_atac_peakcounts_ranges_scrna_union.csv"
ranges.table <- read.csv2(peak_path, sep = ",") %>%
  filter(Chr %in% paste0("chr", 1:22))
rownames(ranges.table) <- ranges.table$X
gr <- GRanges(seqnames = ranges.table$Chr,
              ranges = IRanges(ranges.table$Start, ranges.table$End),
              mcols = data.frame(peakID = rownames(ranges.table)))
names(gr) <- gr$mcols.peakID
gr


# Find overlaps where variant is inside a peak
hits <- findOverlaps(variant_gr, gr)

# Filter for cases where the phenotype_id matches the peak
matching_hits <- hits[
  variant_gr$peakID[queryHits(hits)] == gr$mcols.peakID[subjectHits(hits)]
]

# Get the indices of matching rows in res
matching_indices <- queryHits(matching_hits)

# Filter the dataframe
res_filtered <- res[matching_indices, ] #PLEASE NOTE THAT THE NUMBER OF VARIANTS DROPS FROM 160K TO 49K 
#(I PRESUME IT IS BECAUSE THE DIFFERENCE IN PEAK CALLING BETWEEN METHODS?)

ggplot(res_filtered, aes(as.numeric(slope), -log(as.numeric(pval_nominal)))) + geom_point()

ggplot(cBPNet_var, aes(as.numeric(logfc.mean), -log(as.numeric(abs_logfc_x_jsd_x_active_allele_quantile.mean.pval)))) + geom_point()


#####plot
# Merge 
cBPNet_var_unique <- cBPNet_var %>%
  distinct(variant_id, .keep_all = TRUE)

var_eff <- left_join(res_filtered, cBPNet_var_unique, by = c("var_id" = "variant_id"))
var_eff$logfc.mean <- as.numeric(var_eff$logfc.mean)
var_eff$slope <- as.numeric(var_eff$slope)
var_eff$slope_se <- as.numeric(var_eff$slope_se)
var_eff$pval_nominal <- as.numeric(var_eff$pval_nominal)

var_eff$significant <- "not_significant"
var_eff <- var_eff %>%
  mutate(significant = "Not Significant") %>%
  mutate(
    significant = ifelse(
      abs_logfc_x_jsd_x_active_allele_quantile.mean.pval < 0.05 & pval_nominal >= 0.05,
      "chromBPnet", significant
    ),
    significant = ifelse(
      abs_logfc_x_jsd_x_active_allele_quantile.mean.pval >= 0.05 & pval_nominal < 0.05,
      "caQTL", significant
    ),
    significant = ifelse(
      abs_logfc_x_jsd_x_active_allele_quantile.mean.pval < 0.05 & pval_nominal < 0.05,
      "Both", significant
    ),
    significant = factor(significant, levels = signif_levels)
  )

table(var_eff$significant)
ggplot(var_eff, aes(y = logfc.mean, x = slope, color = significant)) + geom_point(aes( alpha=0.5), size=0.5)

write.csv(var_eff, "data/CD4T_caQTL_nomical_assoc_inpeak_chrombpnet.csv", row.names = FALSE)

# Set consistent significance levels and custom colors
signif_levels <- c("Both", "chromBPnet", "caQTL", "Not Significant")



################################
## Script to comapre concoordance 
## caQTL vs ChromBpnet (all variants)
## Author: Marliette Matos
## Date: 10/13/2025
################################s

## Read and filter caQTL nominal caQTL for variants within peaks
.libPaths(c("/gchm/R/x86_64-pc-linux-gnu-library/4.4",
            "/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"))
setwd("/gchm/cd4_qtl_paper_figures/figure_4")

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(glue)
  library(GenomicRanges)
  library(IRanges)
  library(readr)
  library(fastmatch)
  library(data.table)
  
})

options(bitmapType = "cairo")

dir.create("data", showWarnings = FALSE)
dir.create("qc",   showWarnings = FALSE)

# ---- QC helpers -------------------------------------------------------------
.qc_log <- list()
log_qc <- function(label, df, keys = NULL) {
  stopifnot(!missing(df))
  out <- tibble(
    step = label,
    n_rows = nrow(df),
    !!!setNames(lapply(keys %||% character(), function(k) {
      if (k %in% names(df)) dplyr::n_distinct(df[[k]], na.rm = TRUE) else NA_integer_
    }), paste0("n_distinct_", keys %||% character()))
  )
  .qc_log[[length(.qc_log) + 1]] <<- out
  # Pretty console message
  msg <- glue("{label}: rows={out$n_rows}, ",
              paste0(names(out)[grepl('^n_distinct_', names(out))], "=",
                     unlist(out[grepl('^n_distinct_', names(out))]), collapse=", "))
  message(msg)
  invisible(out)
}

warn_if_increase <- function(curr, prev, label_curr, label_prev) {
  if (!is.null(prev) && curr > prev) {
    warning(glue("QC: Count increased from '{label_prev}' ({prev}) to '{label_curr}' ({curr}). ",
                 "If unexpected, double-check filtering and joins."))
  }
}

# Finalize QC table to disk
finalize_qc <- function(outfile = file.path("qc", glue("qc_counts_{format(Sys.time(), '%Y%m%d_%H%M%S')}.csv"))) {
  qc_tbl <- dplyr::bind_rows(.qc_log)
  readr::write_csv(qc_tbl, outfile)
  message(glue("QC table written: {outfile}"))
  return(qc_tbl)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Read ChromBPNet results -----------------------------------------------
cBPNet_var <- read.delim2(
  "/gchm/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_scores/average_cd4_tcells_AJ_common_variants.mean.variant_scores.tsv",
  stringsAsFactors = FALSE
)
log_qc("chromBPNet: raw", cBPNet_var, keys = c("variant_id"))

cBPNet_variants <- cBPNet_var$variant_id
n_cBPNet <- length(cBPNet_variants)
message(glue("chromBPNet variants: n={n_cBPNet}, unique={dplyr::n_distinct(cBPNet_variants)}"))

# ---- Read caQTL Tensor nominal results -------------------------------------
dir_nom <- "/gchm/cd4_caQTL_analysis/variant_to_peak_QTL/run_012625_qc_aware_qsmooth_CPM_MAF5_FDR5_1MB/results/006_caQTLs/filtered_qsmooth_norm_cpm_1mb/"

# Optional speed-up for `%in%` on big vectors
use_fastmatch <- requireNamespace("fastmatch", quietly = TRUE)
f_in <- if (use_fastmatch) fastmatch::fmatch else match

# Files (allow .csv/.tsv and optional .gz)
nom_files <- list.files(
  dir_nom,
  pattern = "^cd4_qsmooth_cpm_chromatin_narrowpeaks\\.cis_qtl_pairs\\.chr(\\d|1\\d|2[0-2])\\.(csv|tsv)(\\.gz)?$",
  full.names = TRUE
)
if (length(nom_files) == 0) stop("No files found in ", dir_nom)
message("Reading ", length(nom_files), " nominal association files…")

# Universe (unique + keep in memory once)
cBP_set <- unique(cBPNet_variants)

# Precompile the regex and replacement
re_var <- "([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])"
re_rep <- "chr\\1_\\2_\\3_\\4"  # -> chr10_12345_G_A

dt_list <- lapply(seq_along(nom_files), function(i) {
  f <- nom_files[i]
  dt0 <- fread(f, showProgress = FALSE)
  dt0[, file := basename(f)]
  
  # Per-file raw QC (before variant/filter)
  # guard missing columns to avoid log errors
  if (all(c("variant_id", "phenotype_id") %in% names(dt0))) {
    log_qc("nominal: raw (per-file)", dt0, keys = c("variant_id", "phenotype_id"))
  } else {
    log_qc("nominal: raw (per-file, limited keys)", dt0, keys = intersect(c("variant_id", "phenotype_id"), names(dt0)))
  }
  
  # Standardize and filter in-file
  dt0[, var_id := sub(re_var, re_rep, variant_id)]
  keep <- !is.na(dt0$var_id)
  # membership via fastmatch if available
  keep <- keep & !is.na(f_in(dt0$var_id, cBP_set))
  # optional: also require phenotype_id present
  if ("phenotype_id" %in% names(dt0)) keep <- keep & !is.na(dt0$phenotype_id)
  
  dt <- dt0[keep]
  
  # Per-file QC after var_id
  if (nrow(dt)) {
    log_qc("nominal: add var_id (per-file)", dt, keys = c("var_id", "phenotype_id"))
  }
  
  dt
})

# Bind once
res <- rbindlist(dt_list, use.names = TRUE, fill = TRUE, idcol = "chunk")

# Final QC
log_qc("nominal: filtered to ChromBPNet variants", res, keys = c("var_id", "phenotype_id"))

# Write intermediate
#write.csv(res, "data/CD4T_caQTL_nomical_assoc_inpeak_vars.csv", row.names = FALSE)

# --- Ensure var coords exist & are sane
res_var <- res %>%
  mutate(
    var_chr = str_split_fixed(var_id, "_", 4)[, 1],
    var_pos = suppressWarnings(as.numeric(str_split_fixed(var_id, "_", 4)[, 2]))
  ) %>%
  mutate(
    var_chr = ifelse(grepl("^chr", var_chr), var_chr, paste0("chr", gsub("^chr", "", var_chr)))  ) %>%
  filter(!is.na(var_chr), !is.na(var_pos)) %>%
  mutate(row_id = row_number())

# --- GRanges for variants (point positions)
variant_gr <- GRanges(
  seqnames = res_var$var_chr,
  ranges   = IRanges(start = res_var$var_pos, end = res_var$var_pos)
)
mcols(variant_gr)$row_id       <- res_var$row_id
mcols(variant_gr)$var_id       <- res_var$var_id
mcols(variant_gr)$phenotype_id <- res_var$phenotype_id

# Read peak coordinates (union peak set)
peak_path <- "/gchm/ATAC-seq_analysis/diff_accesibility_ana/results/peak_counts/RAW_cd4_atac_peakcounts_ranges_scrna_union.csv"
ranges.table <- read.csv2(peak_path, sep = ",") %>%
  filter(Chr %in% paste0("chr", 1:22))
rownames(ranges.table) <- ranges.table$X

gr <- GRanges(
  seqnames = ranges.table$Chr,
  ranges   = IRanges(ranges.table$Start, ranges.table$End),
  peakID   = rownames(ranges.table)
)
names(gr) <- gr$peakID
message(glue("GRanges(peaks): n={length(gr)}"))

# --- All peaks within 1kb (including overlaps)
# maxgap counts bp between ranges; with points this is exactly the distance threshold.
hits_1kb <- findOverlaps(
  variant_gr, gr,
  maxgap = 1000,
  ignore.strand = TRUE
)

# Distance for QC (optional but helpful)
dist_1kb <- distance(variant_gr[queryHits(hits_1kb)], gr[subjectHits(hits_1kb)])

nearby_df <- tibble(
  row_id         = mcols(variant_gr)$row_id[queryHits(hits_1kb)],
  nearby_peakID  = mcols(gr)$peakID[subjectHits(hits_1kb)],
  distance_bp    = as.integer(dist_1kb)  # 0 == inside peak
)

# --- Keep only rows whose phenotype_id is among the peaks within 1kb of their variant
res_1kb <- res_var %>%
  left_join(nearby_df, by = "row_id") %>%
  filter(phenotype_id == nearby_peakID) %>%
  select(-nearby_peakID)

# If the same (var_id, phenotype_id) appears multiple times via replicated rows, dedup:
res_1kb <- res_1kb %>%
  group_by(var_id, phenotype_id) %>%
  slice_min(distance_bp, with_ties = FALSE) %>%
  ungroup()

# # --- Quick QC
# n_any_peak <- sum(distance_bp == 0, na.rm = TRUE)
# message(sprintf("Pairs within 1kb: %d / %d (%.1f%%)",
#                 nrow(res_1kb), nrow(res),
#                 100 * nrow(res_1kb) / nrow(res)))
# 
# # Optional: keep ALL peaks within 1kb per variant (expanded table)
# # (Useful if you want to explore alternative peak assignments.)
# res_all_peaks_1kb <- res %>%
#   left_join(nearby_df, by = "row_id") %>%
#   filter(!is.na(distance_bp))
# 




caqtl_eff <- res_1kb %>% group_by(var_id) %>% summarise(mean.beta=mean(slope),
                                                    mean.slope_se=mean(slope_se),
                                                    mean.pvalue=mean(pval_nominal)
                                                    )


# ---- Plots (optional quick looks) ------------------------------------------
# ggplot(res_filtered, aes(as.numeric(slope), -log(as.numeric(pval_nominal)))) + geom_point()
# ggplot(cBPNet_var, aes(as.numeric(logfc.mean), -log(as.numeric(abs_logfc_x_jsd_x_active_allele_quantile.mean.pval)))) + geom_point()

# ---- Merge with ChromBPNet per-variant effects ------------------------------
cBPNet_var_unique <- cBPNet_var %>% distinct(variant_id, .keep_all = TRUE)
log_qc("chromBPNet: unique variants", cBPNet_var_unique, keys = c("variant_id"))

var_eff <- left_join(caqtl_eff, cBPNet_var_unique, by = c("var_id" = "variant_id"))
log_qc("merge: caQTL x chromBPNet", var_eff, keys = c("var_id", "phenotype_id"))

# Numeric coercions
num_cols <- c("logfc.mean", "mean.beta", "mean.slope_se", "mean.pvalue")
for (cc in intersect(num_cols, names(var_eff))) {
  var_eff[[cc]] <- suppressWarnings(as.numeric(var_eff[[cc]]))
}

# Define significance levels BEFORE using them
signif_levels <- c("Both", "chromBPnet", "caQTL", "Not Significant")

var_eff <- var_eff %>%
  mutate(
    significant = "Not Significant",
    significant = ifelse(abs_logfc_x_jsd_x_active_allele_quantile.mean.pval < 0.05 & mean.pvalue >= 0.05, "chromBPnet", significant),
    significant = ifelse(abs_logfc_x_jsd_x_active_allele_quantile.mean.pval >= 0.05 & mean.pvalue < 0.05, "caQTL", significant),
    significant = ifelse(abs_logfc_x_jsd_x_active_allele_quantile.mean.pval < 0.05 & mean.pvalue < 0.05, "Both", significant),
    significant = factor(significant, levels = signif_levels)
  )
log_qc("post-classify significance", var_eff, keys = c("var_id", "significant"))

# Table printout for quick QC
sig_tab <- table(var_eff$significant, useNA = "ifany")
message("Significance breakdown:\n", capture.output(print(sig_tab)) |> paste(collapse = "\n"))

# Example plot (you can tune later)
var_eff$mean.pvalue<- as.numeric(var_eff$mean.pvalue)
var_eff$mean.beta<- as.numeric(var_eff$mean.beta)
var_eff$logfc.mean<- as.numeric(var_eff$logfc.mean)
var_eff$significant <- factor(var_eff$significant, 
                                levels = c("Both", "caQTL", "chromBPnet", "Not Significant"))

head(var_eff)
library(tidyr)

# Define significance level order and custom colors
signif_levels <- c("Both", "caQTL", "chromBPnet", "Not Significant")
custom_colors <- c(
  "Both" = "#e4b57a",
  "caQTL" = "#9CCB86",
  "chromBPnet" = "#0076C0",
  "Not Significant" = "gray90"
)

# Make sure `significant` is a factor with correct levels
var_eff <- var_eff %>%
  mutate(significant = factor(significant, levels = signif_levels))

# Count number of variants per group
category_counts <- var_eff %>%
  count(significant) %>%
  complete(significant = signif_levels, fill = list(n = 0))

# Create legend labels with counts
color_labels <- category_counts %>%
  arrange(match(significant, signif_levels)) %>%
  mutate(label = paste0(significant, " (n=", n, ")")) %>%
  pull(label)


# subset
both_df <- dplyr::filter(var_eff, significant == "Both")

# R² for Both
r_both  <- cor(both_df$mean.beta, both_df$logfc.mean, use = "complete.obs")
r2_both <- r_both^2
n_both  <- sum(complete.cases(both_df$mean.beta, both_df$logfc.mean))
r2_lab  <- sprintf("Both: R^2 = %.3f (n = %d)", r2_both, n_both)

set.seed(42)
var_eff <- var_eff[sample(nrow(var_eff)), ]

ca_vareff_plot_r2 <- ggplot( var_eff, aes(x = mean.beta, y = logfc.mean, fill = significant)) +
  ggrastr::geom_point_rast(shape = 21, size = 2, alpha = 0.6, stroke = 0.1, color = "gray10") +
  geom_smooth(data = both_df, aes(x = mean.beta, y = logfc.mean), method = "lm", se = FALSE,
              color = custom_colors["Both"], size = 0.8) +
  # ← diagonal reference
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  annotate("label", x = Inf, y = Inf, label = r2_lab, hjust = 1.02, vjust = 1.2, size = 3,
           label.size = 0, fill = scales::alpha(custom_colors["Both"], 0.15)) +
  scale_fill_manual(values = custom_colors, labels = color_labels) +
  guides(fill = guide_legend(override.aes = list(size = 3), nrow = 2)) +
  labs(x = "caQTL (β)", y = "ChromBPNet log(aFC)", fill = "Significance:",
       title = "Variant Allelic Effect Correlation") +
  theme_light() +
  theme(legend.title = element_text(face = "bold"),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ca_vareff_plot_r2

write.csv(var_eff, "data/CD4T_caQTL_nomical_assoc_1kb_radious_chrombpnet_V2.csv", row.names = FALSE)


ggsave("~/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/general_figures/cbpnet_vs_caqtl_beta_1kb.pdf", ca_vareff_plot_r2, width =5, height = 5)

# ---- Save QC table ----------------------------------------------------------
qc_table <- finalize_qc()
print(qc_table)

length((cBPNet_var_unique %>% filter(abs_logfc_x_jsd_x_active_allele_quantile.mean.pval < 0.05))$variant_id)


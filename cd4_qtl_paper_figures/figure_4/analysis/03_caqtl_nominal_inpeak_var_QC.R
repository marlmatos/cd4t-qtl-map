################################
## Script to comapre concoordance 
## caQTL vs ChromBpnet (only variants in their caPeaks)
## Author: Marliette Matos
## Date: 10/13/2025
################################

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

# Remove ASCII & Unicode whitespace (incl. non-breaking space), keep case
clean_ws <- function(x) {
  x <- gsub("\u00A0", " ", x, fixed = TRUE)  # NBSP -> space
  x <- trimws(x)                              # trim ends
  gsub("\\s+", "", x)                         # remove all remaining whitespace
}

# ---- Read ChromBPNet results -----------------------------------------------
cBPNet_var <- read.delim2(
  "/gchm/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_scores/average_cd4_tcells_AJ_common_variants.mean.variant_scores.tsv",
  stringsAsFactors = FALSE
)
log_qc("chromBPNet: raw", cBPNet_var, keys = c("variant_id"))

cBPNet_variants <- cBPNet_var$variant_id
cBP_set <- unique(na.omit(clean_ws(cBPNet_variants)))

n_cBPNet <- length(cBP_set)
message(glue("chromBPNet variants: n={n_cBPNet}, unique={dplyr::n_distinct(cBPNet_variants)}"))

# ---- Read caQTL Tensor nominal results -------------------------------------
dir_nom <- "/gchm/cd4_caQTL_analysis/variant_to_peak_QTL/run_012625_qc_aware_qsmooth_CPM_MAF5_FDR5_1MB/results/006_caQTLs/filtered_qsmooth_norm_cpm_1mb"

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
log_qc("nominal: filtered to ChromBPNet variants", res, keys = c("var_id", "phenotype_id"))

# Write intermediate
write.csv(res, "data/CD4T_caQTL_nomical_assoc_inpeak_vars.csv", row.names = FALSE)

# ---- Keep only associations where the variant lies within its associated peak
# Parse var coords
res <- res %>%
  mutate(
    var_chr = str_split_fixed(var_id, "_", 4)[, 1],
    var_pos = suppressWarnings(as.numeric(str_split_fixed(var_id, "_", 4)[, 2]))
  )
log_qc("nominal: add var coords", res, keys = c("var_id", "phenotype_id"))

# Make GRanges for variants (tag with phenotype_id as peakID we expect to match)
variant_gr <- GRanges(
  seqnames = res$var_chr,
  ranges   = IRanges(start = res$var_pos, end = res$var_pos),
  peakID   = res$phenotype_id
)
message(glue("GRanges(variants): n={length(variant_gr)}"))

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

# Overlaps: variant inside any peak
hits <- findOverlaps(variant_gr, gr)
message(glue("Overlaps (variant in ANY peak): n_hits={length(hits)}"))

# Keep only those where the nominal association's phenotype_id matches the peak ID
matching_hits <- hits[ variant_gr$peakID[queryHits(hits)] == gr$peakID[subjectHits(hits)] ]
message(glue("Overlaps (variant in ITS associated peak): n_hits={length(matching_hits)}"))

matching_indices <- queryHits(matching_hits)
res_filtered <- res[matching_indices, ]
log_qc("nominal: variant within its associated peak", res_filtered, keys = c("var_id", "phenotype_id"))



# ---- Plots (optional quick looks) ------------------------------------------
# ggplot(res_filtered, aes(as.numeric(slope), -log(as.numeric(pval_nominal)))) + geom_point()
# ggplot(cBPNet_var, aes(as.numeric(logfc.mean), -log(as.numeric(abs_logfc_x_jsd_x_active_allele_quantile.mean.pval)))) + geom_point()

# ---- Merge with ChromBPNet per-variant effects ------------------------------
cBPNet_var_unique <- cBPNet_var %>% distinct(variant_id, .keep_all = TRUE)
log_qc("chromBPNet: unique variants", cBPNet_var_unique, keys = c("variant_id"))

var_eff <- left_join(res_filtered, cBPNet_var_unique, by = c("var_id" = "variant_id"))
log_qc("merge: caQTL x chromBPNet", var_eff, keys = c("var_id", "phenotype_id"))

# Numeric coercions
num_cols <- c("logfc.mean", "slope", "slope_se", "pval_nominal")
for (cc in intersect(num_cols, names(var_eff))) {
  var_eff[[cc]] <- suppressWarnings(as.numeric(var_eff[[cc]]))
}

# Define significance levels BEFORE using them
signif_levels <- c("Both", "chromBPnet", "caQTL", "Not Significant")

var_eff <- var_eff %>%
  mutate(
    significant = "Not Significant",
    significant = ifelse(abs_logfc_x_jsd_x_active_allele_quantile.mean.pval < 0.05 & pval_nominal >= 0.05, "chromBPnet", significant),
    significant = ifelse(abs_logfc_x_jsd_x_active_allele_quantile.mean.pval >= 0.05 & pval_nominal < 0.05, "caQTL", significant),
    significant = ifelse(abs_logfc_x_jsd_x_active_allele_quantile.mean.pval < 0.05 & pval_nominal < 0.05, "Both", significant),
    significant = factor(significant, levels = signif_levels)
  )
log_qc("post-classify significance", var_eff, keys = c("var_id", "significant"))

# Table printout for quick QC
sig_tab <- table(var_eff$significant, useNA = "ifany")
message("Significance breakdown:\n", capture.output(print(sig_tab)) |> paste(collapse = "\n"))

# ggplot(var_eff, aes(y = logfc.mean, x = slope, color = significant)) +
#   geom_point(alpha = 0.5, size = 0.5)

write.csv(var_eff, "data/CD4T_caQTL_nomical_assoc_inpeak_chrombpnet_V2.csv", row.names = FALSE)

# ---- Save QC table ----------------------------------------------------------
qc_table <- finalize_qc()
print(qc_table)



.libPaths(c("/gchm/R/x86_64-pc-linux-gnu-library/4.4", "/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"))
library(dplyr)
library(ggplot2)
library(GenomicRanges)
options(bitmapType = "cairo")
library(BiocParallel)

register(MulticoreParam(workers = 4))  # Use all CPUs allocated

### Load and format input ######################################################
finemapping_dir <- "/gcgls/marlis_pj/coloc/SuSiE_finemap_credible_sets/CD4T_chromatin/"

# Grab every file that matches “…chr{N}_credible_sets.txt”
cred_files <- list.files(
  finemapping_dir,
  pattern = "^CD4T_chromatin_chr[0-9XYM]+_credible_sets\\.txt$",
  full.names = TRUE
)

if (length(cred_files) == 0)
  stop("No credible-set files found in ", finemapping_dir)

message("Reading ", length(cred_files), " credible-set files…")

# Read them all, preserving column classes, then row-bind
res <- do.call(
  rbind,
  lapply(cred_files, function(f) {
    read.delim2(f, stringsAsFactors = FALSE)
  })
)


# Standardise columns ----------------------------------------------------------
res$variant_pos <- as.numeric(sub(".*:(\\d+)\\[.*", "\\1", res$variant_id))
res$chr         <- paste0("chr", res$chr)   # add UCSC-style prefix

peak_coords <- read.delim2("/gchm/cd4_caQTL_analysis/variant_to_peak_QTL/run_012625_qc_aware_qsmooth_CPM_MAF5_FDR5_1MB/results/003_inputs/filtered_qsmooth_norm_cpm/filtered_qsmooth_norm_cpm_cd4_atac_processed_peaks.bed")
peak_coordsGR <- GRanges(peak_coords$X.chr, IRanges(peak_coords$start, peak_coords$end), peak_name = peak_coords$peak_name)
names(peak_coordsGR) <- peak_coordsGR$peak_name

# Subset only to peaks tested
rownames(peak_coords) <- peak_coords$peak_name
peak_counts <- peak_coords %>% select(-c(X.chr, start, end, peak_name))
peak_coordsGR <- peak_coordsGR[names(peak_coordsGR) %in% rownames(peak_counts)]

# caQTL peaks from credible sets
relevant_peaks <- unique(res$peak)
caQTL_peaksGR <- peak_coordsGR[names(peak_coordsGR) %in% relevant_peaks]

message("Starting correlation analysis on ", length(caQTL_peaksGR), " lead peaks...")

### Parallelized correlation computation ###
correlation_results_list <- bplapply(seq_along(caQTL_peaksGR), function(i) {
  lead_peak <- caQTL_peaksGR[i]
  lead_name <- names(lead_peak)
  
  if (!lead_name %in% rownames(peak_counts)) return(NULL)
  
  windowGR <- GRanges(
    seqnames = seqnames(lead_peak),
    ranges = IRanges(start = start(lead_peak) - 2e6, end = end(lead_peak) + 2e6)
  )
  
  hits <- findOverlaps(windowGR, peak_coordsGR)
  nearby_names <- setdiff(names(peak_coordsGR[subjectHits(hits)]), lead_name)
  
  if (length(nearby_names) == 0) return(NULL)
  
  res <- lapply(nearby_names, function(cor_name) {
    if (!cor_name %in% rownames(peak_counts)) return(NULL)
    
    x <- as.numeric(peak_counts[lead_name, ])
    y <- as.numeric(peak_counts[cor_name, ])
    
    ct <- suppressWarnings(cor.test(x, y, method = "pearson"))
    
    data.frame(
      leadPeak = lead_name,
      corPeak = cor_name,
      cor = ct$estimate,
      pval = ct$p.value
    )
  })
  
  do.call(rbind, res)
})

# Combine results
correlation_results <- do.call(rbind, correlation_results_list)
correlation_results$FDR <- p.adjust(correlation_results$pval, method = "fdr")

### Save output ###
write.csv(
  correlation_results,
  "/gchm/cd4_qtl_paper_figures/figure_1/data/caPeak_coacc_correlation_results_allchrs.csv",
  row.names = FALSE
)

message("Finished correlation analysis.")

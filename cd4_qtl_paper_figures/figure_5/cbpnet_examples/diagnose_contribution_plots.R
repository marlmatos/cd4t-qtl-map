# Diagnostic script to check what's actually in your contribution BigWigs

diagnose_contribution_bigwigs <- function(atrr_ref_bw, atrr_alt_bw, region_zoom) {
  
  cat("\n================================================================================\n")
  cat("CONTRIBUTION BIGWIG DIAGNOSTIC\n")
  cat("================================================================================\n\n")
  
  # Get the region
  region_gr <- str_to_gr(region_zoom)
  
  cat("Region:", as.character(region_gr), "\n")
  cat("Width:", width(region_gr), "bp\n\n")
  
  # Extract scores
  ref_scores <- as.numeric(mcols(atrr_ref_bw)$score)
  alt_scores <- as.numeric(mcols(atrr_alt_bw)$score)
  
  cat("=== REF Allele Scores ===\n")
  cat("Total positions:", length(ref_scores), "\n")
  cat("Non-zero positions:", sum(ref_scores != 0, na.rm = TRUE), "\n")
  cat("Finite positions:", sum(is.finite(ref_scores)), "\n")
  cat("Range: [", min(ref_scores, na.rm = TRUE), ",", max(ref_scores, na.rm = TRUE), "]\n")
  cat("Mean:", mean(ref_scores, na.rm = TRUE), "\n")
  cat("Median:", median(ref_scores, na.rm = TRUE), "\n")
  cat("Sum of absolute values:", sum(abs(ref_scores), na.rm = TRUE), "\n")
  
  cat("\nFirst 20 REF scores:\n")
  print(head(ref_scores, 20))
  
  cat("\n=== ALT Allele Scores ===\n")
  cat("Total positions:", length(alt_scores), "\n")
  cat("Non-zero positions:", sum(alt_scores != 0, na.rm = TRUE), "\n")
  cat("Finite positions:", sum(is.finite(alt_scores)), "\n")
  cat("Range: [", min(alt_scores, na.rm = TRUE), ",", max(alt_scores, na.rm = TRUE), "]\n")
  cat("Mean:", mean(alt_scores, na.rm = TRUE), "\n")
  cat("Median:", median(alt_scores, na.rm = TRUE), "\n")
  cat("Sum of absolute values:", sum(abs(alt_scores), na.rm = TRUE), "\n")
  
  cat("\nFirst 20 ALT scores:\n")
  print(head(alt_scores, 20))
  
  # Check overlap with region
  cat("\n=== Overlap with region ===\n")
  ref_overlap <- atrr_ref_bw[atrr_ref_bw %over% region_gr]
  alt_overlap <- atrr_alt_bw[atrr_alt_bw %over% region_gr]
  
  cat("REF entries overlapping region:", length(ref_overlap), "\n")
  cat("ALT entries overlapping region:", length(alt_overlap), "\n")
  
  if (length(ref_overlap) > 0) {
    ref_overlap_scores <- as.numeric(mcols(ref_overlap)$score)
    cat("\nREF scores in region:\n")
    cat("  Non-zero:", sum(ref_overlap_scores != 0, na.rm = TRUE), "\n")
    cat("  Range: [", min(ref_overlap_scores, na.rm = TRUE), ",", 
        max(ref_overlap_scores, na.rm = TRUE), "]\n")
    cat("  Sum:", sum(ref_overlap_scores, na.rm = TRUE), "\n")
    cat("  First 10 values:", paste(head(ref_overlap_scores, 10), collapse = ", "), "\n")
  }
  
  if (length(alt_overlap) > 0) {
    alt_overlap_scores <- as.numeric(mcols(alt_overlap)$score)
    cat("\nALT scores in region:\n")
    cat("  Non-zero:", sum(alt_overlap_scores != 0, na.rm = TRUE), "\n")
    cat("  Range: [", min(alt_overlap_scores, na.rm = TRUE), ",", 
        max(alt_overlap_scores, na.rm = TRUE), "]\n")
    cat("  Sum:", sum(alt_overlap_scores, na.rm = TRUE), "\n")
    cat("  First 10 values:", paste(head(alt_overlap_scores, 10), collapse = ", "), "\n")
  }
  
  # Check your ylim calculation
  cat("\n=== Your ylim calculation ===\n")
  all_scores <- c(ref_scores, alt_scores)
  all_scores <- all_scores[is.finite(all_scores)]
  
  ymax <- max(all_scores) * 1.05
  ymin <- min(all_scores) * 1.05
  
  cat("ymin:", ymin, "\n")
  cat("ymax:", ymax, "\n")
  cat("Range:", ymax - ymin, "\n")
  
  if (abs(ymax) < 1e-6 && abs(ymin) < 1e-6) {
    cat("\n⚠️  WARNING: Both ymin and ymax are near zero!\n")
    cat("   Your plots will be empty or nearly empty.\n")
  }
  
  if (ymax - ymin < 1e-10) {
    cat("\n⚠️  WARNING: Y-axis range is extremely small!\n")
    cat("   The plot will appear flat.\n")
  }
  
  # Test what ggseqlogo would see
  cat("\n=== Testing ggseqlogo input ===\n")
  if (length(ref_overlap) > 0) {
    region_start <- start(region_gr)
    L <- width(region_gr)
    
    scores_per_base <- numeric(L)
    pos_idx <- as.integer(start(ref_overlap)) - region_start + 1L
    keep <- !is.na(pos_idx) & pos_idx >= 1L & pos_idx <= L & is.finite(ref_overlap_scores)
    
    if (any(keep)) {
      pos_idx_keep <- pos_idx[keep]
      vals_keep <- ref_overlap_scores[keep]
      
      # Apply ylim clipping if provided
      vals_keep <- pmin(vals_keep, ymax)
      vals_keep <- pmax(vals_keep, ymin)
      
      agg <- tapply(vals_keep, pos_idx_keep, max, na.rm = TRUE)
      scores_per_base[as.integer(names(agg))] <- as.numeric(agg)
      
      cat("Scores after processing:\n")
      cat("  Non-zero positions:", sum(scores_per_base != 0), "/", L, "\n")
      cat("  Range: [", min(scores_per_base), ",", max(scores_per_base), "]\n")
      cat("  Max absolute value:", max(abs(scores_per_base)), "\n")
      
      if (max(abs(scores_per_base)) < 1e-6) {
        cat("\n⚠️  PROBLEM: Maximum absolute value is too small for ggseqlogo!\n")
        cat("   ggseqlogo will render letters with near-zero height.\n")
      }
    } else {
      cat("⚠️  No valid positions after filtering!\n")
    }
  }
  
  cat("\n================================================================================\n\n")
}

# Quick fix function that scales values if they're too small
fix_contribution_plot <- function(bw, region, genome,
                                  track_label = "Contributions",
                                  facet_label = NULL,
                                  ylim = NULL,
                                  range_label = NULL,
                                  scale_factor = NULL) {  # NEW: auto-scale if needed
  
  region_gr <- str_to_gr(region)
  L <- width(region_gr)
  
  contrib_filt <- bw[bw %over% region_gr]
  start0 <- start(region_gr)
  scores_per_base <- numeric(L)
  
  if (length(contrib_filt) > 0) {
    pos_idx <- as.integer(start(contrib_filt)) - start0 + 1L
    keep <- !is.na(pos_idx) & pos_idx >= 1L & pos_idx <= L & is.finite(contrib_filt$score)
    
    if (any(keep)) {
      pos_idx <- pos_idx[keep]
      vals <- as.numeric(contrib_filt$score)[keep]
      
      # Clip to ylim if provided
      if (!is.null(ylim)) {
        vals <- pmin(vals, ylim[2])
        vals <- pmax(vals, ylim[1])
      }
      
      agg <- tapply(vals, pos_idx, max, na.rm = TRUE)
      scores_per_base[as.integer(names(agg))] <- as.numeric(agg)
    }
  }
  
  # AUTO-SCALE if values are too small
  max_val <- max(abs(scores_per_base))
  if (is.null(scale_factor)) {
    if (max_val > 0 && max_val < 0.01) {
      # Scale up to make values visible
      scale_factor <- 0.1 / max_val
      message("Auto-scaling contributions by factor ", round(scale_factor, 2), 
              " to make them visible")
      scores_per_base <- scores_per_base * scale_factor
      
      # Update ylim and range_label if they exist
      if (!is.null(ylim)) {
        ylim <- ylim * scale_factor
      }
      if (!is.null(range_label)) {
        range_label <- paste0(range_label, " (×", round(scale_factor, 1), ")")
      }
    }
  } else {
    scores_per_base <- scores_per_base * scale_factor
    if (!is.null(ylim)) {
      ylim <- ylim * scale_factor
    }
    if (!is.null(range_label)) {
      range_label <- paste0(range_label, " (×", round(scale_factor, 1), ")")
    }
  }
  
  # Get sequence
  region_seq_str <- as.character(Biostrings::getSeq(genome, region_gr))
  bases <- toupper(strsplit(region_seq_str, "", fixed = TRUE)[[1]])
  
  # Create one-hot matrix
  base_map <- c("A", "C", "G", "T")
  col_idx <- match(bases, base_map)
  mat_ohe <- matrix(0, nrow = L, ncol = 4, dimnames = list(NULL, base_map))
  valid <- !is.na(col_idx)
  if (any(valid)) {
    mat_ohe[cbind(seq_len(L)[valid], col_idx[valid])] <- scores_per_base[valid]
  }
  contribs_ohe <- t(mat_ohe)
  
  # Build plot
  plot <- ggseqlogo::ggseqlogo(contribs_ohe, method = "custom", seq_type = "dna") +
    ggplot2::geom_hline(yintercept = 0, color = "gray90") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = "Genomic Position (bp)", y = track_label) +
    ggplot2::guides(y = "none", fill = "none") +
    BPCells:::trackplot_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  
  if (!is.null(ylim)) {
    plot <- plot + ggplot2::scale_y_continuous(limits = ylim, expand = c(0, 0))
  }
  
  if (!is.null(range_label) && !is.null(ylim)) {
    plot <- plot +
      ggplot2::annotate(
        "text",
        x = 1,
        y = ylim[2],
        label = range_label,
        vjust = 1.5,
        hjust = -0.1,
        size = 11 * .8 / ggplot2::.pt
      )
  }
  
  trackplot <- BPCells:::wrap_trackplot(
    plot,
    ggplot2::unit(0.6, "null"),
    takes_sideplot = FALSE
  )
  
  if (!is.null(facet_label)) {
    trackplot <- BPCells:::set_trackplot_label(trackplot, labels = facet_label)
  }
  
  trackplot
}

cat("
================================================================================
CONTRIBUTION PLOT DIAGNOSTICS LOADED
================================================================================

To diagnose your contribution BigWigs:
--------------------------------------
diagnose_contribution_bigwigs(
  atrr_ref_bw = atrr_ref_bw,
  atrr_alt_bw = atrr_alt_bw,
  region_zoom = region_zoom
)

To create plots with auto-scaling if needed:
--------------------------------------------
Ref_contribs <- fix_contribution_plot(
  bw = atrr_ref_bw,
  region = region_zoom,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  ylim = ylim_shared,
  range_label = range_label
)

Alt_contribs <- fix_contribution_plot(
  bw = atrr_alt_bw,
  region = region_zoom,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  ylim = ylim_shared,
  range_label = range_label
)

================================================================================
")
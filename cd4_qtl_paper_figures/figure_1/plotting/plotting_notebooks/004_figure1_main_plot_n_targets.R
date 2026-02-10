### Script for statistical analysis and plotting of Figures 1 A-G
# Author: Marliette Matos
# Date: 01/17/2025

#NUMBER OF EGENES/CAPEAKS PER LOCUS ###
source("~/cd4_qtl_paper_figures/utils/color_pallete_helper.R")
options(bitmapType = "cairo")

#path to summary df with eQTL & caQTL results
summary_df<-"/gchm/cd4_qtl_paper_figures/figure_1/data/eqtl_caqtl_finemapping_coloc_all.tsv"


# ======================================================================
# ======================================================================
# 1) NUMBER OF FEATURES PER LOCUS
# ======================================================================
# ======================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(GenomicRanges)
  library(IRanges)
  library(igraph)
  library(purrr)
  library(scales)
  library(cowplot)
})

# ----------------------------
# Diagnostics plots
# ----------------------------
plot_cs_size_hist <- function(cs_counts, xmax = 100, binwidth = 5, xlab = "Number of variants per credible set") {
  cs_counts <- cs_counts %>%
    filter(!is.na(n_variants_per_cs), is.finite(n_variants_per_cs))
  
  s <- summary(cs_counts$n_variants_per_cs)
  n_cs <- nrow(cs_counts)
  
  stats_label <- paste0(
    "N CS = ", comma(n_cs), "\n",
    "Min = ", s[["Min."]], "\n",
    "Q1  = ", s[["1st Qu."]], "\n",
    "Median = ", s[["Median"]], "\n",
    "Mean   = ", round(mean(cs_counts$n_variants_per_cs), 2), "\n",
    "Q3  = ", s[["3rd Qu."]], "\n",
    "Max = ", s[["Max."]]
  )
  
  ggplot2::ggplot(cs_counts, aes(x = n_variants_per_cs)) +
    geom_histogram(binwidth = binwidth, boundary = 0, closed = "left") +
    coord_cartesian(xlim = c(1, xmax)) +
    labs(x = xlab, y = "Number of credible sets") +
    annotate(
      "label",
      x = xmax * 0.65,
      y = Inf,
      label = stats_label,
      vjust = 1.1,
      hjust = 0,
      label.size = 0.2
    ) +
    geom_vline(xintercept = unname(s[["Median"]]), linetype = "dashed") +
    geom_vline(xintercept = mean(cs_counts$n_variants_per_cs), linetype = "dotted") +
    theme_classic()
}

plot_jaccard_hist <- function(out_interval, binwidth = 0.02, title = "Jaccard overlap between CS intervals") {
  ggplot2::ggplot(out_interval, aes(jaccard)) +
    geom_histogram(binwidth = binwidth, boundary = 0, closed = "left") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    labs(x = "Jaccard index (overlap/union)", y = "Number of CS pairs", title = title) +
    theme_classic()
}

plot_stack_prop_two_groups <- function(binned, title) {
  ggplot2::ggplot(binned, aes(x = qtl_type, y = n, fill = bin)) +
    geom_col(width = 0.7, color = "grey20", linewidth = 0.2, position = "fill") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Proportion", title = title, fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right") +
    coord_flip()
}

# -------------------------------------------------------
# Core: for ONE QTL type (eQTL OR caQTL)
# - collapse
# - CS size diag
# - filter bottom q_keep CS by size
# - build CS intervals
# - compute Jaccard pairs + diag
# - cluster to loci (connected components) at threshold T
# - count targets per locus
# - bin targets per locus
# -------------------------------------------------------
analyze_one_qtl_type <- function(
    df,
    qtl_type_name,
    cs_col,
    target_col,
    chr_col,
    pos_col,
    q_keep = 0.75,
    jaccard_T = 0.3,
    breaks = c(1, 2, 11, 51, 100, Inf),
    labels = c("1", "2–10", "11–50", "51–99", "100+"),
    cs_xmax = 500,
    cs_binwidth = 5
) {
  stopifnot(length(breaks) == length(labels) + 1)
  
  dat <- df %>%
    transmute(
      cs_id  = as.character(.data[[cs_col]]),
      target = as.character(.data[[target_col]]),
      chr    = as.character(.data[[chr_col]]),
      pos    = as.integer(.data[[pos_col]])
    ) %>%
    filter(!is.na(cs_id), !is.na(target), !is.na(chr), !is.na(pos)) %>%
    distinct()
  
  # DIAG 1: CS sizes (pre-filter)
  cs_counts_all <- dat %>%
    distinct(cs_id, chr, pos) %>%
    dplyr::count(cs_id, name = "n_variants_per_cs")
  
  p_cs_size <- plot_cs_size_hist(
    cs_counts_all,
    xmax = cs_xmax,
    binwidth = cs_binwidth,
    xlab = paste0(qtl_type_name, ": variants per credible set")
  )
  
  # Filter to bottom q_keep CS by #variants (global within this QTL type)
  cs_sizes <- dat %>%
    distinct(cs_id, chr, pos) %>%
    dplyr::count(cs_id, name = "n_var")
  
  q_cut <- quantile(cs_sizes$n_var, probs = q_keep, na.rm = TRUE)
  
  dat_keep <- dat %>%
    left_join(cs_sizes, by = "cs_id") %>%
    filter(n_var <= q_cut) %>%
    dplyr::select(cs_id, target, chr, pos)
  
  # Build CS ranges
  cs_ranges <- dat_keep %>%
    group_by(cs_id, chr) %>%
    summarise(
      start = min(pos, na.rm = TRUE),
      end   = max(pos, na.rm = TRUE),
      n_var = n(),
      .groups = "drop"
    )
  
  gr <- GRanges(
    seqnames = cs_ranges$chr,
    ranges   = IRanges(start = cs_ranges$start, end = cs_ranges$end),
    cs_id    = cs_ranges$cs_id,
    n_var    = cs_ranges$n_var
  )
  
  # CS-vs-CS overlaps and Jaccard
  hits <- findOverlaps(gr, gr, ignore.strand = TRUE)
  pairs <- as.data.frame(hits) %>% filter(queryHits < subjectHits)
  
  if (nrow(pairs) == 0) {
    out_interval <- tibble(
      cs1 = character(), cs2 = character(), chr = character(),
      overlap_bp = integer(), union_bp = integer(),
      pct_overlap_of_smaller = numeric(), jaccard = numeric()
    )
    locus_map <- tibble(cs_id = cs_ranges$cs_id, locus_id = seq_along(cs_ranges$cs_id))
  } else {
    q <- gr[pairs$queryHits]
    s <- gr[pairs$subjectHits]
    
    overlap_bp <- width(pintersect(ranges(q), ranges(s)))
    len1 <- width(q)
    len2 <- width(s)
    union_bp <- len1 + len2 - overlap_bp
    jacc <- overlap_bp / union_bp
    
    out_interval <- tibble(
      cs1 = mcols(q)$cs_id,
      cs2 = mcols(s)$cs_id,
      chr = as.character(seqnames(q)),
      overlap_bp = overlap_bp,
      union_bp = union_bp,
      pct_overlap_of_smaller = overlap_bp / pmin(len1, len2),
      jaccard = jacc
    )
    
    # Cluster CS into loci by Jaccard threshold
    edges <- out_interval %>%
      filter(jaccard >= jaccard_T) %>%
      dplyr::select(cs1, cs2)
    
    g <- if (nrow(edges) == 0) make_empty_graph() else graph_from_data_frame(edges, directed = FALSE)
    mem <- components(g)$membership  # named: cs_id -> component
    
    all_cs <- unique(cs_ranges$cs_id)
    max_id <- if (length(mem)) max(mem) else 0L
    missing <- setdiff(all_cs, names(mem))
    
    locus_map <- tibble(
      cs_id    = c(names(mem), missing),
      locus_id = c(as.integer(mem), seq.int(max_id + 1L, max_id + length(missing)))
    )
  }
  
  # DIAG 2: Jaccard histogram (post-filter, what you actually used)
  p_jaccard <- plot_jaccard_hist(
    out_interval,
    binwidth = 0.02,
    title = paste0(qtl_type_name, ": Jaccard(CS interval overlap), T = ", jaccard_T)
  )
  
  # Count targets per locus
  targets_per_locus <- dat_keep %>%
    distinct(cs_id, target) %>%
    left_join(locus_map, by = "cs_id") %>%
    group_by(locus_id) %>%
    summarise(n_targets = n_distinct(target), .groups = "drop")
  
  # Bin targets per locus
  binned <- targets_per_locus %>%
    mutate(
      bin = cut(n_targets, breaks = breaks, labels = labels,
                right = FALSE, include.lowest = TRUE)
    ) %>%
    dplyr::count(bin, name = "n") %>%
    tidyr::complete(bin = factor(labels, levels = labels), fill = list(n = 0L)) %>%
    mutate(
      qtl_type = qtl_type_name,
      prop = n / sum(n),
      bin  = factor(bin, levels = rev(labels))
    )
  
  list(
    qtl_type = qtl_type_name,
    binned = binned,
    contingency_row = setNames(as.integer(binned$n), as.character(binned$bin)),
    targets_per_locus = targets_per_locus,   # <--- ADD THIS LINE
    diag_cs_size = p_cs_size,
    diag_jaccard = p_jaccard,
    out_interval = out_interval,
    locus_map = locus_map
  )
}

# -------------------------------------------------------
# Global eQTL vs caQTL comparison (single 2xB test only)
# -------------------------------------------------------
compare_two_qtl_types_global <- function(binned_eqtl, binned_caqtl, labels) {
  binned_all <- dplyr::bind_rows(binned_eqtl, binned_caqtl) %>%
    dplyr::mutate(
      bin      = factor(as.character(bin), levels = rev(labels)),
      qtl_type = factor(qtl_type, levels = unique(qtl_type))
    )
  
  tab <- binned_all %>%
    dplyr::select(qtl_type, bin, n) %>%
    tidyr::pivot_wider(names_from = bin, values_from = n, values_fill = 0) %>%
    dplyr::arrange(qtl_type)
  
  mat <- as.matrix(tab[, -1, drop = FALSE])
  rownames(mat) <- as.character(tab$qtl_type)
  
  # One global test: chi-square if expected counts ok; otherwise Fisher
  chisq0 <- suppressWarnings(stats::chisq.test(mat))
  exp_ok <- all(chisq0$expected >= 5)
  
  global <- if (exp_ok) {
    chisq0
  } else {
    # Fisher can be slow for many bins; keep defaults unless you want to tune workspace/simulate.p.value
    stats::fisher.test(mat)
  }
  
  # Optional diagnostics if chi-square was used (what bins drive the association)
  diagnostics <- if (inherits(global, "htest") && identical(global$method, chisq0$method)) {
    list(
      expected = chisq0$expected,
      std_resid = chisq0$stdres,
      contrib = (chisq0$observed - chisq0$expected)^2 / chisq0$expected
    )
  } else {
    NULL
  }
  
  list(
    binned_all   = binned_all,
    contingency  = mat,
    global_test  = global,
    diagnostics  = diagnostics,
    plot = plot_stack_prop_two_groups(
      binned_all,
      title = "Targets per locus: binned distribution (eQTL vs caQTL)"
    )
  )
}


# =======================================================
# ===================== RUN IT ==========================
# =======================================================

# ---- Load / filter your eQTL and caQTL universes ----
# You control these filters. The key is: eqtl_df0 contains only eQTL CS universe,

eqtl_only <- c("coloc", "eQTL_only")
caqtl_only <- c("coloc", "caQTL_only") 

df0 <- fread(summary_df)

eqtl_df0 <- df0 %>% filter(coloc_status %in% eqtl_only)
caqtl_df0 <- df0 %>% filter(coloc_status %in% caqtl_only)

# ---- Common binning parameters (shared so distributions are comparable) ----
brk  <- c(1, 2, 11,  Inf)
labs <- c("1", "2–10", "11+")

# ---- Run per-type analyses (Jaccard computed globally within type) ----
eqtl_res <- analyze_one_qtl_type(
  df = eqtl_df0,
  qtl_type_name = "eQTL",
  cs_col = "finemapped_cs_eqtl",
  target_col = "gene",
  chr_col = "chromosome",
  pos_col = "variant_pos.y",
  q_keep = 0.75,
  jaccard_T = 0.3,
  breaks = brk,
  labels = labs,
  cs_xmax = 100,
  cs_binwidth = 5
)

caqtl_res <- analyze_one_qtl_type(
  df = caqtl_df0,
  qtl_type_name = "caQTL",
  cs_col = "finemapped_cs_caqtl",
  target_col = "peak",
  chr_col = "chromosome",
  pos_col = "variant_pos.x",   
  q_keep = 0.75,
  jaccard_T = 0.3,
  breaks = brk,
  labels = labs,
  cs_xmax = 100,
  cs_binwidth = 5
)

# ---- Compare bin proportions: eQTL vs caQTL ----
cmp <- compare_two_qtl_types_global(eqtl_res$binned, caqtl_res$binned, labels = labs)

# ---- Print diagnostics and main results ----
s.pt1<-eqtl_res$diag_cs_size
s.pt2<-eqtl_res$diag_jaccard
s.pt3<-caqtl_res$diag_cs_size
s.pt4<-caqtl_res$diag_jaccard

print(cmp$plot)
cmp$global_test
cmp$contingency

p <- cmp$global_test$p.value
format(p, scientific = TRUE)
#"2.833199e-78"

p_labeled <- cmp$binned_all %>%
  group_by(qtl_type) %>%
  mutate(prop = n / sum(n),
         pct  = percent(prop, accuracy = 0.1)) %>%
  ungroup() %>%
  ggplot(aes(x = qtl_type, y = n, fill = bin)) +
  geom_col(position = "fill", width = 0.7, color = "grey20", linewidth = 0.2) +
  geom_text(
    aes(label = ifelse(prop > 0.02, pct, "")),  # hide tiny slivers
    position = position_fill(vjust = 0.5),
    size = 3
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Proportion of loci", fill = NULL, title = "Targets per locus (binned)") +
  coord_flip() +
  theme_minimal(base_size = 12)

print(p_labeled)


##supplements


supplements<-plot_grid(s.pt1, s.pt2, s.pt3, s.pt4)
supplements2<-plot_grid(s.pt1,  s.pt3, ncol=1)

supplements
p_labeled

ggsave("~/cd4_qtl_paper_figures/figure_1/plotting/plots_oct2025/caQTLs_characteristics_summary_n_targets_suppl_1.pdf",supplements,
       width = 8, height = 6)
ggsave("~/cd4_qtl_paper_figures/figure_1/plotting/plots_oct2025/caQTLs_characteristics_summary_n_targets_suppl_4.pdf",supplements2,
       width = 5, height = 8)
ggsave("~/cd4_qtl_paper_figures/figure_1/plotting/plots_oct2025/caQTLs_characteristics_n_targets.pdf",p_labeled,
        width = 5, height = 3)



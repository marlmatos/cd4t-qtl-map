#!/usr/bin/env Rscript

# ==============================================================================
# Motif family promoter enrichment (Fisher's exact test)
# ------------------------------------------------------------------------------
# Question:
#   Are variants with a hit in TF family F enriched in promoters vs distal,
#   relative to the universe of all tested variants?
#
# Unit of analysis: variant (each variant counted at most once per family)
# Universe: all tested variants mapped to ChromBPNet ATAC peaks (promoter/distal)
# Hits: base_df rows collapsed to variant_id x family
# Test: 2x2 Fisher per family + BH correction
# Output: results table + forest plot of top families
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
})

# -----------------------------
# Config
# -----------------------------
paths <- list(
  motif_hits     = "~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap_hit_caller/motif_hits_summary_annotations.csv",
  variant_scores = "~/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_scores/average_cd4_tcells_AJ_common_variants.mean.variant_scores.tsv",
  peaks_bed      = "~/cd4_chrombpnet/data/inputs/peaks/merged_cd4_samples_peaks_no_blacklist.sorted.bed",
  gencode_bed    = "~/resources/genome/hg38_gencode_raw/gencode.v44.primary_assembly.annotation.genes.bed"
)

top_n     <- 25
min_hits  <- 4  # minimum total hit variants per family to report
out_table <- "~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/TF_motif_family_promoter_enrich_fisher_results.tsv"
out_plot  <- "~/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/variant_motif_overlap/tf_motif_family_promoter_enrich_fisher_forestplot_top25.pdf"

# -----------------------------
# Helpers
# -----------------------------
recode_region <- function(x) {
  x <- str_to_lower(as.character(x))
  case_when(
    str_detect(x, "promot") ~ "Promoter",
    str_detect(x, "^dist")  ~ "Distal",
    TRUE ~ NA_character_
  )
}

sig_stars <- function(q) {
  case_when(
    q < 1e-4 ~ "****",
    q < 1e-3 ~ "***",
    q < 1e-2 ~ "**",
    q < 5e-2 ~ "*",
    TRUE ~ ""
  )
}

fisher_one <- function(a, b, c, d) {
  # matrix:
  #          Promoter  Distal
  # Hit         a       b
  # NoHit       c       d
  m <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  ft <- fisher.test(m)
  tibble(
    odds_ratio = unname(ft$estimate),
    conf_low   = unname(ft$conf.int[1]),
    conf_high  = unname(ft$conf.int[2]),
    p_value    = ft$p.value
  )
}

# -----------------------------
# 1) Load data
# -----------------------------
base_df <- fread(paths$motif_hits)

cbpnet_all <- fread(paths$variant_scores) %>% distinct()

# -----------------------------
# 2) Peak annotation: promoter vs distal (TSS overlap)
# -----------------------------
bp_peaks <- fread(paths$peaks_bed, header = FALSE)

bp_peaks_gr <- GRanges(
  seqnames  = bp_peaks$V1,
  ranges    = IRanges(bp_peaks$V2, bp_peaks$V3),
  peak_name = bp_peaks$V4
)

gencode <- fread(paths$gencode_bed, sep = "\t", header = FALSE, data.table = FALSE) %>%
  transmute(
    chr        = V1,
    strand     = V6,
    gene_id    = V7,
    gene       = V4,
    tss_1based = if_else(V6 == "+", V2 + 1L, V3)
  ) %>%
  distinct(chr, gene, .keep_all = TRUE)

tss_gr <- GRanges(
  seqnames = gencode$chr,
  ranges   = IRanges(start = gencode$tss_1based, end = gencode$tss_1based),
  strand   = gencode$strand
)

ov_tss <- countOverlaps(bp_peaks_gr, tss_gr, ignore.strand = TRUE)
mcols(bp_peaks_gr)$annotation <- if_else(ov_tss > 0, "promoter", "distal")

# -----------------------------
# 3) Map ALL tested variants to peak_annotation (universe)
# -----------------------------
vars_gr_all <- cbpnet_all %>%
  separate(variant_id, into = c("chr", "pos", "ref", "alt"), sep = "_", remove = FALSE) %>%
  mutate(pos = as.integer(pos)) %>%
  {
    GRanges(
      seqnames   = .$chr,
      ranges     = IRanges(start = .$pos, end = .$pos),
      strand     = "*",
      variant_id = .$variant_id
    )
  }

ov <- findOverlaps(vars_gr_all, bp_peaks_gr, ignore.strand = TRUE)

ann_vec <- rep(NA_character_, length(vars_gr_all))
ann_vec[queryHits(ov)] <- mcols(bp_peaks_gr)$annotation[subjectHits(ov)]
mcols(vars_gr_all)$peak_annotation <- ann_vec

universe <- tibble(
  variant_id = vars_gr_all$variant_id,
  region     = recode_region(mcols(vars_gr_all)$peak_annotation)
) %>%
  filter(!is.na(region)) %>%
  distinct(variant_id, region)

N_prom <- sum(universe$region == "Promoter")
N_dist <- sum(universe$region == "Distal")

# -----------------------------
# 4) Hits collapsed to variant x family (variant is the unit)
# -----------------------------
hits_vxf <- base_df %>%
  filter(!is.na(family), family != "") %>%
  distinct(variant_id, family) %>%
  inner_join(universe, by = "variant_id")  # restrict to tested variants w/ region

# counts per family x region
fam_counts <- hits_vxf %>%
  count(family, region, name = "n_hit") %>%
  pivot_wider(names_from = region, values_from = n_hit, values_fill = 0) %>%
  mutate(
    Promoter = if (!"Promoter" %in% names(.)) 0L else Promoter,
    Distal   = if (!"Distal"   %in% names(.)) 0L else Distal,
    a = Promoter,
    b = Distal,
    c = N_prom - a,
    d = N_dist - b,
    total_hits = a + b
  ) %>%
  filter(total_hits >= min_hits)

# -----------------------------
# 5) Fisher test per family + BH
# -----------------------------
res <- fam_counts %>%
  rowwise() %>%
  mutate(stats = list(fisher_one(a, b, c, d))) %>%
  unnest(stats) %>%
  ungroup() %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
    stars = sig_stars(q_value),
    log2_or = log2(odds_ratio)
  ) %>%
  arrange(q_value)

fwrite(res, out_table, sep = "\t")

# -----------------------------
# 6) Forest plot (top families)
# -----------------------------
plot_df <- res %>%
  slice_head(n = top_n) %>%
  mutate(family = factor(family, levels = rev(family)))

plot <- ggplot(plot_df, aes(x = odds_ratio, y = family)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.2) +
  geom_point(size = 2.4) +
  geom_text(aes(label = stars), hjust = -0.25, vjust = 0.3, size = 4) +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.25))) +
  labs(
    x = "Odds ratio (Promoter vs Distal)",
    y = NULL,
    title = "TF family promoter enrichment among motif-hit variants",
    subtitle = paste0("Fisher's exact test per family; universe = all tested variants; top ", top_n,
                      " families; stars show BH q-value")
  ) +
  theme_classic(base_size = 12)

ggsave(out_plot, p, width = 7.5, height = 5.5, useDingbats = FALSE)

message("Universe: Promoter N=", N_prom, " | Distal N=", N_dist)
message("Wrote: ", out_table)
message("Wrote: ", out_plot)


################################
## Quantifying accessibility in ChromBPNet peaks
## Author: Marliette Matos
## Date: 10/13/2025
################################
library(ggplot2)
library(readr)
library(dplyr)
library(cowplot)
library(GenomicRanges)
library(edgeR)
library(rlang)
library(cowplot)
library(tidyr)

source("~/cd4_qtl_paper_figures/utils/color_pallete_helper.R")

#read finemapping QTL results
summary_df<-read_delim("/gchm/cd4_qtl_paper_figures/figure_1/data/eqtl_caqtl_finemapping_coloc_all.tsv")
#read chromBPNet results
cBPNet_var<-read_delim("~/cd4_qtl_paper_figures/figure_4/data/motif_variant_overlap/cd4_top_cpbnet_variants_motif_overlap_moreTFS_boolean.tsv")
#read featurecounts for inout peaks for chrombpnet 
featcounts<-read_delim("~/cd4_chrombpnet/data/inputs/peaks/cd4_atac_8_samples_narrowPeaks.featureCounts.txt")
featcoords<-read_delim("~/cd4_chrombpnet/data/inputs/peaks/cd4_atac_8_samples_annotation.txt")

#make granges for chrombpnet peaks

featcoords_gr <- makeGRangesFromDataFrame(
  featcoords,
  seqnames.field = "Chr",
  start.field    = "Start",
  end.field      = "End",
  strand.field   = "Strand",
  keep.extra.columns = TRUE
)

featcoords_gr
### what fraction of chromBPNet variants are any type of QTLs?
##### since a single ChromBPNet variant can show up in multiple credible sets and multiple QTL types we are going global QTL appearance
#--- summary_df contains the all finemapped variants for both eQTLs and caQTLs
# ChromBPNet variants universe
chrombpnet_vars <- cBPNet_var$variant_id 
## ---- ca side: all variants that are in any caQTL credible set ----
ca_all <- summary_df %>%
  mutate(
    var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])",
                  "chr\\1_\\2_\\3_\\4",
                  variant_id.x)
  ) %>%
  semi_join(cBPNet_var, by = c("var_id" = "variant_id")) %>%
  pull(var_id) %>%
  unique()

## ---- eq side: all variants that are in any eQTL credible set ----
eq_all <- summary_df %>%
  mutate(
    var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])",
                  "chr\\1_\\2_\\3_\\4",
                  variant_id.y)
  ) %>%
  semi_join(cBPNet_var, by = c("var_id" = "variant_id")) %>%
  pull(var_id) %>%
  unique()

rm(summary_df)

qtl_any <- union(ca_all, eq_all)

# find 

variant_status <- tibble(
  variant_id = chrombpnet_vars
) %>%
  mutate(
    is_qtl = variant_id %in% qtl_any,
    category = ifelse(is_qtl, "QTL_shared", "ChromBPNet_specific")
  )

variant_status
fwrite(variant_status, "~/cd4_qtl_paper_figures/figure_4/data/chrombpnet_QTL_status.csv")

## make granges for variants
variant_status <- variant_status %>%
  tidyr::separate(
    col  = variant_id,
    into = c("Chr", "Pos", "Ref", "Alt"),
    sep  = "_",
    remove = FALSE
  ) %>%
  mutate(Pos = as.integer(Pos))

variant_status_gr <- makeGRangesFromDataFrame(
  variant_status,
  seqnames.field      = "Chr",
  start.field         = "Pos",
  end.field           = "Pos",
  keep.extra.columns  = TRUE
)

# Build a mapping GeneID -> coverage (1 row per peak id)
peak_cov_map <- featcounts %>%
  select(GeneID, atac_cov = merged_cd4_atac.bam) %>%
  distinct(GeneID, .keep_all = TRUE)

# Match peaks in GRanges to coverage map
idx <- match(mcols(featcoords_gr)$GeneID, peak_cov_map$GeneID)

# If you expect *every* peak to have coverage, enforce it; otherwise skip this stopifnot
# stopifnot(!any(is.na(idx)))

mcols(featcoords_gr)$atac_cov <- peak_cov_map$atac_cov[idx]

####----strip the suffix from the peak name (because they are narrow peaks)
strip_peak_suffix <- function(x) {
  sub("([0-9]+)[A-Za-z]$", "\\1", x)
}
# in featcoords_gr
#mcols(featcoords_gr)$peak_base <- strip_suffix(mcols(featcoords_gr)$GeneID)

#find overlaps between variants and cpbnet peaks
# Variant ↔ peak overlaps (SNPs in peaks)

# ----------------------------
#  Long table: variant × peak overlaps
# ----------------------------
hits <- findOverlaps(variant_status_gr, featcoords_gr, ignore.strand = TRUE)
q <- queryHits(hits)
s <- subjectHits(hits)

var_df <- as.data.frame(mcols(variant_status_gr))[q, , drop = FALSE] %>%
  mutate(
    var_chr = as.character(seqnames(variant_status_gr))[q],
    var_pos = start(variant_status_gr)[q]
  )

peak_df0 <- as.data.frame(mcols(featcoords_gr))[s, , drop = FALSE] %>%
  mutate(
    peak_chr   = as.character(seqnames(featcoords_gr))[s],
    peak_start = start(featcoords_gr)[s],
    peak_end   = end(featcoords_gr)[s],
    peak_width = width(featcoords_gr)[s]
  )

ov_tbl <- bind_cols(var_df, peak_df0) %>%
  mutate(
    var_offset0 = var_pos - peak_start,
    var_frac    = var_offset0 / peak_width
  )

stopifnot(nrow(ov_tbl) == length(hits))

# ----------------------------
# 2) Collapse: one row per variant × (collapsed) peak
# ----------------------------
ov_tbl2 <- ov_tbl %>%
  mutate(
    peak_base = strip_peak_suffix(GeneID)
  ) %>%
  group_by(variant_id, peak_base) %>%
  summarise(
    peak_chr   = dplyr::first(peak_chr),
    peak_start = dplyr::first(peak_start),
    peak_end   = dplyr::first(peak_end),
    peak_width = dplyr::first(peak_width),
    
    var_chr    = dplyr::first(var_chr),
    var_pos    = dplyr::first(var_pos),
    Ref        = dplyr::first(Ref),
    Alt        = dplyr::first(Alt),
    is_qtl     = dplyr::first(is_qtl),
    category   = dplyr::first(category),
    
    atac_cov   = dplyr::first(atac_cov),
    
    n_rows_collapsed = dplyr::n(),
    GeneID_collapsed = paste(unique(GeneID), collapse = ";"),
    .groups = "drop"
  )

# ----------------------------
# 3) Peak-level table (each peak counted once)
#    Rule: peaks with BOTH categories are counted as QTL_shared ("QTL_shared_any")
# ----------------------------
ov_tbl2 <- ov_tbl2 %>%
  mutate(
    category = trimws(category),
    category = recode(category,
                      ChromBPNet_specific = "ChromBPNet_spec",
                      ChromBPNet_spec     = "ChromBPNet_spec",
                      QTL_shared          = "QTL_shared"
    )
  )

peak_group <- ov_tbl2 %>%
  distinct(peak_base, category) %>%
  group_by(peak_base) %>%
  summarise(
    has_qtl  = any(category == "QTL_shared"),
    has_spec = any(category == "ChromBPNet_spec"),
    peak_group = case_when(
      has_qtl              ~ "QTL_shared_any",        # includes “Both”
      !has_qtl & has_spec  ~ "ChromBPNet_spec_only",
      TRUE                 ~ "None"                   # catch-all for unexpected categories
    ),
    .groups = "drop"
  )

table(peak_group$peak_group, useNA = "ifany")

peak_cov <- ov_tbl2 %>%
  distinct(peak_base, atac_cov, peak_width) %>%
  mutate(
    atac_per_bp       = atac_cov / peak_width,
    atac_per_bp_log10 = log10(atac_per_bp + 1)
  )

peak_df <- peak_cov %>%
  left_join(peak_group, by = "peak_base") %>%
  filter(!is.na(peak_group)) %>%
  mutate(peak_group = factor(peak_group, levels = c("ChromBPNet_spec_only", "QTL_shared_any")))

# ----------------------------
# 4) Test + plot
# ----------------------------
tt <- t.test(atac_per_bp_log10 ~ peak_group, data = peak_df)
p_lab <- paste0("p = ", signif(tt$p.value, 3))

y_max <- max(peak_df$atac_per_bp_log10, na.rm = TRUE)

pt <- ggplot(peak_df, aes(x = peak_group, y = atac_per_bp_log10, fill = peak_group)) +
  geom_violin() +
  annotate("text", x = 1.5, y = y_max * 1.05, label = p_lab, size = 6) +
  labs(x = NULL, y = "log10(mean ATAC counts per bp)") +
  theme_bw() +
  scale_fill_manual(values = variant_cmap3, drop = FALSE) +
  theme(legend.position = "none")

pt




ggsave("~/cd4_qtl_paper_figures/figure_4/plotting/general_figures/cbpnet_qtl_variant_coverage.pdf", pt, width =4, height = 5)


##########################
##########################
#########################
#Comparing other metrics
# ----------------------------
# 0) Read ChromBPNet metrics
# ----------------------------
cBPNet_var <- read_delim(
  "~/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_scores/average_cd4_tcells_AJ_common_variants.mean.variant_scores_IPS05.tsv"
)

# ----------------------------
# 1) Build a VARIANT-level table (no peak duplication)
#    Use ONE row per variant_id, with category
# ----------------------------
variant_df <- ov_tbl2 %>%
  distinct(variant_id, category)

# If you prefer to count "Both" peaks as QTL_shared at the peak-level,
# that DOES NOT change variant category. So keep variant category as-is.

# Standardize category labels (important!)
variant_df <- variant_df %>%
  mutate(
    category = trimws(category),
    category = recode(category,
                      ChromBPNet_specific = "ChromBPNet_spec",
                      ChromBPNet_spec     = "ChromBPNet_spec",
                      QTL_shared          = "QTL_shared"
    ),
    category = factor(category, levels = c("ChromBPNet_spec", "QTL_shared"))
  )

# Join ChromBPNet scores (now 1 row per variant)
variant_metrics <- variant_df %>%
  left_join(cBPNet_var, by = "variant_id")

# ----------------------------
# 2) Helpers
# ----------------------------
p_to_asterisk <- function(p) {
  if (is.na(p)) return("n/a")
  if (p < 1e-4) return("****")
  if (p < 1e-3) return("***")
  if (p < 1e-2) return("**")
  if (p < 5e-2) return("*")
  return("n.s.")
}

metrics <- c(
  "abs_logfc.mean",
  "jsd.mean",
  "abs_logfc_x_jsd.mean",
  "active_allele_quantile.mean",
  "abs_logfc_x_active_allele_quantile.mean",
  "jsd_x_active_allele_quantile.mean",
  "abs_logfc_x_jsd_x_active_allele_quantile.mean",
  "abs_quantile_change.mean"
)

pseudo <- 1e-4
pt.list <- list()

# ----------------------------
# 3) Plot each metric
# ----------------------------
for (metric in metrics) {
  
  if (!metric %in% names(variant_metrics)) {
    message("Skipping missing metric: ", metric)
    next
  }
  
  df_metric <- variant_metrics %>%
    transmute(
      category,
      value = .data[[metric]],
      value_log = log10(value + pseudo)
    ) %>%
    filter(!is.na(value_log), !is.na(category))
  
  # t-test on transformed values (2 groups)
  tt <- t.test(value_log ~ category, data = df_metric)
  
  asterisk <- p_to_asterisk(tt$p.value)
  p_lab    <- paste0("p = ", signif(tt$p.value, 3))
  
  y_max <- max(df_metric$value_log, na.rm = TRUE)
  
  pt <- ggplot(df_metric, aes(x = category, y = value_log, fill = category)) +
    geom_violin() +
    scale_fill_manual(values = variant_cmap4, drop = FALSE) +
    annotate("text", x = 1.5, y = y_max * 1.10, label = asterisk,
             size = 3.5, fontface = "bold") +
    annotate("text", x = 1.5, y = y_max * 1.02, label = p_lab,
             size = 3.3) +
    labs(
      title = metric,
      x = NULL,
      y = paste0("log10(value + ", pseudo, ")")
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title      = element_text(size = 20),
      axis.text.y     = element_text(size = 15),
      axis.text.x     = element_text(size = 15),
      plot.title      = element_text(size = 12)
    )
  
  pt.list[[metric]] <- pt
}

pt.metrics <- cowplot::plot_grid(
  plotlist = pt.list,
  ncol = 4,
  align = "hv"
)

pt.metrics

ggsave("~/cd4_qtl_paper_figures/figure_4/plotting/general_figures/cbpnet_qtl_variant_bpnet_metric.pdf", pt.metrics, width =15, height = 10)



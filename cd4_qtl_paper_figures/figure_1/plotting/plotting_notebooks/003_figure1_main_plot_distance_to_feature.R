#!/usr/bin/env Rscript
### Continuous distance analysis for eQTL vs caQTL + caQTL categories
### Author: Marliette Matos
### Date: 01/17/2025 

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(GenomicRanges)
  library(IRanges)
  library(readr)
  library(ggpubr)
  library(rstatix)
  library(patchwork)
})

# ------------------------------------------------------------
# 0) Inputs
# ------------------------------------------------------------
source("cd4_qtl_paper_figures/utils/color_pallete_helper.R")

summary_df   <- "/gchm/cd4_qtl_paper_figures/figure_1/data/eqtl_caqtl_finemapping_coloc_all.tsv"
gencode_bed  <- "~/resources/genome/hg38_gencode_raw/gencode.v44.primary_assembly.annotation.genes.bed"
peak_ranges  <- "/gchm/ATAC-seq_analysis/diff_accesibility_ana/results/peak_counts/RAW_cd4_atac_peakcounts_ranges_scrna_union.csv"

out_pdf <- "~/cd4_qtl_paper_figures/figure_1/plotting/plots_oct2025/mean_distance_to_egene_capeak_stats.pdf"

# groups used earlier
eqtls  <- c("coloc", "eQTL_only")
caqtls <- c("coloc", "caQTL_only")

# ------------------------------------------------------------
# 1) Load + split summary table
# ------------------------------------------------------------
df0 <- fread(summary_df) %>% as_tibble()

eqtl_df0  <- df0 %>% filter(coloc_status %in% eqtls)
caqtl_df0 <- df0 %>% filter(coloc_status %in% caqtls)
rm(df0)

# ------------------------------------------------------------
# 2) Build GENCODE TSS table (per gene)
# ------------------------------------------------------------
bed <- fread(gencode_bed, sep = "\t", header = FALSE, showProgress = FALSE)

setnames(
  bed,
  old = names(bed)[1:min(6, ncol(bed))],
  new = c("chr", "start0", "end", "name", "score", "strand")[1:min(6, ncol(bed))]
)
if (!"strand" %in% names(bed)) bed$strand <- "*"
bed[, gene_id := if (ncol(bed) >= 7) as.character(bed[[7]]) else NA_character_]

gencode_df <- bed %>%
  transmute(
    chr,
    start_1b   = as.integer(start0) + 1L,
    end_1b     = as.integer(end),
    gene       = as.character(name),
    strand     = as.character(strand),
    gene_id    = gene_id,
    tss_1based = if_else(strand == "+", start_1b, end_1b)
  ) %>%
  filter(str_detect(chr, "^chr(\\d+|X|Y)$")) %>%
  group_by(gene, chr) %>%
  slice_max(end_1b - start_1b, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  distinct(gene, .keep_all = TRUE)

# ------------------------------------------------------------
# 3) eQTL: per-CS mean distance to eGene TSS
# ------------------------------------------------------------
tss_dist_df <- eqtl_df0 %>%
  select(gene, finemapped_cs_eqtl, chr.y, variant_pos.y) %>%
  distinct() %>%
  left_join(gencode_df %>% select(gene, chr_tss = chr, strand, tss_1based), by = "gene") %>%
  filter(!is.na(tss_1based), chr.y == chr_tss) %>%
  mutate(
    signed_tss_dist = if_else(strand == "+", variant_pos.y - tss_1based, tss_1based - variant_pos.y),
    abs_tss_dist    = abs(signed_tss_dist)
  ) %>%
  group_by(finemapped_cs_eqtl) %>%
  summarise(mean_tss_dist = mean(abs_tss_dist), .groups = "drop")

# ------------------------------------------------------------
# 4) caQTL: per-CS mean distance to caPeak center (by category)
# ------------------------------------------------------------

# peak ranges table -> GRanges
ranges.table <- read.csv2(peak_ranges, sep = ",") %>%
  filter(Chr %in% paste0("chr", 1:22))
rownames(ranges.table) <- ranges.table$X

peak_gr <- GRanges(
  seqnames = ranges.table$Chr,
  ranges   = IRanges(ranges.table$Start, ranges.table$End),
  mcols    = data.frame(peakID = rownames(ranges.table))
)
names(peak_gr) <- peak_gr$mcols.peakID
peak_gr$peakID <- names(peak_gr)

# normalize caQTL category labels
caqtl_df0 <- caqtl_df0 %>%
  mutate(
    caqtl_category = ifelse(
      cs_type_any_var %in% c("in_corr_Peak", "in_other_Peak"),
      "in_other_Peak",
      cs_type_any_var
    ),
    caqtl_category = factor(caqtl_category,
                            levels = c("no_Peak_overlap", "in_other_Peak", "in_caPeak"))
  )

caqtl_dist_df <- caqtl_df0 %>%
  filter(!is.na(finemapped_cs_caqtl), !is.na(variant_id.x),
         !is.na(peak), !is.na(chr.x), !is.na(variant_pos.x)) %>%
  distinct(finemapped_cs_caqtl, peak, variant_id.x, chr.x, variant_pos.x, caqtl_category) %>%
  left_join(as.data.frame(peak_gr), by = c("peak" = "peakID")) %>%
  filter(!is.na(start), !is.na(end), as.character(seqnames) == chr.x) %>%
  mutate(
    peak_center  = start + (end - start) / 2,
    abs_distance = abs(variant_pos.x - peak_center)
  ) %>%
  group_by(finemapped_cs_caqtl, caqtl_category) %>%
  summarise(mean_peak_dist = mean(abs_distance), .groups = "drop")

# ------------------------------------------------------------
# 5) Continuous combined df + log2 transform
# ------------------------------------------------------------
df_eqtl <- tss_dist_df %>%
  transmute(
    type = "eQTL",
    caqtl_category = NA_character_,
    dist = mean_tss_dist
  )

df_caqtl <- caqtl_dist_df %>%
  transmute(
    type = "caQTL",
    caqtl_category = as.character(caqtl_category),
    dist = mean_peak_dist
  )

plot_df <- bind_rows(df_eqtl, df_caqtl) %>%
  mutate(
    type = factor(type, levels = c("eQTL", "caQTL")),
    caqtl_category = factor(caqtl_category,
                            levels = c("no_Peak_overlap", "in_other_Peak", "in_caPeak")),
    log2dist = log2(dist + 1)
  )

# ------------------------------------------------------------
# 6) Stats + plots
# ------------------------------------------------------------

# ---- (A) eQTL vs caQTL (all pooled) ----
test1 <- wilcox.test(log2dist ~ type, data = plot_df, exact = FALSE)

p_lab <- if (test1$p.value == 0) {
  "Wilcoxon: p < 1e-300"
} else {
  paste0("Wilcoxon: p = ", format.pval(test1$p.value, digits = 3, eps = 1e-300))
}

plot_df2 <- plot_df %>%
  mutate(
    type_lab = recode(as.character(type),
                      "eQTL" = "eQTLs",
                      "caQTL" = "caQTLs"),
    type_lab = factor(type_lab, levels = c("eQTLs", "caQTLs"))
  )

ymax1 <- max(plot_df2$log2dist, na.rm = TRUE)

p1 <- ggplot(plot_df2, aes(x = type_lab, y = log2dist, fill = type_lab)) +
  geom_violin(trim = TRUE, linewidth = 0.35) +
  geom_boxplot(width = 0.16, outlier.size = 0.3, linewidth = 0.35, color = "black") +
  scale_fill_manual(values = qtl_cmap, guide = "none") +
  annotate("text", x = 1.5, y = ymax1 + 0.35, label = p_lab, size = 5) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = NULL, y = "log2(mean distance to feature + 1)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(size = 16, color = "black"),
    axis.text.y  = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.line    = element_line(linewidth = 0.4),
    axis.ticks   = element_line(linewidth = 0.4)
  )

# ---- (B) within caQTLs by category ----
caqtl_only <- plot_df %>%
  filter(type == "caQTL") %>%
  mutate(
    caqtl_category = factor(as.character(caqtl_category),
                            levels = rev(c("no_Peak_overlap", "in_other_Peak", "in_caPeak")))
  )

# pairwise Wilcoxon + BH (for brackets)
pw_tbl <- caqtl_only %>%
  pairwise_wilcox_test(log2dist ~ caqtl_category, p.adjust.method = "BH") %>%
  mutate(
    p_label = ifelse(p.adj == 0, "BH p < 1e-300",
                     paste0("BH p = ", format.pval(p.adj, digits = 3, eps = 1e-300)))
  )

# y positions for brackets
ymax2 <- max(caqtl_only$log2dist, na.rm = TRUE)
step  <- 0.35

pw_tbl <- pw_tbl %>%
  mutate(
    group1 = as.character(group1),
    group2 = as.character(group2),
    y.position = ymax2 + step * row_number()
  )

kw <- kruskal.test(log2dist ~ caqtl_category, data = caqtl_only)
kw_lab <- if (kw$p.value == 0) {
  "Kruskal–Wallis: p < 1e-300"
} else {
  paste0("Kruskal–Wallis: p = ", format.pval(kw$p.value, digits = 3, eps = 1e-300))
}

p2 <- ggplot(caqtl_only, aes(x = caqtl_category, y = log2dist, fill = caqtl_category)) +
  geom_violin(trim = TRUE, linewidth = 0.35) +
  geom_boxplot(width = 0.16, outlier.size = 0.3, linewidth = 0.35, color = "black") +
  scale_fill_manual(values = caQTL_category_cmap2, guide = "none") +
  stat_pvalue_manual(
    pw_tbl,
    label = "p_label",
    y.position = "y.position",
    xmin = "group1",
    xmax = "group2",
    tip.length = 0.01,
    size = 5,
    inherit.aes = FALSE
  ) +
  annotate(
    "text",
    x = 2,
    y = ymax2 + step * (nrow(pw_tbl) + 1) + 0.05,
    label = kw_lab,
    size = 5
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  labs(x = NULL, y = "log2(mean distance + 1)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(size = 16, color = "black"),
    axis.text.y  = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    axis.line    = element_line(linewidth = 0.4),
    axis.ticks   = element_line(linewidth = 0.4)
  )

# ------------------------------------------------------------
# 7) Save combined figure
# ------------------------------------------------------------
final_plot <- p1 + p2 + plot_layout(heights = c(0.8, 1))
ggsave(out_pdf, final_plot, width = 5, height = 10)

# optional: print objects if running interactively
print(final_plot)

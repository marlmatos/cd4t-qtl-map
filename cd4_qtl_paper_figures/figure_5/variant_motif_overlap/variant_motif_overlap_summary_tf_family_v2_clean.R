## =========================
## Script summarizing variant-motif overlap
## =========================
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)
library(GenomicRanges)
library(scales)
library(data.table)
library(purrr)
library(stringr)
library(forcats)
library(RColorBrewer)
library(ggseqlogo)
options(bitmapType = "cairo")
source("~/cd4_qtl_paper_figures/utils/color_pallete_helper.R")

## =========================
## 1) Load motif–variant hits
## =========================

unique_hits <- read_tsv(
  "~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap_hit_caller/filtered_variant_hits.tsv"
)

## Order TF families by # of variants with motif hits (ChromBPNet universe)
tf_counts <- unique_hits %>%
  group_by(family) %>%
  summarise(
    total_hits = n_distinct(variant_loc),   # unique variants per family
    .groups = "drop"
  ) %>%
  arrange(desc(total_hits))

tf_order  <- tf_counts$family
n_colors  <- length(tf_order)
my_colors <- colorRampPalette(brewer.pal(8, "Blues"))(n_colors)

## =========================
## 2) Load ChromBPNet predictions
## =========================

var_prediction <- read_tsv(
  "~/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_scores/average_cd4_tcells_AJ_common_variants.mean.variant_scores_IPS05.tsv"
) %>%
  mutate(var_chr_pos = paste(chr, pos, sep = "_")) %>%
  # keep only variants present in unique_hits
  semi_join(
    unique_hits %>%
      transmute(var_chr_pos = paste(seqnames, variant_loc, sep = "_")),
    by = "var_chr_pos"
  ) %>%
  dplyr::select(variant_id, chr, pos, allele1, allele2, logfc.mean, jsd.mean, var_chr_pos) %>%
  group_by(var_chr_pos) %>%
  filter(n() == 1) %>%        # drop positions with multiple prediction rows
  ungroup()

## =========================
## 3) QTL overlap annotation
## =========================

summary_df <- read_delim(
  "/gpfs/commons/home/mmatos/cd4_qtl_paper_figures/figure_1/data/eqtl_caqtl_finemapping_coloc_all.tsv"
)

chrombpnet_vars <- var_prediction$variant_id

## caQTL side
ca_all <- summary_df %>%
  mutate(
    var_id = gsub(
      "([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])",
      "chr\\1_\\2_\\3_\\4",
      variant_id.x
    )
  ) %>%
  filter(var_id %in% chrombpnet_vars) %>%
  pull(var_id) %>%
  unique()

## eQTL side
eq_all <- summary_df %>%
  mutate(
    var_id = gsub(
      "([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])",
      "chr\\1_\\2_\\3_\\4",
      variant_id.y
    )
  ) %>%
  filter(var_id %in% chrombpnet_vars) %>%
  pull(var_id) %>%
  unique()

rm(summary_df)

qtl_any <- union(ca_all, eq_all)

variant_status <- tibble(
  variant_id = chrombpnet_vars
) %>%
  mutate(
    is_qtl  = variant_id %in% qtl_any,
    category = if_else(is_qtl, "QTL_shared", "ChromBPNet_specific")
  )

var_prediction <- var_prediction %>%
  left_join(variant_status, by = "variant_id")

## =========================
## 4) Peak annotation (promoter / distal)
## =========================

## ATAC peaks used by ChromBPNet
bp.peaks <- fread(
  "~/cd4_chrombpnet/data/inputs/peaks/merged_cd4_samples_peaks_no_blacklist.sorted.bed"
)
bp.peaksGR <- GRanges(
  seqnames = bp.peaks$V1,
  ranges   = IRanges(bp.peaks$V2, bp.peaks$V3),
  peak_name = bp.peaks$V4
)

## Gencode gene TSSs
gencode <- fread(
  "~/resources/genome/hg38_gencode_raw/gencode.v44.primary_assembly.annotation.genes.bed",
  sep = "\t", data.table = FALSE
) %>%
  transmute(
    chr        = V1,
    strand     = V6,
    gene_id    = V7,
    gene       = V4,
    tss_1based = if_else(V6 == "+", V2 + 1L, V3)
  ) %>%
  distinct(chr, gene, .keep_all = TRUE)

tss_gr <- with(
  gencode,
  GRanges(
    seqnames = chr,
    ranges   = IRanges(start = tss_1based, end = tss_1based),
    strand   = strand
  )
)

ov_tss <- countOverlaps(bp.peaksGR, tss_gr, ignore.strand = TRUE)
mcols(bp.peaksGR)$promoter   <- ov_tss > 0
mcols(bp.peaksGR)$annotation <- if_else(mcols(bp.peaksGR)$promoter, "promoter", "distal")

## Map variant_id → peak_annotation
vars_GR <- variant_status %>%
  separate(
    col   = variant_id,
    into  = c("chr", "pos", "ref", "alt"),
    sep   = "_",
    remove = FALSE
  ) %>%
  mutate(pos = as.integer(pos)) %>%
  {
    GRanges(
      seqnames   = .$chr,
      ranges     = IRanges(start = .$pos, end = .$pos),
      strand     = "*",
      is_qtl     = .$is_qtl,
      category   = .$category,
      ref        = .$ref,
      alt        = .$alt,
      variant_id = .$variant_id
    )
  }

hits <- findOverlaps(vars_GR, bp.peaksGR, ignore.strand = TRUE)

ann_vec  <- rep(NA_character_, length(vars_GR))
name_vec <- rep(NA_character_, length(vars_GR))

ann_vec[queryHits(hits)]  <- mcols(bp.peaksGR)$annotation[subjectHits(hits)]
name_vec[queryHits(hits)] <- mcols(bp.peaksGR)$peak_name[subjectHits(hits)]

mcols(vars_GR)$peak_annotation <- ann_vec
mcols(vars_GR)$peak_name       <- name_vec

var_tss <- tibble(
  variant_id     = vars_GR$variant_id,
  peak_annotation = vars_GR$peak_annotation
)

var_prediction <- var_prediction %>%
  left_join(var_tss, by = "variant_id")

## =========================
## 5) Base per-variant–motif table (for plots)
## =========================

base_df <- unique_hits %>%
  mutate(var_chr_pos = paste(seqnames, variant_loc, sep = "_")) %>%
  left_join(var_prediction, by = "var_chr_pos") %>%
  filter(!if_any(c(logfc.mean, jsd.mean), is.na)) %>%
  distinct(variant_id, motif_name, .keep_all = TRUE) %>%    # 1 row per variant–motif
  mutate(
    TF_family = factor(family, levels = rev(tf_order))
  )

## =========================
## 6) Allelic importance Δ (allele2 - allele1)
## =========================

allelic_ratio_df <- unique_hits %>%
  mutate(var_chr_pos = paste(seqnames, variant_loc, sep = "_")) %>%
  left_join(var_prediction, by = "var_chr_pos") %>%
  filter(!if_any(c(logfc.mean, jsd.mean), is.na)) %>%
  dplyr::select(variant_id, motif_name, family, allele, hit_importance) %>%
  group_by(variant_id, motif_name, family, allele) %>%
  summarise(
    hit_importance = mean(hit_importance),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = allele,          # assumes columns "allele1", "allele2"
    values_from = hit_importance
  ) %>%
  mutate(
    allele1    = replace_na(allele1, 0),
    allele2    = replace_na(allele2, 0),
    allele_diff = if_else(
      allele1 == 0 & allele2 == 0,
      NA_real_,
      allele2 - allele1
    )
  ) %>%
  filter(!is.na(allele_diff)) %>%
  mutate(
    TF_family = factor(family, levels = rev(tf_order))
  )

## =========================
## 7) Per-family master summary table
## =========================
write_csv(base_df, "~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap_hit_caller/motif_hits_summary_annotations.csv")

## Core family-level stats from base_df
family_core <- base_df %>%
  group_by(family) %>%
  summarise(
    n_var_motif   = n(),                               # variant–motif instances
    n_variants    = n_distinct(variant_id),            # distinct variants with motifs
    median_logfc  = median(logfc.mean, na.rm = TRUE),
    n_qtl         = sum(is_qtl, na.rm = TRUE),
    frac_qtl      = n_qtl / n_variants,
    n_promoter    = sum(peak_annotation == "promoter", na.rm = TRUE),
    n_distal      = sum(peak_annotation == "distal",   na.rm = TRUE),
    frac_promoter = n_promoter / n_variants,
    frac_distal   = n_distal   / n_variants,
    .groups = "drop"
  )

## Allelic importance summary
allelic_summary <- allelic_ratio_df %>%
  group_by(family) %>%
  summarise(
    median_allele_diff = median(allele_diff, na.rm = TRUE),
    .groups = "drop"
  )

## Master summary table
family_summary <- family_core %>%
  left_join(allelic_summary,   by = "family") %>%
  left_join(tf_counts,   by = "family") %>%
  mutate(
    TF_family = factor(family, levels = rev(tf_order))
  ) %>%
  arrange(TF_family)

write_tsv(
  family_summary,
  "~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap_hit_caller/motif_family_summary_table.tsv"
)


write_tsv(
  var_prediction,
  "~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap_hit_caller/unique_variants_overlapping_motifs.tsv"
)

family_summary<-read_tsv("~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap_hit_caller/motif_family_summary_table.tsv")

mean(family_summary$frac_qtl)
## =========================
## 8) Example: plots re-using family_summary / base_df
## =========================

label_plot <- family_summary %>%
  ggplot(aes(y = TF_family, x = 2, label = TF_family)) +
  geom_text(hjust = 1, size = 4) +
  theme_void() +
  xlim(c(0, 2)) +
  coord_cartesian(clip = "off") +  # prevent cropping of text
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )



#just plot the weight matrix for a representiave motif pattern for each family
# 1) representative motif -> family map (fixed a missing comma after TYY1)
rep_map <- c(
  "ATF1_CREB1"              = "AP-1/ATF/CREB (bZIP)",
  "CTCF"                    = "CTCF",
  "PAX5"                    = "PAX",
  "KLF12"                   = "KLF (C2H2 ZF)",
  "RUNX1_BCL11B"            = "Mixed: RUNX/BTB-ZF",
  "NRF1"                    = "NRF (bZIP)",
  "NFYA_B_C"                = "NF-Y (CCAAT)",
  "ETS1"                    = "ETS",
  "TYY1"                    = "YY (C2H2 ZF)",
  "IRF1_IRF2"               = "IRF",
  "ZBTB33_KAISO"            = "BTB/POZ ZF",
  "PUTATIVE_ETS1_IKZF3_ERG" = "Mixed: ETS/IKAROS",
  "ZNF76_ZN143_1"           = "ZNF (C2H2 ZF)",
  "TFE3"                    = "MiT/TFE (bHLH-ZIP)",
  "RFX3_RFX5_RFX1"          = "RFX",
  "PUTATIVE_RORG"           = "Nuclear receptor (ROR)",
  "TCF7_LEF1"               = "TCF/LEF (HMG)",
  "NFKB1"                   = "NF-κB (Rel)",
  "SP2_SP3"                 = "Sp/KLF (C2H2 ZF)",
  "PUTATIVE_ETS_IRF"        = "Mixed: ETS/IRF"
)
motif_list<-readRDS("~/cd4_chrombpnet/plotting/motifs_denovo_fwd_wcm.rds")

#motif_list <- motif_list[tf_order[tf_order %in% names(motif_list)]]

# --- 2) build family-ordered list by TF names, then rename to families ---
# invert mapping: family -> representative motif name
family_to_motif <- setNames(names(rep_map), unname(rep_map))
family_to_motif <- family_to_motif[tf_order[tf_order %in% names(family_to_motif)]]

# desired family order: keep your tf_order, then any extra families from rep_map
families_in_map <- names(family_to_motif)
family_order <- c(intersect(tf_order, families_in_map),
                  setdiff(families_in_map, tf_order))

motifs_needed <- unname(family_to_motif[family_order])
present <- motifs_needed %in% names(motif_list)

# --- optional: fallback if representative motif is missing but you have tf_map ---
logo_strength <- function(mat) sum(abs(mat))  # good for CWMs
if (any(!present) && exists("tf_map")) {
  missing_fams <- family_order[!present]
  # candidates among *available* motifs that belong to the same family
  # tf_map is motif -> family
  available <- intersect(names(tf_map), names(motif_list))
  cand_tbl <- tibble(
    motif  = available,
    family = unname(tf_map[available]),
    strength = map_dbl(available, ~ logo_strength(motif_list[[.x]]))
  )
  for (fam in missing_fams) {
    cand <- cand_tbl %>% filter(family == fam)
    if (nrow(cand) > 0) {
      best <- cand %>% slice_max(strength, n = 1, with_ties = FALSE) %>% pull(motif)
      motifs_needed[family_order == fam] <- best
      present[family_order == fam] <- TRUE
      message("Fallback for '", fam, "': using motif '", best, "'.")
    }
  }
}

# drop families still missing
if (any(!present)) {
  warning("Dropped families without available motif matrices: ",
          paste(family_order[!present], collapse = ", "))
}
family_order_kept <- family_order[present]
motifs_kept <- motifs_needed[present]

rep_motif_list <- motif_list[motifs_kept]
#names(rep_motif_list) <- family_order_kept  # rename entries to FAMILY labels

# --- 3) plot ---
# If matrices are CWMs -> method = "custom"; if PWMs -> method = "prob"
p <- ggseqlogo(
  rep_motif_list,
  method = "custom",
  ncol   = 1
) +
  theme(
    axis.text.x   = element_blank(),
    axis.text.y   = element_blank(),
    strip.text    = element_blank(),
    panel.spacing = unit(0.001, "lines")
  )

p
## (a) # variants per family
p_variants <- family_summary %>%
  ggplot(aes(y = TF_family, x = n_variants)) +
  geom_col(color = "black", fill = "#0071ba",linewidth = 0.2, width = 0.75) +
  geom_text(aes(label = n_variants), hjust = -0.3, size = 4) +
  labs(x = "# Variants") +
  theme_classic() +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = 9),
    axis.title   = element_text(size = 10),
    plot.title   = element_text(size = 12),
    legend.position = "none"
  ) +
  scale_x_continuous(expand = c(0, 0))

## (b) Allelic logFC boxplot (per variant–motif)
fc_boxplot <- ggplot(base_df, aes(x = logfc.mean, y = TF_family)) +
  geom_violin( fill = "#0071ba", color = "gray30", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "#c72e29", linewidth = 0.6) +
  labs(x = "Allelic logFC (Alt/Ref)") +
  theme_classic() +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = 9),
    axis.title   = element_text(size = 10),
    plot.title   = element_text(size = 12),
    legend.position = "none"
  )

## (c) Hit importance Δ boxplot
importance_boxplot <- ggplot(allelic_ratio_df,
                             aes(x = allele_diff, y = TF_family)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "gray30") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", linewidth = 0.4) +
  labs(x = "Global Hit Importance (allele2 - allele1)") +
  theme_classic() +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = 9),
    axis.title   = element_text(size = 10),
    plot.title   = element_text(size = 12),
    legend.position = "none"
  )

## (d) QTL overlap stacked (QTL vs ChromBPNet-specific)
sum_df_long <- family_summary %>%
  transmute(
    TF_family,
    n_qtl,
    n_nonqtl = n_variants - n_qtl
  ) %>%
  pivot_longer(
    cols      = c(n_qtl, n_nonqtl),
    names_to  = "var_type",
    values_to = "n"
  ) %>%
  mutate(
    bucket = recode(var_type,
                    n_qtl    = "QTL_shared",
                    n_nonqtl = "ChromBPNet_spec")
  )

motif_qtl.plt <- ggplot(sum_df_long, aes(x = TF_family, y = n, fill = bucket)) +
  geom_col(position = "fill", width = 0.75, color = "black", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = variant_cmap4, drop = FALSE) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = 9),
    axis.title   = element_text(size = 10),
    plot.title   = element_text(size = 12),
    legend.position = "none") +
  labs(x = NULL, y = "QTL overlap", fill = "Variant set")

## (e) Promoter vs distal stacked
stack_df <- family_summary %>%
  transmute(
    TF_family,
    n_promoter,
    n_distal
  ) %>%
  pivot_longer(
    cols      = c(n_promoter, n_distal),
    names_to  = "annotation",
    values_to = "n"
  ) %>%
  mutate(
    annotation = recode(annotation,
                        n_promoter = "promoter",
                        n_distal   = "distal")
  )

motif_tss.plt <- ggplot(stack_df,
                        aes(x = TF_family, y = n, fill = annotation)) +
  geom_col(position = "fill", width = 0.75, color = "black", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = c("promoter" = "#c72e29",
                               "distal"   = "#eaadad")) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = 9),
    axis.title   = element_text(size = 10),
    plot.title   = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(
    x    = NULL,
    y    = "Proportion of motif-overlapping variants",
    fill = "Genomic annotation"
  )

## replace motif_tss.plt with enrichment results

TF_promoter_enrich<-fread("~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/TF_motif_family_promoter_enrich_fisher_results.tsv") %>% select(c("family","odds_ratio", "log2_or", "conf_low", "conf_high", "q_value", "stars"))

colnames(TF_promoter_enrich)<-c("family", "fisher_odds_ratio" , "fisher_log2_or" ,  "fisher_conf_low" , "fisher_conf_high" ,     "fisher_q_value" , "fisher_stars")

enrich_df <-family_summary %>% left_join(TF_promoter_enrich, by="family") %>% select(c("TF_family", "fisher_odds_ratio" , "fisher_conf_low", "fisher_conf_high", "fisher_q_value", "fisher_stars"))

p.enrich <- ggplot(enrich_df, aes(x = fisher_odds_ratio, y = TF_family)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = fisher_conf_low, xmax = fisher_conf_high), height = 0.2) +
  geom_point(size = 1.4) +
  geom_text(aes(label = fisher_stars), hjust = -0.25, vjust = 0.3, size = 4) +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.25))) +
  labs(
    x = "Odds ratio (Promoter vs Distal)",
    y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = 9),
    axis.title   = element_text(size = 10),
    plot.title   = element_text(size = 12),
    legend.position = "none"
  )


grid<-plot_grid(label_plot, p, p_variants, fc_boxplot, p.enrich, motif_qtl.plt, rel_widths = c(.8, 1,0.8, 0.8, 1,0.8), nrow = 1,
                align = "h",
                axis = "tb",
                label_size = 10,
                labels =c("TF Family", "De Novo Pattern", "Ovelapping Variants", "CA Allelic Effect", "CRE Overlap", "QTL overlap"))
grid


ggsave("~/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/variant_motif_overlap/plots/variant_motif_overlap_summary_hitcaller_v2.pdf", grid, width =8, height = 5)







# ============================================================
# Annotate peaks + genes, derive relative position features,
# and run enrichment tests (rel_pos, caqtl_category, ChromBPNet)
# Minimal changes; mostly re-ordered + de-duplicated libraries.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(data.table)
  library(purrr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(sandwich)  # vcovCL
  library(lmtest)    # coeftest
  library(rtracklayer)
  library(GenomicFeatures)
  library(txdbmaker)
})

options(bitmapType = "cairo")

# ------------------------------------------------------------
# 0) Inputs assumed to exist
# ------------------------------------------------------------
# - findr_res_sig : Findr results (must contain: type, phenotype_A, phenotype_B, locus_id, caqtl_category, ...)
setwd("~/cd4_qtl_paper_figures/figure_3")

findr_res<-fread("~/cd4_qtl_paper_figures/figure_3/data/findr_MEDIATION_All_results_collapsed_not_filtered.tsv")
findr_res_sig<-findr_res %>% filter(fdr_sig_5==TRUE)

## add information about the type of caQTL category
findr_res_sig <- findr_res_sig %>%
  mutate(
    caqtl_finemapped_cs = vapply(
      strsplit(finemapped_cs_coloc, "_"),
      function(x) paste(x[1:7], collapse = "_"),
      character(1)
    )
  )

caqttls<-read_csv("~/cd4_qtl_paper_figures/figure_1/data/CD4T_combined_finemapping_annotated.csv") %>% dplyr::select(finemapped_cs, cs_type_any_var)
caqttls<-caqttls[caqttls$finemapped_cs %in% findr_res_sig$caqtl_finemapped_cs, ] %>% distinct() %>%
  mutate(caqtl_category=ifelse(cs_type_any_var=="in_corr_Peak", "in_other_Peak", cs_type_any_var))

findr_res_sig<-findr_res_sig %>% left_join(caqttls, join_by("caqtl_finemapped_cs"=="finemapped_cs"))


# ============================================================
# 1) Add peak coordinates to findr_res_sig
# ============================================================

peak_coords <- read_tsv(
  "/gcgl/mmatos/cd4_aging_project/data/ATAC-seq/atac_preprocessing/merged_library/peak_calling/MACS3/BAMPE/peaks_102024/cd4_atac_padded_summits_peaks.bed",
  col_names = NULL,
  show_col_types = FALSE
)

peaks_names <- read.delim2(
  "/gchm/cd4_caQTL_analysis/variant_to_peak_QTL/run_012625_qc_aware_qsmooth_CPM_MAF5_FDR5_1MB/results/003_inputs/filtered_qsmooth_norm_cpm/cd4_atac_processed_peaks_coordinates.bed"
)$peak_name

peak_coords <- peak_coords[peak_coords$X4 %in% peaks_names, ]
colnames(peak_coords) <- c(
  "peak_chr", "peak_start", "peak_end", "peak_name",
  "X5", "X6", "X7", "X8", "X9", "X10", "X11"
)

peaks_forward <- findr_res_sig$phenotype_A[findr_res_sig$type == "forward"]
peaks_rev     <- findr_res_sig$phenotype_B[findr_res_sig$type == "reverse"]
all_peaks     <- union(peaks_forward, peaks_rev)

peak_coords <- peak_coords[peak_coords$peak_name %in% all_peaks, ]

# unified peak identifier in findr_res_sig
findr_res_sig <- findr_res_sig %>%
  mutate(
    peak_id = if_else(
      type == "forward",
      phenotype_A,   # chromatin peak in forward model
      phenotype_B    # chromatin peak in reverse model
    )
  )

# join peak coordinates back into findr_res_sig
peak_coords_min <- peak_coords %>%
  dplyr::select(peak_name, peak_chr, peak_start, peak_end)

findr_res_sig <- findr_res_sig %>%
  left_join(peak_coords_min, by = c("peak_id" = "peak_name"))

# ============================================================
# 2) Add gene coordinates to findr_res_sig (Gencode v44 GTF)
# ============================================================

gtf_path <- "~/resources/genome/hg38_gencode_raw/gencode.v44.primary_assembly.annotation.gtf"

gtf  <- import(gtf_path)
txdb <- txdbmaker::makeTxDbFromGFF(gtf_path)

gene_annot <- gtf[gtf$type == "gene"]

gene_map <- as_tibble(mcols(gene_annot)) %>%
  transmute(
    gene_id_full = gene_id,
    gene_id      = sub("\\..*$", "", gene_id),
    gene_name    = gene_name
  ) %>%
  distinct()

gene_gr <- genes(txdb)

gene_coords <- as_tibble(gene_gr) %>%
  transmute(
    gene_id     = sub("\\..*$", "", gene_id),
    gene_chr    = as.character(seqnames),
    gene_start  = start,
    gene_end    = end,
    gene_strand = as.character(strand),
    gene_tss    = if_else(gene_strand == "+", gene_start, gene_end)
  )

gene_table <- gene_map %>%
  inner_join(gene_coords, by = "gene_id")

# unified gene identifier in findr_res_sig
findr_res_sig <- findr_res_sig %>%
  mutate(
    gene_symbol = if_else(
      type == "forward",
      phenotype_B,
      phenotype_A
    )
  )

gene_table <- gene_table[gene_table$gene_name %in% findr_res_sig$gene_symbol, ]

# collapse multiple transcripts per gene_name+chr
genes_collapsed <- gene_table %>%
  group_by(gene_name, gene_chr) %>%
  summarise(
    gene_start = min(gene_start, na.rm = TRUE),
    gene_end   = max(gene_end,   na.rm = TRUE),
    .idx       = which.min(replace(gene_start, is.na(gene_start), Inf)),
    gene_strand = gene_strand[.idx],
    gene_tss    = gene_tss[.idx],
    .groups = "drop"
  ) %>%
  dplyr::select(-.idx)

findr_res_sig <- findr_res_sig %>%
  left_join(
    genes_collapsed,
    by = c(
      "gene_symbol" = "gene_name",
      "peak_chr"    = "gene_chr"
    )
  )

# ============================================================
# 3) Relative peak position features + locus model + triad ID
# ============================================================

findr_res_sig <- findr_res_sig %>%
  mutate(
    peak_mid = (peak_start + peak_end) / 2,
    signed_dist_TSS = case_when(
      is.na(gene_tss)    ~ NA_real_,
      gene_strand == "+" ~ as.numeric(peak_mid - gene_tss),
      gene_strand == "-" ~ as.numeric(gene_tss - peak_mid),
      TRUE               ~ NA_real_
    ),
    abs_dist_TSS = abs(signed_dist_TSS),
    
    peak_in_gene_body =
      !is.na(gene_start) & !is.na(gene_end) &
      peak_mid >= gene_start &
      peak_mid <= gene_end,
    
    peak_in_promoter =
      !is.na(gene_tss) &
      peak_mid >= (gene_tss - 2000L) &
      peak_mid <= (gene_tss + 2000L),
    
    rel_pos = case_when(
      peak_in_promoter    ~ "promoter",
      peak_in_gene_body   ~ "intragenic",
      signed_dist_TSS < 0 ~ "upstream",
      signed_dist_TSS > 0 ~ "downstream",
      TRUE                ~ "unknown"
    )
  )

# triad string (consistent forward/reverse orientation)
findr_res_sig <- findr_res_sig %>%
  mutate(
    triad = ifelse(
      type == "forward",
      paste0(locus_id, "_", phenotype_A, "_", phenotype_B),
      paste0(locus_id, "_", phenotype_B, "_", phenotype_A)
    )
  )

# locus_model (based on counts of forward vs reverse)
locus_models <- findr_res_sig %>%
  count(locus_id, type) %>%
  pivot_wider(names_from = type, values_from = n, values_fill = 0) %>%
  mutate(
    locus_model = case_when(
      forward > reverse ~ "var->chr->gene",
      reverse > forward ~ "var->gene->chr",
      TRUE              ~ "ambiguous"
    )
  )

findr_res_sig <- findr_res_sig %>%
  left_join(locus_models, by = "locus_id")

# ============================================================
# 4) Simple promoter Fisher test (forward vs reverse)
# ============================================================

tab_prom <- findr_res_sig %>%
  distinct(triad, type, peak_in_promoter) %>%
  count(type, peak_in_promoter) %>%
  pivot_wider(names_from = peak_in_promoter, values_from = n, values_fill = 0)

prom_mat <- as.matrix(tab_prom[, -1])
rownames(prom_mat) <- tab_prom$type
prom_mat <- prom_mat[, c("TRUE", "FALSE")]

print(fisher.test(prom_mat))

# ============================================================
# 5) Relative position enrichment (GLM with cluster SE by locus)
# ============================================================

rel_levels <- c("downstream", "intragenic", "promoter", "upstream")

relpos_enrich <- map_dfr(rel_levels, function(pos) {
  
  df_pos <- findr_res_sig %>%
    filter(rel_pos != "unknown") %>%
    mutate(is_pos = rel_pos == pos) %>%
    distinct(triad, locus_id, locus_model, is_pos) %>%
    filter(locus_model %in% c("var->chr->gene", "var->gene->chr")) %>%
    mutate(
      locus_model = factor(locus_model),
      locus_model = relevel(locus_model, ref = "var->gene->chr"),
      is_pos_num  = as.numeric(is_pos)
    )
  
  if (n_distinct(df_pos$locus_model) < 2) {
    return(tibble(
      rel_pos    = pos,
      odds_ratio = NA_real_,
      conf_low   = NA_real_,
      conf_high  = NA_real_,
      p_value    = NA_real_
    ))
  }
  
  fit <- glm(is_pos_num ~ locus_model, data = df_pos, family = binomial())
  
  vcov_locus <- vcovCL(fit, cluster = ~ locus_id)
  ct <- coeftest(fit, vcov = vcov_locus)
  
  beta <- ct["locus_modelvar->chr->gene", "Estimate"]
  se   <- ct["locus_modelvar->chr->gene", "Std. Error"]
  p    <- ct["locus_modelvar->chr->gene", "Pr(>|z|)"]
  
  OR <- exp(beta)
  CI <- exp(beta + c(-1, 1) * 1.96 * se)
  
  tibble(
    rel_pos    = pos,
    odds_ratio = OR,
    conf_low   = CI[1],
    conf_high  = CI[2],
    p_value    = p
  )
}) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  mutate(
    sig = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ ""
    ),
    or_ci_lab = ifelse(
      is.na(odds_ratio),
      NA_character_,
      sprintf("OR=%.2f (%.2f–%.2f)%s", odds_ratio, conf_low, conf_high, sig)
    )
  )

gg_relpos <- ggplot(relpos_enrich, aes(x = rel_pos, y = odds_ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_pointrange(aes(ymin = conf_low, ymax = conf_high),
                  position = position_dodge(width = 0.4),
                  size = 1.2) +
  geom_text(aes(y = conf_high * 1.15, label = or_ci_lab),
            hjust = 0, size = 3.4, na.rm = TRUE) +
  scale_y_log10() +
  labs(
    x = "Chromatin peak position",
    y = "Odds ratio (var→chr→gene vs var→gene→chr, log scale)",
    title = "Enrichment of positional classes in var→chr→gene loci"
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(5.5, 60, 5.5, 5.5)) +
  coord_flip(clip = "off")

print(gg_relpos)

write_csv(relpos_enrich, "~/cd4_qtl_paper_figures/figure_3/data/findr_mediation_relative_pos_locus_enrich.csv")

ggsave(
  "~/cd4_qtl_paper_figures/figure_3/plots/dec2025/glm_test_model_vs_rel_peak_position_v3.pdf",
  gg_relpos, width = 14, height = 12, units = "cm", dpi = 300
)

# ============================================================
# 6) Category enrichment function (GLM with cluster SE by locus)
# ============================================================

run_enrich_by_category <- function(df,
                                   category_col,
                                   levels = NULL,
                                   unknown_values = c("unknown", NA),
                                   model_keep = c("var->chr->gene", "var->gene->chr")) {
  
  if (is.null(levels)) {
    levels <- df %>%
      pull(.data[[category_col]]) %>%
      unique() %>%
      sort() %>%
      as.character()
  }
  
  df0 <- df %>%
    filter(!is.na(.data[[category_col]])) %>%
    filter(!(.data[[category_col]] %in% setdiff(unknown_values, NA))) %>%
    filter(locus_model %in% model_keep)
  
  map_dfr(levels, function(cat) {
    
    df_cat <- df0 %>%
      mutate(is_cat = .data[[category_col]] == cat) %>%
      distinct(triad, locus_id, locus_model, is_cat) %>%
      mutate(
        locus_model = factor(locus_model),
        locus_model = relevel(locus_model, ref = "var->gene->chr"),
        is_cat_num  = as.numeric(is_cat)
      )
    
    if (n_distinct(df_cat$locus_model) < 2) {
      return(tibble(
        category   = cat,
        odds_ratio = NA_real_,
        conf_low   = NA_real_,
        conf_high  = NA_real_,
        p_value    = NA_real_
      ))
    }
    
    fit <- glm(is_cat_num ~ locus_model, data = df_cat, family = binomial())
    
    vc <- vcovCL(fit, cluster = ~ locus_id)
    ct <- coeftest(fit, vcov. = vc)
    
    term <- "locus_modelvar->chr->gene"
    beta <- ct[term, "Estimate"]
    se   <- ct[term, "Std. Error"]
    p    <- ct[term, "Pr(>|z|)"]
    
    OR <- exp(beta)
    CI <- exp(beta + c(-1, 1) * 1.96 * se)
    
    tibble(
      category   = cat,
      odds_ratio = OR,
      conf_low   = CI[1],
      conf_high  = CI[2],
      p_value    = p
    )
  })
}

caqtl_enrich <- run_enrich_by_category(
  df = findr_res_sig,
  category_col = "caqtl_category"
) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  mutate(
    sig = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ ""
    ),
    or_ci_lab = ifelse(
      is.na(odds_ratio),
      NA_character_,
      sprintf("OR=%.2f (%.2f–%.2f)%s", odds_ratio, conf_low, conf_high, sig)
    )
  )

gg_caqtl <- ggplot(caqtl_enrich, aes(x = category, y = odds_ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_pointrange(aes(ymin = conf_low, ymax = conf_high), size = 1.2) +
  geom_text(aes(y = conf_high * 1.15, label = or_ci_lab),
            hjust = 0, size = 3.4, na.rm = TRUE) +
  scale_y_log10() +
  coord_flip(clip = "off") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(5.5, 60, 5.5, 5.5)) +
  labs(
    x = "caQTL category",
    y = "Odds ratio (var→chr→gene vs var→gene->chr, log scale)",
    title = "Enrichment of caQTL categories in var→chr→gene loci"
  )

print(gg_caqtl)
write_csv(caqtl_enrich, "~/cd4_qtl_paper_figures/figure_3/data/findr_mediation_caqtl_category_locus_enrich.csv")

ggsave(
  "~/cd4_qtl_paper_figures/figure_3/plots/dec2025/glm_test_model_vs_caqtl-categry_v3.pdf",
  gg_caqtl, width = 14, height = 12, units = "cm", dpi = 300
)

# ============================================================
# 7) Locus-based ChromBPNet enrichment (GLM)
# ============================================================

# Load ChromBPNet peak sets
upset_memberships <- readRDS(
  "~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/upset_credible_sets.rds"
) %>%
  filter(
    coloc_status   == "coloc",
    chromBPNet     == "TRUE",
    caqtl_category != "no_Peak_overlap"
  ) %>%
  mutate(
    peak_name = str_extract(variant_id, "^cd4_atac_summits_peak_[0-9]+[a-z]*")
  )

cpbnet_peaks <- upset_memberships %>% pull(peak_name) %>% unique()
cpbnet_peaks_TF <- upset_memberships %>%
  filter(chromBPNet_motif == TRUE) %>%
  pull(peak_name) %>%
  unique()

# locus-level dataset: one row per (locus_id, type)
locus_df <- findr_res_sig %>%
  filter(!is.na(peak_id), nzchar(peak_id)) %>%
  mutate(
    type = factor(type, levels = c("reverse", "forward")),
    has_cpnet    = peak_id %in% cpbnet_peaks,
    has_cpnet_tf = peak_id %in% cpbnet_peaks_TF
  ) %>%
  group_by(locus_id, type) %>%
  summarise(
    any_cpnet    = any(has_cpnet),
    any_cpnet_tf = any(has_cpnet_tf),
    n_triads     = n(),
    n_peaks      = n_distinct(peak_id),
    .groups = "drop"
  )

fit_any <- glm(any_cpnet ~ type, data = locus_df, family = binomial())
fit_tf  <- glm(any_cpnet_tf ~ type, data = locus_df, family = binomial())

tidy_or <- function(fit, label) {
  broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "typeforward") %>%
    transmute(
      feature    = label,
      odds_ratio = estimate,
      conf_low   = conf.low,
      conf_high  = conf.high,
      p_value    = p.value
    )
}

glm_locus_enrich <- bind_rows(
  tidy_or(fit_any, "ChromBPNet var (any)"),
  tidy_or(fit_tf,  "ChromBPNet var + TF motif")
) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ ""
    ),
    or_ci_lab = sprintf("OR=%.2f (%.2f–%.2f)%s", odds_ratio, conf_low, conf_high, sig)
  )

print(glm_locus_enrich)

gg_cbp_locus <- ggplot(glm_locus_enrich, aes(x = feature, y = odds_ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_pointrange(aes(ymin = conf_low, ymax = conf_high), size = 1.2) +
  geom_text(aes(y = conf_high * 1.15, label = or_ci_lab),
            hjust = 0, size = 3.2, na.rm = TRUE) +
  scale_y_log10() +
  labs(
    x = NULL,
    y = "Odds ratio (forward vs reverse, log scale)",
    title = "Locus-based enrichment of ChromBPNet peaks in forward vs reverse loci"
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(5.5, 70, 5.5, 5.5)) +
  coord_flip(clip = "off")

print(gg_cbp_locus)

ggsave(
  "~/cd4_qtl_paper_figures/figure_3/plots/dec2025/glm_locus_enrich_model_vs_cbpnet.pdf",
  gg_cbp_locus,
  width = 14, height = 12, units = "cm", dpi = 300
)

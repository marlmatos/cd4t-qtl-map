# ============================================================
# Chromatin state summaries (eQTL vs caQTL chi-square)
# + caQTL category enrichment (Fisher OR forest)
# + 75% CS-size filtering
# Author: Marliette Matos
# ============================================================

.libPaths(c("/gchm/R/x86_64-pc-linux-gnu-library/4.4",
            "/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(forcats)
  library(ggplot2)
  library(scales)
  library(purrr)
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(readr)
})

# -----------------------------
# 0) Paths / inputs
# -----------------------------
summary_df      <- "/gchm/cd4_qtl_paper_figures/figure_1/data/eqtl_caqtl_finemapping_coloc_all.tsv"
chrom_bed_path  <- "~/cd4_qtl_paper_figures/resources/ENOCDE_chromatin_state/ENCFF651QPF.bed.gz"

# what defines eQTL / caQTL in your table
eqtl_only  <- c("coloc", "eQTL_only")
caqtl_only <- c("coloc", "caQTL_only")

# CS-size filter: keep bottom q (e.g., 0.75)
q_keep <- 0.75

# caQTL category collapse rule
caqtl_cat_levels <- c("no_Peak_overlap", "in_other_Peak", "in_caPeak")

# chromatin state groups (collapsed labels) + plotting order
cs_state_levels <- c(
  "Promoter",
  "Enhancer_active",
  "Enhancer_weak",
  "Transcribed",
  "Repressed_polycomb",
  "ZNF_Het_Quies"
)

# -----------------------------
# 1) Read data
# -----------------------------
dt <- data.table::fread(summary_df)

eqtl_df0  <- dt %>% dplyr::filter(coloc_status %in% eqtl_only)
caqtl_df0 <- dt %>% dplyr::filter(coloc_status %in% caqtl_only)
rm(dt)

# -----------------------------
# 2) Plot helper (dodged + per-state q label)
# -----------------------------
plot_dodged_states_with_q <- function(plot_df,
                                      group_col, state_col, prop_col,
                                      q_tbl,
                                      title = "", subtitle = "",
                                      y_pad = 0.06) {
  ann <- plot_df %>%
    group_by(.data[[state_col]]) %>%
    summarise(max_prop = max(.data[[prop_col]]), .groups = "drop") %>%
    transmute(cs_state = .data[[state_col]], max_prop) %>%
    left_join(q_tbl %>% transmute(cs_state, q_label), by = "cs_state") %>%
    mutate(y = pmin(max_prop + y_pad, 0.98))

  ymax <- min(1, max(ann$y, na.rm = TRUE) + 0.05)

  ggplot(plot_df, aes(x = .data[[state_col]], y = .data[[prop_col]], fill = .data[[group_col]])) +
    geom_col(position = position_dodge(width = 0.85), width = 0.75,
             color = "grey20", linewidth = 0.2) +
    geom_text(
      aes(label = scales::percent(.data[[prop_col]], accuracy = 0.1)),
      position = position_dodge(width = 0.85),
      vjust = -0.25, size = 3
    ) +
    geom_text(
      data = ann,
      aes(x = cs_state, y = y, label = q_label),
      inherit.aes = FALSE,
      size = 3
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, ymax),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = "Chromatin state (collapsed)", y = "Proportion (within group)", fill = NULL,
         title = title, subtitle = subtitle) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

# -----------------------------
# 3) Core helpers: filtering + GRanges + annotation + majority vote
# -----------------------------
filter_bottom_cs <- function(df, cs_col, var_col, q = 0.75) {
  cs_sym  <- rlang::sym(cs_col)
  var_sym <- rlang::sym(var_col)

  cs_sizes <- df %>%
    distinct(!!cs_sym, !!var_sym) %>%
    dplyr::count(!!cs_sym, name = "n_variants_per_cs")

  q_cut <- as.numeric(stats::quantile(cs_sizes$n_variants_per_cs, q, na.rm = TRUE))

  keep_cs <- cs_sizes %>%
    dplyr::filter(n_variants_per_cs <= q_cut) %>%
    dplyr::select(!!cs_sym)

  df_filt <- df %>%
    dplyr::semi_join(keep_cs, by = cs_col)

  list(df = df_filt, cs_sizes = cs_sizes, q_cut = q_cut)
}

make_variants_gr <- function(df, var_id_col, chr_col, pos_col) {
  u <- df %>%
    distinct(
      variant_id = .data[[var_id_col]],
      chr        = .data[[chr_col]],
      pos        = .data[[pos_col]]
    ) %>%
    filter(!is.na(variant_id), !is.na(chr), !is.na(pos))

  gr <- GRanges(seqnames = u$chr,
                ranges   = IRanges(start = u$pos, end = u$pos))
  mcols(gr)$variant_id <- u$variant_id
  gr
}

collapse_chrom_state <- function(x) {
  case_when(
    x %in% c("TssA", "TssFlnkU", "TssFlnkD", "TssFlnk", "TssBiv") ~ "Promoter",
    x %in% c("EnhA1", "EnhA2", "EnhG1", "EnhG2", "EnhBiv")       ~ "Enhancer_active",
    x %in% c("EnhWk")                                            ~ "Enhancer_weak",
    x %in% c("Tx", "TxWk")                                       ~ "Transcribed",
    x %in% c("ReprPC", "ReprPCWk")                               ~ "Repressed_polycomb",
    x %in% c("ZNF/Rpts", "Het", "Quies")                         ~ "ZNF_Het_Quies",
    TRUE                                                         ~ NA_character_
  )
}

annotate_variants_with_chromstate <- function(variantsGR, chrom_bed_path) {
  states <- rtracklayer::import.bed(chrom_bed_path)
  ov <- findOverlaps(variantsGR, states, ignore.strand = TRUE)

  variantsGR$T_cell_chromatin_state <- NA_character_
  variantsGR$T_cell_chromatin_state[queryHits(ov)] <- states$name[subjectHits(ov)]

  tibble(
    variant_id = mcols(variantsGR)$variant_id,
    T_cell_chromatin_state = variantsGR$T_cell_chromatin_state,
    chrom_state_group = collapse_chrom_state(variantsGR$T_cell_chromatin_state)
  )
}

# majority vote per CS (optionally stratified by some grouping columns like caqtl_category)
majority_vote_cs <- function(df, cs_col, label_col, by_cols = NULL, levels_ordered = NULL) {
  cs_sym  <- rlang::sym(cs_col)
  lab_sym <- rlang::sym(label_col)

  tmp <- df %>%
    dplyr::mutate(.cs = !!cs_sym, .lab = !!lab_sym)

  if (!is.null(levels_ordered)) {
    tmp <- tmp %>%
      dplyr::mutate(.lab = factor(.lab, levels = levels_ordered, ordered = TRUE))
  }

  if (is.null(by_cols)) {
    out <- tmp %>%
      dplyr::count(.cs, .lab, name = "N") %>%
      dplyr::arrange(.cs, dplyr::desc(N), .lab) %>%
      dplyr::group_by(.cs) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(!!cs_sym := .cs, cs_label = as.character(.lab))
  } else {
    by_syms <- rlang::syms(by_cols)
    out <- tmp %>%
      dplyr::count(!!!by_syms, .cs, .lab, name = "N") %>%
      dplyr::arrange(!!!by_syms, .cs, dplyr::desc(N), .lab) %>%
      dplyr::group_by(!!!by_syms, .cs) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(!!!by_syms, !!cs_sym := .cs, cs_label = as.character(.lab))
  }

  out
}

# -----------------------------
# 4) Chi-square helpers (chromatin state distributions)
# -----------------------------
chisq_global <- function(count_df, group_col, state_col, count_col) {
  tab <- count_df %>%
    transmute(
      g = as.character(.data[[group_col]]),
      s = as.character(.data[[state_col]]),
      n = .data[[count_col]]
    ) %>%
    pivot_wider(names_from = s, values_from = n, values_fill = 0) %>%
    arrange(g)

  mat <- as.matrix(tab[, -1, drop = FALSE])
  rownames(mat) <- tab$g
  suppressWarnings(chisq.test(mat))
}

# per-state 2x2 (state vs rest), BH across states
per_state_2x2 <- function(count_df, group_col, state_col, count_col) {
  tab <- count_df %>%
    transmute(
      g = as.character(.data[[group_col]]),
      s = as.character(.data[[state_col]]),
      n = .data[[count_col]]
    ) %>%
    pivot_wider(names_from = s, values_from = n, values_fill = 0) %>%
    arrange(g)

  mat <- as.matrix(tab[, -1, drop = FALSE])
  rownames(mat) <- tab$g
  states_here <- colnames(mat)

  out <- map_dfr(states_here, function(st) {
    in_st <- mat[, st]
    other <- rowSums(mat) - in_st
    m2 <- cbind(in_state = in_st, other = other)
    tt <- suppressWarnings(chisq.test(m2))

    tibble(
      cs_state = st,
      p_value  = tt$p.value,
      method   = tt$method
    )
  }) %>%
    mutate(
      p_adj_BH = p.adjust(p_value, method = "BH"),
      q_label  = paste0("q = ", format.pval(p_adj_BH, digits = 2, eps = 1e-300))
    ) %>%
    arrange(p_value)

  out
}

# ============================================================
# 5) Build CS-level chromatin-state labels for eQTL + caQTL
# ============================================================

# ---- eQTL: filter CS sizes -> annotate variants -> CS majority
eqtl_f <- filter_bottom_cs(eqtl_df0,
                           cs_col  = "finemapped_cs_eqtl",
                           var_col = "variant_id.y",
                           q = q_keep)
eqtl_df_filt <- eqtl_f$df
message("eQTL CS-size cutoff (", q_keep*100, "%): <= ", eqtl_f$q_cut, " variants per CS")

eqtl.variantsGR <- make_variants_gr(eqtl_df_filt,
                                    var_id_col = "variant_id.y",
                                    chr_col    = "chr.y",
                                    pos_col    = "variant_pos.y")

eqtl_var_state_tbl <- annotate_variants_with_chromstate(eqtl.variantsGR, chrom_bed_path)

eqtl_state_df <- eqtl_df_filt %>%
  dplyr::select(finemapped_cs_eqtl, variant_id = variant_id.y) %>%
  left_join(eqtl_var_state_tbl %>% dplyr::select(variant_id, chrom_state_group),
            by = "variant_id") %>%
  distinct(finemapped_cs_eqtl, variant_id, chrom_state_group) %>%
  mutate(chrom_state_group = factor(chrom_state_group, levels = cs_state_levels, ordered = TRUE))

eqtl_state_major <- majority_vote_cs(
  eqtl_state_df,
  cs_col = "finemapped_cs_eqtl",
  label_col = "chrom_state_group",
  by_cols = NULL,
  levels_ordered = cs_state_levels
) %>%
  dplyr::rename(cs_state = cs_label)

# ---- caQTL: normalize categories -> filter CS sizes -> annotate -> CS majority by category
caqtl_df0 <- caqtl_df0 %>%
  mutate(
    caqtl_category = ifelse(cs_type_any_var %in% c("in_corr_Peak", "in_other_Peak"),
                            "in_other_Peak",
                            cs_type_any_var),
    caqtl_category = factor(caqtl_category, levels = caqtl_cat_levels)
  )

caqtl_f <- filter_bottom_cs(caqtl_df0,
                            cs_col  = "finemapped_cs_caqtl",
                            var_col = "variant_id.x",
                            q = q_keep)
caqtl_df_filt <- caqtl_f$df
message("caQTL CS-size cutoff (", q_keep*100, "%): <= ", caqtl_f$q_cut, " variants per CS")

caqtl.variantsGR <- make_variants_gr(caqtl_df_filt,
                                     var_id_col = "variant_id.x",
                                     chr_col    = "chr.x",
                                     pos_col    = "variant_pos.x")

caqtl_var_state_tbl <- annotate_variants_with_chromstate(caqtl.variantsGR, chrom_bed_path)

caqtl_state_df <- caqtl_df_filt %>%
  dplyr::select(finemapped_cs_caqtl, caqtl_category, variant_id = variant_id.x) %>%
  left_join(caqtl_var_state_tbl %>% dplyr::select(variant_id, chrom_state_group),
            by = "variant_id") %>%
  distinct(finemapped_cs_caqtl, caqtl_category, variant_id, chrom_state_group) %>%
  mutate(chrom_state_group = factor(chrom_state_group, levels = cs_state_levels, ordered = TRUE))

caqtl_state_major <- majority_vote_cs(caqtl_state_df,
                                      cs_col = "finemapped_cs_caqtl",
                                      label_col = "chrom_state_group",
                                      by_cols = c("caqtl_category"),
                                      levels_ordered = cs_state_levels) %>%
  dplyr::rename(cs_state = cs_label)

# ============================================================
# 6) eQTL vs caQTL_all chi-square + plot (chromatin states)
# ============================================================

eqtl_cs_state_tbl <- eqtl_state_major %>%
  transmute(group = factor("eQTL", levels = c("eQTL","caQTL_all")),
            cs_state = factor(cs_state, levels = cs_state_levels))

caqtl_all_cs_state_tbl <- caqtl_state_major %>%
  transmute(group = factor("caQTL_all", levels = c("eQTL","caQTL_all")),
            cs_state = factor(cs_state, levels = cs_state_levels))

eqtl_vs_caqtl_state_counts <- bind_rows(eqtl_cs_state_tbl, caqtl_all_cs_state_tbl) %>%
  count(group, cs_state, name = "n_cs") %>%
  complete(group = factor(c("eQTL","caQTL_all"), levels = c("eQTL","caQTL_all")),
           cs_state = factor(cs_state_levels, levels = cs_state_levels),
           fill = list(n_cs = 0L)) %>%
  group_by(group) %>%
  mutate(prop = n_cs / sum(n_cs)) %>%
  ungroup()

chi_eqtl_vs_caqtl_global <- chisq_global(eqtl_vs_caqtl_state_counts,
                                         group_col = "group",
                                         state_col = "cs_state",
                                         count_col = "n_cs")

chi_eqtl_vs_caqtl_perstate <- per_state_2x2(eqtl_vs_caqtl_state_counts,
                                            group_col = "group",
                                            state_col = "cs_state",
                                            count_col = "n_cs") %>%
  mutate(cs_state = factor(cs_state, levels = cs_state_levels))

plot_eqtl_vs_caqtl_states <- plot_dodged_states_with_q(
  plot_df = eqtl_vs_caqtl_state_counts %>%
    mutate(cs_state = factor(cs_state, levels = cs_state_levels)),
  group_col = "group",
  state_col = "cs_state",
  prop_col  = "prop",
  q_tbl     = chi_eqtl_vs_caqtl_perstate,
  title     = "Chromatin state distribution: eQTL vs caQTL (collapsed)",
  subtitle  = "State labels show BH-adjusted q-values from per-state (state vs rest) 2x2 chi-square."
)

# ============================================================
# 7) caQTL category enrichment (Fisher OR forest)
#    (uses CS-level majority chromatin state per category)
# ============================================================

# one row per CS (already CS-majority), so distinct is optional but safe
df_enrich <- caqtl_state_major %>%
  distinct(caqtl_category, finemapped_cs_caqtl, cs_state)

N <- nrow(df_enrich)

cat_tot   <- df_enrich %>% count(caqtl_category, name = "n_cat")
state_tot <- df_enrich %>% count(cs_state, name = "n_state")

grid <- df_enrich %>%
  count(caqtl_category, cs_state, name = "a") %>%
  left_join(cat_tot, by = "caqtl_category") %>%
  left_join(state_tot, by = "cs_state") %>%
  mutate(
    b = n_cat - a,
    c = n_state - a,
    d = N - a - b - c
  )

fisher_one <- function(a,b,c,d){
  ft <- fisher.test(matrix(c(a,b,c,d), nrow = 2, byrow = TRUE))
  tibble(
    odds_ratio = unname(ft$estimate),
    conf_low   = unname(ft$conf.int[1]),
    conf_high  = unname(ft$conf.int[2]),
    p_value    = ft$p.value
  )
}

res <- grid %>%
  mutate(out = pmap(list(a,b,c,d), fisher_one)) %>%
  unnest(out) %>%
  group_by(caqtl_category) %>%
  mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(caqtl_category, q_value)

# helper: q-value -> asterisks
q_to_stars <- function(q) {
  dplyr::case_when(
    is.na(q) ~ "",
    q < 1e-3 ~ "***",
    q < 1e-2 ~ "**",
    q < 5e-2 ~ "*",
    TRUE     ~ ""
  )
}

plot_df <- res %>%
  mutate(
    sig   = q_value < 0.05,
    stars = q_to_stars(q_value),
    cs_state = factor(cs_state, levels = cs_state_levels)
  ) %>%
  mutate(
    caqtl_category = factor(
      caqtl_category,
      levels = c("in_caPeak","in_other_Peak","no_Peak_overlap")
    )
  )

# put stars just to the right of the CI, but not past the axis max
star_pos <- plot_df %>%
  group_by(caqtl_category) %>%
  summarise(
    x_star = pmin(6.0, max(conf_high, na.rm = TRUE) * 1.05),
    .groups = "drop"
  )

plot_df <- plot_df %>% left_join(star_pos, by = "caqtl_category")

p_forest <- ggplot(plot_df, aes(y = cs_state, x = odds_ratio)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high),
                 height = 0.18, linewidth = 0.6) +
  geom_point(aes(size = -log10(q_value), alpha = sig), shape = 16) +
  geom_text(aes(x = x_star, label = stars),
            vjust = 0.35, size = 4, na.rm = TRUE) +
  scale_x_log10(
    limits = c(0.1, 6.0),
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5),
    labels = c("0.1","0.2","0.5","1","2","5")
  ) +
  coord_cartesian(clip = "off") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
  scale_size_continuous(range = c(1.6, 4.8), name = expression(-log[10](q))) +
  facet_wrap(~ caqtl_category, ncol = 1) +
  labs(
    x = "Odds ratio (log scale)", y = NULL,
    title = "Chromatin-state enrichment by caQTL category",
    subtitle = "Points show OR; bars show 95% CI. Asterisks denote BH q-value (*<0.05, **<0.01, ***<0.001)."
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text = element_text(size = 10),
    panel.spacing = unit(0.9, "lines"),
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  )

# ============================================================
# 8) Outputs 
# ============================================================

# Print key test objects/tables
chi_eqtl_vs_caqtl_global
chi_eqtl_vs_caqtl_perstate
res

# Show plots 
plot_eqtl_vs_caqtl_states
p_forest

# Save 
out_dir <- "~/cd4_qtl_paper_figures/figure_1/plotting/plots_oct2025/chromstate_minimal"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write_tsv_safe <- function(x, path) {
  x2 <- tryCatch(as.data.frame(x), error = function(e) x)
  readr::write_tsv(x2, path)
}

save_plot_both <- function(p, stem, width = 8, height = 5, dpi = 300) {
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = width, height = height, dpi = dpi)
}

# save tables
saveRDS(chi_eqtl_vs_caqtl_global, file.path(out_dir, "chi_eqtl_vs_caqtl_global.rds"))
write_tsv_safe(
  tibble::tibble(
    statistic = unname(chi_eqtl_vs_caqtl_global$statistic),
    df        = unname(chi_eqtl_vs_caqtl_global$parameter),
    p_value   = unname(chi_eqtl_vs_caqtl_global$p.value),
    method    = unname(chi_eqtl_vs_caqtl_global$method)
  ),
  file.path(out_dir, "chi_eqtl_vs_caqtl_global.tsv")
)
write_tsv_safe(chi_eqtl_vs_caqtl_perstate, file.path(out_dir, "chi_eqtl_vs_caqtl_perstate.tsv"))
write_tsv_safe(res, file.path(out_dir, "caqtl_category_chromstate_enrichment_fisher.tsv"))

# save plots
save_plot_both(plot_eqtl_vs_caqtl_states, "chrom_state_dodged_eqtl_vs_caqtlall_with_q", width = 9, height = 5)
save_plot_both(p_forest, "caqtl_chromstate_enrichment_forest", width = 6, height = 8)

message("Saved outputs to: ", out_dir)
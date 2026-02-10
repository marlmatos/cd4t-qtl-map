# ======================================================================
# 4) coloc percent
# + pairwise category comparisons (coloc vs not_coloc) with BH q-values
# + stacked bar (optional) + barplot-with-brackets (main)
# Author: Marliette Matos
# ======================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(purrr)
  library(tibble)
  library(rlang)
  library(forcats)
  library(readr)
  
})


# -----------------------------
# Global style
# -----------------------------
BAR_WIDTH <- 0.6
BASE_SIZE <- 12

# -----------------------------
# Plot helpers
# -----------------------------
stack_prop_onebar <- function(df, fill_var, title, legend_title, label_min = 0.03) {
  ggplot(df, aes(x = "All CS", y = n, fill = {{ fill_var }})) +
    geom_col(width = BAR_WIDTH, color = "grey20", linewidth = 0.2, position = "fill") +
    geom_text(
      aes(label = ifelse(n / sum(n) >= label_min, percent(n / sum(n), accuracy = 0.1), "")),
      position = position_fill(vjust = 0.5),
      size = 3, color = "grey20"
    ) +
    scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = NULL, title = title, fill = legend_title) +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme_minimal(base_size = BASE_SIZE) +
    theme(axis.ticks.x = element_blank(), legend.position = "right") +
    coord_flip()
}

stack_prop_by_group <- function(df, x_var, fill_var, title, legend_title, label_min = 0.03) {
  ggplot(df, aes(x = {{ x_var }}, y = prop, fill = {{ fill_var }})) +
    geom_col(width = BAR_WIDTH, color = "grey20", linewidth = 0.2, position = "fill") +
    geom_text(
      aes(label = ifelse(prop >= label_min, percent(prop, accuracy = 0.1), "")),
      position = position_fill(vjust = 0.5),
      size = 3, color = "grey20"
    ) +
    scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = NULL, title = title, fill = legend_title) +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme_minimal(base_size = BASE_SIZE) +
    theme(axis.ticks.x = element_blank(), legend.position = "right") +
    coord_flip()
}

# -----------------------------
# Stats helper: ONE comparison (global 2xK) for coloc vs not-coloc across categories
# -----------------------------
coloc_only_global <- function(plot_df,
                              cat_col = "caqtl_category",
                              status_col = "coloc_status",
                              n_col = "n",
                              coloc_level = "coloc") {
  
  coloc_counts <- plot_df %>%
    dplyr::mutate(is_coloc = (.data[[status_col]] == coloc_level)) %>%
    dplyr::group_by(.data[[cat_col]], is_coloc) %>%
    dplyr::summarise(n = sum(.data[[n_col]]), .groups = "drop") %>%
    tidyr::complete(!!rlang::sym(cat_col), is_coloc, fill = list(n = 0L)) %>%
    tidyr::pivot_wider(names_from = is_coloc, values_from = n, values_fill = 0) %>%
    dplyr::rename(n_not_coloc = `FALSE`, n_coloc = `TRUE`) %>%
    dplyr::mutate(
      total = n_coloc + n_not_coloc,
      prop_coloc = dplyr::if_else(total > 0, n_coloc / total, NA_real_)
    )
  
  # 2 x K contingency: rows = categories, cols = (coloc, not_coloc)
  mat <- as.matrix(coloc_counts[, c("n_coloc", "n_not_coloc"), drop = FALSE])
  rownames(mat) <- as.character(coloc_counts[[cat_col]])
  storage.mode(mat) <- "integer"
  
  # One global test: chi-square if expected ok, else Fisher
  chisq0 <- suppressWarnings(stats::chisq.test(mat, correct = FALSE))
  exp_ok <- all(chisq0$expected >= 5)
  
  global <- if (exp_ok) chisq0 else stats::fisher.test(mat)
  
  list(
    coloc_counts = coloc_counts,
    contingency  = mat,
    global_test  = global
  )
}

# -----------------------------
# Inputs
# -----------------------------
summary_df <- "/gchm/cd4_qtl_paper_figures/figure_1/data/eqtl_caqtl_finemapping_coloc_all.tsv"

eqtl_only  <- c("coloc", "eQTL_only")
caqtl_only <- c("coloc", "caQTL_only")

# Palette
qtl_cmap_coloc2 <- c(
  "eQTL_only"  = "#c7197c",
  "caQTL_only" = "#9ccb86",
  "coloc"      = "#f39c12"
)

# -----------------------------
# Read and subset
# -----------------------------
dt <- data.table::fread(summary_df)
eqtl_df0  <- dt %>% filter(coloc_status %in% eqtl_only)
caqtl_df0 <- dt %>% filter(coloc_status %in% caqtl_only)
rm(dt)

# ======================================================================
# A) eQTL coloc percent (one stacked bar)
# ======================================================================
plot_df_eqtl <- eqtl_df0 %>%
  distinct(finemapped_cs_eqtl, coloc_status) %>%
  dplyr::count(coloc_status, name = "n") %>%
  mutate(
    coloc_status = factor(coloc_status, levels = c("coloc", "eQTL_only")),
    prop  = n / sum(n),
    label = paste0(n, " (", percent(prop, accuracy = 0.1), ")")
  )

eqtl_coloc_pt <- stack_prop_onebar(
  df           = plot_df_eqtl,
  fill_var     = coloc_status,
  title        = "% eQTL coloc",
  legend_title = "Coloc status"
) +
  geom_text(
    data = plot_df_eqtl,
    aes(label = label),
    position = position_fill(vjust = 0.5),
    size = 3
  ) +
  scale_fill_manual(
    values = qtl_cmap_coloc2[c("coloc", "eQTL_only")],
    drop   = FALSE
  )

eqtl_coloc_pt

# ======================================================================
# B) caQTL coloc percent by category (stacked bars)
# ======================================================================
desired_cats  <- c("in_caPeak","in_other_Peak","no_Peak_overlap")
status_levels <- c("coloc", "caQTL_only")

caqtl_df0 <- caqtl_df0 %>%
  mutate(
    caqtl_category = if_else(cs_type_any_var %in% c("in_other_Peak", "in_corr_Peak"),
                             "in_other_Peak", cs_type_any_var),
    caqtl_category = factor(caqtl_category, levels = desired_cats),
    coloc_status   = factor(coloc_status, levels = status_levels)
  )

plot_df_caqtl_stack <- caqtl_df0 %>%
  distinct(finemapped_cs_caqtl, coloc_status, caqtl_category) %>%
  count(caqtl_category, coloc_status, name = "n") %>%
  complete(caqtl_category, coloc_status, fill = list(n = 0L)) %>%
  group_by(caqtl_category) %>%
  mutate(
    prop  = n / sum(n),
    label = paste0(n, " (", percent(prop, accuracy = 0.1), ")")
  ) %>%
  ungroup() %>%
  mutate(caqtl_category = factor(caqtl_category, levels = rev(desired_cats)))  # flipped order for coord_flip

plot_df_caqtl_stack_lab <- plot_df_caqtl_stack %>% filter(prop >= 0.06)

caqtl_coloc_stack_pt <- stack_prop_by_group(
  df           = plot_df_caqtl_stack,
  x_var        = caqtl_category,
  fill_var     = coloc_status,
  title        = "% caQTL coloc by category",
  legend_title = "Coloc status"
) +
  geom_text(
    data = plot_df_caqtl_stack_lab,
    aes(label = label),
    position = position_fill(vjust = 0.5, reverse = TRUE),
    size = 3
  ) +
  scale_fill_manual(values = qtl_cmap_coloc2[status_levels], drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05)))

caqtl_coloc_stack_pt

# ======================================================================
# C) supplement: % coloc (bars) + pairwise brackets with BH q-values
# ======================================================================
# For testing we want plot_df with counts per (category, status)
plot_df_pairwise <- plot_df_caqtl_stack %>%
  transmute(
    caqtl_category = factor(caqtl_category, levels = desired_cats),
    coloc_status   = factor(coloc_status,   levels = status_levels),
    n              = n
  )

out <- coloc_only_global(plot_df_pairwise)

coloc_counts_tbl   <- out$coloc_counts %>%
  mutate(
    caqtl_category = factor(caqtl_category, levels = desired_cats),
    pct_label = percent(prop_coloc, accuracy = 0.1)
  ) %>%
  arrange(caqtl_category)

global_coloc_tbl <- out$global_test

gt <- out$global_test  # or just global_coloc_tbl if that object is the htest

x2 <- unname(gt$statistic)
df <- unname(gt$parameter)

p_exact <- pchisq(x2, df = df, lower.tail = FALSE)
sprintf("Pearson’s χ²(%d) = %.2f, p = %.2e",
        df, x2, p_exact)


# bracket geometry
brackets_df <- pairwise_coloc_tbl %>%
  mutate(
    g1 = factor(g1, levels = desired_cats),
    g2 = factor(g2, levels = desired_cats),
    x1 = as.numeric(g1),
    x2 = as.numeric(g2),
    span = abs(x2 - x1)
  ) %>%
  arrange(span, x1, x2)

base_y   <- max(coloc_counts_tbl$prop_coloc, na.rm = TRUE) + 0.05
step_y   <- 0.07
tick_h   <- 0.02
text_pad <- 0.012

brackets_df <- brackets_df %>%
  mutate(
    y = base_y + (row_number() - 1) * step_y,
    xmid = (x1 + x2) / 2
  )

bracket_segs <- bind_rows(
  brackets_df %>% transmute(x = x1, xend = x2, y = y, yend = y),
  brackets_df %>% transmute(x = x1, xend = x1, y = y, yend = y - tick_h),
  brackets_df %>% transmute(x = x2, xend = x2, y = y, yend = y - tick_h)
)

ymax <- max(brackets_df$y, na.rm = TRUE) + 0.08

coloc_barplot_with_brackets <- ggplot(coloc_counts_tbl,
                                      aes(x = caqtl_category, y = prop_coloc)) +
  geom_col(width = 0.7, color = "grey20", linewidth = 0.25) +
  geom_text(aes(label = pct_label), vjust = -0.4, size = 3) +
  geom_segment(
    data = bracket_segs,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    linewidth = 0.4,
    color = "grey20"
  ) +
  geom_text(
    data = brackets_df,
    aes(x = xmid, y = y + text_pad, label = q_label),
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, ymax),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "% coloc (within category)",
    title = "Coloc percent across caQTL categories",
    subtitle = "Brackets: pairwise 2×2 chi-square (coloc vs not_coloc); BH-adjusted q-values."
  ) +
  theme_classic(base_size = BASE_SIZE) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

coloc_barplot_with_brackets

# ======================================================================
# SAVE PLOTS + TABLES (PDF + PNG) for coloc percent scripts
# ======================================================================

# ---- Output directory ----
out_dir <- "~/cd4_qtl_paper_figures/figure_1/plotting/plots_oct2025/coloc_percent"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Helpers ----
write_tsv_safe <- function(x, path) {
  if (is.null(x)) return(invisible(NULL))
  x2 <- tryCatch(as.data.frame(x), error = function(e) x)
  readr::write_tsv(x2, path)
  invisible(NULL)
}

save_plot_both <- function(p, stem, width = 7, height = 4, dpi = 300) {
  if (is.null(p)) return(invisible(NULL))
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".pdf")),
                  plot = p, width = width, height = height, device = cairo_pdf)
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".png")),
                  plot = p, width = width, height = height, dpi = dpi)
  invisible(NULL)
}



# eQTL: one stacked bar (coloc vs eQTL_only)
save_plot_both(eqtl_coloc_pt, "eqtl_coloc_onebar", width = 6.5, height = 3.8)

# caQTL: stacked by category (coloc vs caQTL_only)  [optional plot]
save_plot_both(caqtl_coloc_stack_pt, "caqtl_coloc_stacked_by_category", width = 7.5, height = 4.5)

# caQTL: % coloc bars + brackets with BH q-values (MAIN)
save_plot_both(coloc_barplot_with_brackets, "caqtl_coloc_percent_with_pairwise_brackets", width = 4, height = 5)

message("Saved outputs to: ", out_dir)
print(list.files(out_dir))





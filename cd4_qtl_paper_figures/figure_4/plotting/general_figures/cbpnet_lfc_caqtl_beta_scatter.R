library(dplyr)
library(ggrastr)
library(ggplot2)
library(tidyr)
library(scales)
options(bitmapType = "cairo")

### Plotting overlap between ChromBPNet variants and caQTL variants
## This test checks correlation between:
##  - caQTL allelic effect size (beta; here stored in `slope`)
##  - ChromBPNet predicted allelic fold-change (here stored in `logfc.mean`)
## for variants located within peaks.

#####################
# Read ChromBPNet + caQTL merged results
ca_vareff <- read.delim2(
  "/gpfs/commons/home/mmatos/cd4_qtl_paper_figures/figure_2_chrombpnet/data/bpnet_lfc_caqtl_beta/CD4T_caQTL_nomical_assoc_inpeak_chrombpnet.csv",
  sep = ","
)

# Coerce numeric columns
ca_vareff$pval_nominal <- as.numeric(ca_vareff$pval_nominal)
ca_vareff$slope        <- as.numeric(ca_vareff$slope)
ca_vareff$logfc.mean   <- as.numeric(ca_vareff$logfc.mean)

# Define significance level order and custom colors
signif_levels <- c("Both", "caQTL", "chromBPnet", "Not Significant")
custom_colors <- c(
  "Both"            = "#e4b57a",
  "caQTL"           = "#9CCB86",
  "chromBPnet"      = "#0076C0",
  "Not Significant" = "gray90"
)

# Make sure `significant` is a factor with correct levels
ca_vareff <- ca_vareff %>%
  mutate(significant = factor(significant, levels = signif_levels))

head(ca_vareff)

# Count number of variants per group (for legend labels)
category_counts <- ca_vareff %>%
  count(significant) %>%
  complete(significant = signif_levels, fill = list(n = 0))

# Create legend labels with counts
color_labels <- category_counts %>%
  arrange(match(significant, signif_levels)) %>%
  mutate(label = paste0(significant, " (n=", n, ")")) %>%
  pull(label)

# ---------------------------
# Correlations (overall + by category)
# ---------------------------
cor_stats <- function(df, method = "pearson") {
  x <- df$slope
  y <- df$logfc.mean
  ok <- is.finite(x) & is.finite(y)
  n  <- sum(ok)
  
  if (n < 3) {
    return(data.frame(n = n, r = NA_real_, r2 = NA_real_, p = NA_real_))
  }
  
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method))
  r  <- unname(ct$estimate)
  
  data.frame(
    n  = n,
    r  = r,
    r2 = r^2,
    p  = ct$p.value
  )
}

# Overall correlation (all variants)
overall <- cor_stats(ca_vareff) %>%
  mutate(significant = "All variants") %>%
  dplyr::select(significant, everything())

# By significance class
by_group <- ca_vareff %>%
  group_by(significant) %>%
  group_modify(~ as_tibble(cor_stats(.x))) %>%
  ungroup()

cor_table <- bind_rows(overall, by_group) %>%
  mutate(
    significant = factor(significant, levels = c("All variants", signif_levels)),
    p_label = if_else(is.na(p), "NA", format.pval(p, digits = 2, eps = 1e-300)),
    line = sprintf("%s: R²=%.3f (r=%.3f), p=%s, n=%d", significant, r2, r, p_label, n)
  )

print(cor_table)

#save correlation table
write.csv(
  cor_table,
  "~/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/general_figures/cbpnet_vs_caqtl_beta_correlations.csv",
  row.names = FALSE
)

# Multi-line caption text (includes overall + each class)
stats_text <- paste(cor_table$line, collapse = "\n")

# ---------------------------
# Subset for "Both" (used for regression line only)
# ---------------------------
both_df <- dplyr::filter(ca_vareff, significant == "Both")

# Shuffle to reduce overplotting bias
set.seed(42)
ca_vareff <- ca_vareff[sample(nrow(ca_vareff)), ]

# ---------------------------
# Plot
# ---------------------------
ca_vareff_plot_r2 <- ggplot(ca_vareff, aes(x = slope, y = logfc.mean, fill = significant)) +
  ggrastr::geom_point_rast(shape = 21, size = 2, alpha = 0.6, stroke = 0.1, color = "gray10") +
  geom_smooth(
    data = both_df,
    aes(x = slope, y = logfc.mean),
    method = "lm",
    se = FALSE,
    color = custom_colors["Both"],
    size = 0.8
  ) +
  # diagonal reference
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  scale_fill_manual(values = custom_colors, labels = color_labels) +
  guides(fill = guide_legend(override.aes = list(size = 3), nrow = 2)) +
  labs(
    x = "caQTL (β)",
    y = "ChromBPNet log(aFC)",
    fill = "Significance:",
    title = "Variant Allelic Effect Correlation",
    caption = stats_text
  ) +
  theme_light() +
  theme(
    legend.title = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0)
  )

ca_vareff_plot_r2

ggsave(
  "~/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/general_figures/cbpnet_vs_caqtl_beta.pdf",
  ca_vareff_plot_r2,
  width = 5,
  height = 5
)


### looking at concoordance
dir_conc <- function(df) {
  x <- df$slope
  y <- df$logfc.mean
  ok <- is.finite(x) & is.finite(y) & x != 0 & y != 0   # drop zeros (no direction)
  n  <- sum(ok)
  if (n < 1) return(data.frame(n_dir = n, prop_concordant = NA_real_))
  data.frame(
    n_dir = n,
    prop_concordant = mean(sign(x[ok]) == sign(y[ok]))
  )
}

# Overall
dir_overall <- dir_conc(ca_vareff) %>%
  mutate(significant = "All variants")

# By group
dir_by_group <- ca_vareff %>%
  group_by(significant) %>%
  group_modify(~ as_tibble(dir_conc(.x))) %>%
  ungroup()

dir_table <- bind_rows(dir_overall, dir_by_group)
print(dir_table)


write.csv(
  dir_table,
  "~/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/general_figures/cbpnet_vs_caqtl_beta_concoorance.csv",
  row.names = FALSE
)

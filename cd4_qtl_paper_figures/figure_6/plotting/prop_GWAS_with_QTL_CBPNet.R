## ----------------------------------------------------------------------
## Setup
## ----------------------------------------------------------------------
.libPaths(c(
  "/gchm/R/x86_64-pc-linux-gnu-library/4.4",
  "/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"
))
options(bitmapType = "cairo")

library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(cowplot)

setwd("~/cd4_qtl_paper_figures/figure_6")
source("~/cd4_qtl_paper_figures/utils/color_pallete_helper.R")

## ----------------------------------------------------------------------
## User options
## ----------------------------------------------------------------------
TRAIT_SET_MODE <- "union"   # "union" or "intersect"

## ----------------------------------------------------------------------
## Helper: safe pretty labels
## ----------------------------------------------------------------------
label_trait <- function(trait_code) {
  dplyr::recode(trait_code, !!!trait_labels, .default = trait_code)
}

## ----------------------------------------------------------------------
## 1) QTL–GWAS colocalization summary (plot2 inputs)
## ----------------------------------------------------------------------

coloc_table <- fread(
  "~/cd4_qtl_paper_figures/figure_6/data/CD4T_coloc_summary_table_v2.csv"
)

coloc_table$trait <- gsub("_preprocessed$", "", coloc_table$trait)
coloc_table$trait_study <- coloc_table$trait

coloc_table_summary <- coloc_table %>%
  separate(
    trait_study,
    into   = c("year", "accession", "trait", "ethnicity"),
    sep    = "_",
    remove = FALSE
  ) %>%
  group_by(trait) %>%
  summarise(
    median_cond_indep = median(total_cond_indep, na.rm = TRUE),
    q1_cond_indep     = quantile(total_cond_indep, 0.25, na.rm = TRUE),
    q3_cond_indep     = quantile(total_cond_indep, 0.75, na.rm = TRUE),
    n_studies         = n(),
    total_cond_indep  = sum(total_cond_indep,  na.rm = TRUE),
    just_eqtl         = sum(just_eqtl,         na.rm = TRUE),
    just_caqtl        = sum(just_caqtl,        na.rm = TRUE),
    both              = sum(both,              na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    trait_label       = label_trait(trait),
    total_coloc       = just_eqtl + just_caqtl + both,
    frac_coloc        = ifelse(total_cond_indep > 0, total_coloc / total_cond_indep, 0),
    frac_caqtl_only_within_coloc = ifelse(total_coloc > 0, (just_caqtl + both) / total_coloc, 0)
  )

## ----------------------------------------------------------------------
## 2) ChromBPNet–GWAS membership summary (plot3 inputs)
## ----------------------------------------------------------------------

all_results <- fread(
  "~/cd4_qtl_paper_figures/figure_6/data/chrombpnet_immune_gwas_membership_all_results.csv"
)

trait_summary <- all_results %>%
  separate(
    gwas_id,
    into   = c("year", "accession", "trait", "ethnicity"),
    sep    = "_",
    remove = FALSE
  ) %>%
  group_by(trait) %>%
  summarise(
    total_cs_all.mean          = mean(total_cs_all,     na.rm = TRUE),
    total_cs_all.median        = median(total_cs_all,   na.rm = TRUE),
    q1_cond_indep              = quantile(total_cs_all, 0.25, na.rm = TRUE),
    q3_cond_indep              = quantile(total_cs_all, 0.75, na.rm = TRUE),
    total_cs_all               = sum(total_cs_all,      na.rm = TRUE),
    total_cs_eligible          = sum(total_cs_eligible, na.rm = TRUE),
    cs_with_cbpnet             = sum(cs_with_cbpnet,    na.rm = TRUE),
    cs_cbpnet_specific         = sum(cs_cbpnet_specific, na.rm = TRUE),
    cs_cbpnet_shared           = sum(cs_cbpnet_shared,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    trait_label                = label_trait(trait),
    prop_cs_with_cbpnet_pooled = ifelse(total_cs_all > 0, cs_with_cbpnet / total_cs_all, 0)
  )

## ----------------------------------------------------------------------
## 3) Harmonize trait sets across panels (union or intersect)
## ----------------------------------------------------------------------

if (TRAIT_SET_MODE == "intersect") {
  common_traits <- intersect(coloc_table_summary$trait, trait_summary$trait)
  
  coloc_table_summary <- coloc_table_summary %>% filter(trait %in% common_traits)
  trait_summary       <- trait_summary       %>% filter(trait %in% common_traits)
  
} else if (TRAIT_SET_MODE == "union") {
  all_traits <- union(coloc_table_summary$trait, trait_summary$trait)
  
  coloc_table_summary <- tibble(trait = all_traits) %>%
    left_join(coloc_table_summary, by = "trait") %>%
    mutate(
      trait_label = label_trait(trait),
      across(
        c(median_cond_indep, q1_cond_indep, q3_cond_indep, n_studies,
          total_cond_indep, just_eqtl, just_caqtl, both,
          total_coloc, frac_coloc, frac_caqtl_only_within_coloc),
        ~ replace_na(., 0)
      )
    )
  
  trait_summary <- tibble(trait = all_traits) %>%
    left_join(trait_summary, by = "trait") %>%
    mutate(
      trait_label = label_trait(trait),
      across(
        c(total_cs_all.mean, total_cs_all.median, q1_cond_indep, q3_cond_indep,
          total_cs_all, total_cs_eligible, cs_with_cbpnet,
          cs_cbpnet_specific, cs_cbpnet_shared, prop_cs_with_cbpnet_pooled),
        ~ replace_na(., 0)
      )
    )
} else {
  stop("TRAIT_SET_MODE must be 'union' or 'intersect'")
}

## ----------------------------------------------------------------------
## 4) FILTER TO ONLY TRAITS WITH QTL–GWAS COLOC > 0
##     (and apply same trait set to ChromBPNet panel)
## ----------------------------------------------------------------------

keep_traits <- coloc_table_summary %>%
  filter(total_coloc > 0) %>%
  pull(trait)

coloc_table_summary <- coloc_table_summary %>%
  filter(trait %in% keep_traits)

trait_summary <- trait_summary %>%
  filter(trait %in% keep_traits)

## ----------------------------------------------------------------------
## 5) Trait ordering (shared across panels AFTER FILTER)
## ----------------------------------------------------------------------

trait_order <- coloc_table_summary %>%
  arrange(total_cond_indep, trait_label) %>%
  pull(trait_label)

## ------------------- Coloc long data -------------------
plot_df <- coloc_table_summary %>%
  select(trait, trait_label, total_cond_indep,
         median_cond_indep, q1_cond_indep, q3_cond_indep,
         just_eqtl, just_caqtl, both) %>%
  pivot_longer(
    cols      = c(just_eqtl, just_caqtl, both),
    names_to  = "category",
    values_to = "count"
  ) %>%
  mutate(
    proportion  = ifelse(total_cond_indep > 0, count / total_cond_indep, 0),
    trait_label = factor(trait_label, levels = trait_order)
  )

avg_total_coloc <- coloc_table_summary %>%
  select(trait_label, total_cond_indep, median_cond_indep, q1_cond_indep, q3_cond_indep) %>%
  mutate(trait_label = factor(trait_label, levels = trait_order)) %>%
  arrange(total_cond_indep)

## ------------------- ChromBPNet long data -------------------
trait_long <- trait_summary %>%
  select(trait, trait_label, total_cs_all,
         cs_cbpnet_specific, cs_cbpnet_shared) %>%
  pivot_longer(
    cols      = c(cs_cbpnet_specific, cs_cbpnet_shared),
    names_to  = "cbpnet_cat",
    values_to = "cs_n"
  ) %>%
  mutate(
    cbpnet_cat = dplyr::recode(
      cbpnet_cat,
      cs_cbpnet_specific = "ChromBPNet_specific",
      cs_cbpnet_shared   = "QTL_shared"
    ),
    prop        = ifelse(total_cs_all > 0, cs_n / total_cs_all, 0),
    trait_label = factor(trait_label, levels = trait_order)
  )

avg_total_cbpnet <- trait_summary %>%
  distinct(trait_label, total_cs_all, q1_cond_indep, q3_cond_indep) %>%
  mutate(trait_label = factor(trait_label, levels = trait_order)) %>%
  arrange(total_cs_all)

trait_summary_lab <- trait_summary %>%
  filter(cs_with_cbpnet > 0) %>%
  mutate(trait_label = factor(trait_label, levels = trait_order))

## ----------------------------------------------------------------------
## 6) PANEL-SPECIFIC SCALING FACTORS
## ----------------------------------------------------------------------

max_total_coloc <- max(avg_total_coloc$total_cond_indep, na.rm = TRUE)
scaling_factor_coloc <- ifelse(max_total_coloc > 0, 1 / max_total_coloc, 1)

max_total_cbpnet <- max(avg_total_cbpnet$total_cs_all, na.rm = TRUE)
scaling_factor_cbpnet <- ifelse(max_total_cbpnet > 0, 1 / max_total_cbpnet, 1)

## ----------------------------------------------------------------------
## 7) Plot 2 – QTL–GWAS colocalization
## ----------------------------------------------------------------------

plot2 <- ggplot(plot_df, aes(x = trait_label, y = proportion, fill = category)) +
  geom_col(color = "gray40", linewidth = 0.1) +
  geom_text(
    data = subset(plot_df, count != 0),
    aes(label = count),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 6
  ) +
  geom_line(
    data = avg_total_coloc,
    aes(x = trait_label, y = total_cond_indep * scaling_factor_coloc, group = 1),
    inherit.aes = FALSE,
    color = "gray50",
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  geom_point(
    data = avg_total_coloc,
    aes(x = trait_label, y = total_cond_indep * scaling_factor_coloc),
    inherit.aes = FALSE,
    shape = 21,
    fill = "gray30",
    size = 4,
    stroke = 0.3,
    alpha = 0.5
  ) +
  scale_fill_manual(values = qtl_cmap_coloc) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0.01, 0),
    name   = "Proportion colocalized GWAS loci (bars)",
    sec.axis = sec_axis(
      ~ . / scaling_factor_coloc,
      breaks = pretty(c(0, max_total_coloc)),
      name   = "Total number loci per trait study (line)"
    )
  ) +
  coord_flip() +
  labs(
    title = "Summary of QTL–GWAS colocalization",
    x     = "GWAS trait"
  ) +
  theme_cowplot() +
  theme(
    axis.text.x        = element_text(size = 14),
    axis.text.y        = element_text(size = 18),
    axis.title         = element_text(size = 14, face = "bold"),
    legend.position    = "bottom",
    legend.justification = "center",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 14),
    axis.text.x.bottom = element_text(size = 12, color = "black"),
    axis.text.x.top    = element_text(size = 12, color = "gray50"),
    axis.title.x.top   = element_text(size = 14, color = "gray50")
  )


coloc_med_by_cat <- plot_df %>%
  group_by(category) %>%
  summarise(
    n_traits     = sum(!is.na(proportion)),
    median_prop  = median(proportion, na.rm = TRUE),
    q1_prop      = quantile(proportion, 0.25, na.rm = TRUE),
    q3_prop      = quantile(proportion, 0.75, na.rm = TRUE),
    iqr_prop     = q3_prop - q1_prop,
    .groups = "drop"
  )

coloc_med_by_cat

coloc_all <- plot_df %>%
  group_by(trait, trait_label) %>%
  summarise(prop_all_coloc = sum(proportion, na.rm = TRUE), .groups = "drop") %>%
  summarise(
    category    = "all_coloc",
    n_traits    = sum(!is.na(prop_all_coloc)),
    median_prop = median(prop_all_coloc, na.rm = TRUE),
    q1_prop     = quantile(prop_all_coloc, 0.25, na.rm = TRUE),
    q3_prop     = quantile(prop_all_coloc, 0.75, na.rm = TRUE),
    iqr_prop    = q3_prop - q1_prop
  )

coloc_all

## ----------------------------------------------------------------------
## 8) Plot 3 – ChromBPNet membership (same filtered traits)
## ----------------------------------------------------------------------

plot3 <- ggplot(trait_long, aes(x = trait_label, y = prop, fill = cbpnet_cat)) +
  geom_col(color = "gray40", linewidth = 0.1) +
  geom_text(
    data = trait_summary_lab,
    aes(x = trait_label, y = prop_cs_with_cbpnet_pooled / 2, label = cs_with_cbpnet),
    inherit.aes = FALSE,
    color = "black",
    size  = 6
  ) +
  scale_fill_manual(
    values = c(
      "ChromBPNet_specific" = "#0071ba",
      "QTL_shared"          = "#ff9f1c"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 0.30),
    breaks = seq(0, 0.30, 0.05),
    expand = c(0.01, 0),
    name   = "Proportion GWAS CS with CBPNet variants (bars)"
  ) +
  coord_flip() +
  theme_cowplot() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x        = element_text(size = 14),
    axis.title.y       = element_blank(),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.title         = element_text(size = 16, face = "bold"),
    legend.position    = "bottom",
    legend.justification = "center",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 14)
  )

med_by_cat <- trait_long %>%
  group_by(cbpnet_cat) %>%
  summarise(
    n_traits = sum(!is.na(prop)),
    median_prop = median(prop, na.rm = TRUE),
    q1_prop = quantile(prop, 0.25, na.rm = TRUE),
    q3_prop = quantile(prop, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

med_by_cat
## ----------------------------------------------------------------------
## 9) Combine panels
## ----------------------------------------------------------------------

plt <- plot_grid(
  plot2, plot3,
  rel_widths = c(1, 0.4),
  align      = "h"
)

plt

ggsave("~/cd4_qtl_paper_figures/figure_6/plotting/plots/cbpnet_gwas_coloc_proportions_by_category.pdf",plt, height = 17, width = 17, dpi=300, units = "in")


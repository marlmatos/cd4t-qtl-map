############
## plotting figure 1- main "caQTL Categories"
## Author: Marliette Matos
#############

library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(GenomicRanges)
library(rcartocolor)
library(tidyr)
options(bitmapType = "cairo")

source("~/cd4_qtl_paper_figures/utils/color_pallete_helper.R")

summary_df<-fread("/gchm/cd4_qtl_paper_figures/figure_1/data/eqtl_caqtl_finemapping_coloc_all.tsv")

#simplifying caQTL categories
summary_df <- summary_df %>%
  mutate(
    caqtl_category = ifelse(
      cs_type_any_var %in% c("in_corr_Peak", "in_other_Peak"),
      "in_other_Peak",
      cs_type_any_var
    )
  )

##plot caqtl categories breakdown
label_map <- c(
  in_caPeak = "Within caQTL peak",
  in_other_Peak = "In other open peak",
  no_Peak_overlap = "No peak overlap"
)

# Plot with custom fill scale

plot_df<-summary_df %>% 
  filter(coloc_status %in% c("coloc", "caQTL_only")) %>% 
  group_by(caqtl_category) %>%
  summarise(n = n_distinct(finemapped_cs)) %>%
  mutate(unit = "caQTLs",
         pct = n / sum(n) * 100) 

fwrite(plot_df, "~/cd4_qtl_paper_figures/figure_1/data/caqtl_categories_pct.csv")

stack_plot<-ggplot(plot_df, aes(x = unit, y = pct, fill = caqtl_category)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 3
  ) +
  scale_fill_manual(
    values = caQTL_category_cmap2,
    labels = label_map,
    name = "Category"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),   # keep axis as percent
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "caQTL Categories",
    x = NULL,
    y = "Percent"
  ) +
  theme_minimal()
stack_plot
ggsave("~/cd4_qtl_paper_figures/figure_1/plotting/plots_oct2025/peaks_vs_cs_plot_stacked_pct.pdf", stack_plot, width = 4, height = 5)


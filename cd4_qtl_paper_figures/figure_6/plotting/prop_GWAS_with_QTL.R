
.libPaths(c("/gchm/R/x86_64-pc-linux-gnu-library/4.4","/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"))
options(bitmapType = "cairo")
library(ggplot2)
library(data.table)
library(tidyverse)
library(forcats)
library(dplyr)
library(cowplot)
setwd("~/cd4_qtl_paper_figures/fiure_6")
source("~/cd4_qtl_paper_figures/utils/color_pallete_helper.R")

data.dir="/gcgl/sghatan/marlis_pj/coloc/"
# folder names
file_names = unlist(fread(paste0(data.dir, "preprocessed_folder_names.txt"), header = F))

count_cs <- function(file_name){
  region <- readRDS(file_name)
  if (isFALSE(region$converged)) 1L else length(region$sets$cs_index)
}

make_table_union <- function(gwas_id) {
  message("Processing: ", gwas_id)
  gwas_dir <- "/gcgnl/finemapping_autoimmune/"
  
  # ----- finemapping (can be missing; just set to NA/0)
  regions <- list.files(
    path = file.path(gwas_dir, "output/nathan_completed", gwas_id),
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  regions <- regions[grepl("_1.5Mb/susie/", regions)]
  total_cs <- if (length(regions) == 0) NA_integer_ else sum(sapply(regions, count_cs))
  
  # ----- coloc paths (either side may be missing)
  eqtl_path  <- file.path(data.dir, "coloc_results/eqtl_gwas_coloc", paste0(gwas_id, "_coloc_results.csv"))
  caqtl_path <- file.path(data.dir, "coloc_results/ca_gwas_coloc", paste0("CD4T_chromatin_", gwas_id, "_coloc_results.csv"))
  
  eqtl_exists  <- file.exists(eqtl_path)
  caqtl_exists <- file.exists(caqtl_path)
  
  # helper: read+dedup if file exists, else empty dt with region_cs col
  read_eqtl <- function(p) {
    if (!file.exists(p)) return(data.table(region_cs = character()))
    dt <- fread(p, showProgress = FALSE)
    if (nrow(dt) == 0) return(data.table(region_cs = character()))
    dt <- as.data.table(dt)
    dt[, region_cs := paste(region, idx1, sep = "_")]
    dt <- dt[pval < 1e-5 & pval_nominal < 1e-3]
    if (nrow(dt) == 0) return(data.table(region_cs = character()))
    # keep best per (variant,gene)
    dt <- dt[order(-PP.H4.abf)]
    dt <- dt[, .SD[1], by = .(eQTL_variant_GRC37, gene)]
    dt
  }
  
  read_caqtl <- function(p) {
    if (!file.exists(p)) return(data.table(region_cs = character()))
    dt <- fread(p, showProgress = FALSE)
    if (nrow(dt) == 0) return(data.table(region_cs = character()))
    dt <- as.data.table(dt)
    dt[, region_cs := paste(region, idx1, sep = "_")]
    dt <- dt[pval < 1e-5 & pval_nominal < 1e-3]
    if (nrow(dt) == 0) return(data.table(region_cs = character()))
    # keep best per (variant,peak)
    dt <- dt[order(-PP.H4.abf)]
    dt <- dt[, .SD[1], by = .(eQTL_variant_GRC37, peak)]
    dt
  }
  
  eqtl_coloc_result  <- read_eqtl(eqtl_path)
  caqtl_coloc_result <- read_caqtl(caqtl_path)
  
  # shared on region_cs
  shared_coloc_result <- caqtl_coloc_result[region_cs %in% eqtl_coloc_result$region_cs]
  shared_coloc_result[, gwas_id := gwas_id]
  
  # Summary counts (safe even if one side empty)
  only_eqtl  <- length(setdiff(unique(eqtl_coloc_result$region_cs),  unique(caqtl_coloc_result$region_cs)))
  only_caqtl <- length(setdiff(unique(caqtl_coloc_result$region_cs), unique(eqtl_coloc_result$region_cs)))
  both       <- length(unique(shared_coloc_result$region_cs))
  
  summary_table <- data.frame(
    trait = gwas_id,
    total_cond_indep = total_cs,
    just_eqtl  = only_eqtl,
    just_caqtl = only_caqtl,
    both = both,
    has_eqtl_file  = eqtl_exists,
    has_caqtl_file = caqtl_exists
  )
  
  list(summary = summary_table, shared = shared_coloc_result, eqtl = eqtl_coloc_result)
}


combined_table <- lapply(file_names, make_table_union)

results2 <- bind_rows(lapply(combined_table, `[[`, "summary"))
fwrite(results2, "~/cd4_qtl_paper_figures/fiure_6/data/CD4T_coloc_summary_table_v2.csv")

all_shared <- bind_rows(lapply(combined_table, `[[`, "shared"))
fwrite(all_shared, "~/cd4_qtl_paper_figures/fiure_6/data/CD4T_all_shared_coloc_v2.csv")

all_eqtls <- bind_rows(lapply(combined_table, function(x) {
  as.data.frame(x$eqtl) %>% mutate(trait = x$summary$trait[1])
}))
fwrite(all_eqtls, "~/cd4_qtl_paper_figures/fiure_6/data/CD4T_all_eqtl_coloc_v2.csv")

results2 = fread("~/cd4_qtl_paper_figures/fiure_6/data/CD4T_coloc_summary_table_v2.csv")

# Clean `trait` column (remove "_preprocessed")
results2$trait <- gsub("_preprocessed$", "", results2$trait)

#assign to `trait_study` before separation if needed
results2$trait_study <- results2$trait

# Create summary table
results2_summary <- results2 %>%
  separate(trait_study, into = c("year", "accession", "trait", "ethnicity"), sep = "_", remove = FALSE) %>%
  group_by(trait) %>%
  summarise(
    mean_cond_indep = mean(total_cond_indep, na.rm = TRUE),
    total_cond_indep = sum(total_cond_indep, na.rm = TRUE),
    just_eqtl        = sum(just_eqtl, na.rm = TRUE),
    just_caqtl       = sum(just_caqtl, na.rm = TRUE),
    both             = sum(both, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    trait_label = recode(trait, !!!trait_labels, .default = trait)  # <- key fix
  )

# take first 10
#results2 = results2[1:20,]

plot_df = pivot_longer(
  results2_summary,
  cols = c(just_eqtl, just_caqtl, both),
  names_to = "category",
  values_to = "count"
) 
plot_df$proportion = plot_df$count / plot_df$total_cond_indep


##plot
# Compute sorting order based on "both" proportion
trait_order <- plot_df %>%
  filter(category == "both") %>%
  arrange((count)) %>%
  pull(trait_label)

#Set trait as a factor with that order
plot_df <- plot_df %>%
  mutate(trait_label = factor(trait_label, levels = trait_order))

# Plot
plot <- plot_df %>%
  ggplot(aes(x = trait_label, y = proportion, fill = category)) +
  geom_bar(stat = "identity", color = "gray30", size = 0.1) +
  geom_text(
    data = subset(plot_df, count != 0),
    aes(label = count),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 6
  ) +
  scale_fill_manual(values = qtl_cmap_coloc) +
  labs(
    title = "Proportion of Conditional Independent Variants",
    y = "Proportion Variants Colocalized",
    x = "GWAS Trait"
  ) +
  coord_flip() + 
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "top",
    legend.justification = "center",
    legend.title = element_blank(),
    legend.text = element_text(size = 16)
  )

plot

ggsave("~/cd4_qtl_paper_figures/figure_2/plotting/plots/gwas_coloc_proportions.pdf", height = 16, width = 7, dpi=300, units = "in")


library(dplyr)
library(ggplot2)

# 1) Keep only traits with any nonzero counts
keep_traits <- plot_df %>%
  group_by(trait_label) %>%
  summarise(has_nonzero = any(count > 0), .groups = "drop") %>%
  filter(has_nonzero) %>%
  pull(trait_label)

plot_df <- plot_df %>%
  filter(trait_label %in% keep_traits)

results2_summary_nz <- results2_summary %>%
  filter(trait_label %in% keep_traits)

# 2) Define order (trait-level!)
trait_order <- results2_summary_nz %>%
  arrange(mean_cond_indep) %>%
  pull(trait_label)

# 3) Apply factor levels consistently
plot_df <- plot_df %>%
  mutate(trait_label = factor(trait_label, levels = trait_order))

avg_total_df <- results2_summary_nz %>%
  select(trait_label, mean_cond_indep) %>%
  mutate(trait_label = factor(trait_label, levels = trait_order))

# 4) Scaling factor (safe)
max_prop  <- max(plot_df$proportion, na.rm = TRUE)
max_total <- max(avg_total_df$mean_cond_indep, na.rm = TRUE)
scaling_factor <- ifelse(max_total > 0, max_prop / max_total, 1)

# 5) Plot
plot2 <- ggplot(plot_df, aes(x = trait_label, y = proportion, fill = category)) +
  geom_col(color = "gray40", linewidth = 0.1) +
  geom_text(
    data = subset(plot_df, count != 0),
    aes(label = count),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 6
  ) +
  geom_point(
    data = avg_total_df,
    aes(x = trait_label, y = mean_cond_indep * scaling_factor),
    inherit.aes = FALSE,
    shape = 21, fill = "gray30", size = 3, stroke = 0.3, alpha = 0.5
  ) +
  geom_line(
    data = avg_total_df,
    aes(x = trait_label, y = mean_cond_indep * scaling_factor, group = 1),
    inherit.aes = FALSE,
    color = "gray50", linewidth = 0.6, linetype = "dashed"
  ) +
  scale_fill_manual(values = qtl_cmap_coloc) +
  scale_y_continuous(
    expand = c(0.01, 0),
    name = "Proportion Variants Colocalized (Bars)",
    sec.axis = sec_axis(~ . / scaling_factor,
                        name = "Average Number of Variants per Trait (Dots)")
  ) +
  coord_flip() +
  labs(title = "Summary of QTL-GWAS Colocalization", x = "GWAS Trait") +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title  = element_text(size = 16, face = "bold"),
    legend.position = "right",
    legend.justification = "center",
    legend.title = element_blank(),
    legend.text  = element_text(size = 16),
    axis.text.x.bottom = element_text(size = 14, color = "black"),
    axis.text.x.top    = element_text(size = 12, color = "gray50"),
    axis.title.x.top   = element_text(size = 14, color = "gray50")
  )

plot2
ggsave("~/cd4_qtl_paper_figures/fiure_6/plotting/plots/gwas_coloc_proportions.pdf",plot2, height = 17, width = 14, dpi=300, units = "in")

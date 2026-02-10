############
## plotting figure 1- main
## Author: Marliette Matos
#############

library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyverse)
options(bitmapType = "cairo")

#######Summarize eQTL_results
eqtl_res <- data.table()
per_chr  <- vector("list", 22)

for (i in 1:22) {
  ## ---- CI eQTLs (row-level, remove only identical rows) ----
  f_ci <- sprintf("/gchm/cd4_QTL_analysis/03_Run_cisQTL_perchr/analysis/cis_eQTLs/all_CD4T_cells_MAF5/allcells.independent_cis_qtl_pairs.chr%d.csv", i)
  tmp  <- fread(f_ci) %>%
    filter(pval_perm < 0.05) %>%
    distinct() %>%  # identical rows only
    mutate(
      chr  = paste0("chr", i),
      eqtl = paste(phenotype_id, "-", variant_id, sep = "")
    )
  
  ## ---- Finemap credible sets: DO NOT collapse before per-CS check ----
  f_fm <- sprintf("/gcgls/marlis_pj/coloc/SuSiE_finemap_credible_sets/All_CD4T_cells/All_CD4T_cells_chr%d_credible_sets.txt", i)
 
  fm_raw <- fread(f_fm)
  
  # normalize variant id if needed (e.g., 'snp' -> 'variant_id')
  if (!"variant_id" %in% names(fm_raw) && "snp" %in% names(fm_raw)) {
    setnames(fm_raw, "snp", "variant_id")
  }
  
  fm_raw <- fm_raw %>%
    mutate(
      unique_id = paste(region, gene, cs),
      eqtl      = paste(gene, "-", variant_id, sep = "")
    )
  
  ## ---- Per-credible-set overlap with CI eQTL pairs ----
  cs_status <- fm_raw %>%
    group_by(unique_id) %>%
    summarise(
      n_vars_in_cs   = n(),
      n_overlap_with_CI = sum(eqtl %in% tmp$eqtl),
      any_in_CI      = n_overlap_with_CI > 0,
      .groups = "drop"
    )
  
  # CS that appeared only after finemapping (no CI variant in that CS)
  cs_only_in_fm <- cs_status %>% filter(!any_in_CI)
  
  ## ---- Flag CI rows that are in ANY credible set (pair-level) ----
  # Build pair-level set AFTER per-CS check to avoid losing CS info
  fm_pairs <- fm_raw %>% distinct(eqtl)
  tmp$finemapping_status <- tmp$eqtl %in% fm_pairs$eqtl
  
  ## ---- Collect outputs ----
  eqtl_res <- rbind(eqtl_res, tmp, fill = TRUE)
  
  per_chr[[i]] <- data.frame(
    chr                = paste0("chr", i),
    # CI row-level
    n_ci_rows          = nrow(tmp),
    n_fm_rows          = sum(tmp$finemapping_status),
    prop_fm_rows       = ifelse(nrow(tmp) == 0, NA_real_, mean(tmp$finemapping_status)),
    # CS-level
    n_cs_total         = n_distinct(cs_status$unique_id),
    n_cs_only_in_fm    = nrow(cs_only_in_fm)
  )
  
}

per_chr_df <- bind_rows(per_chr) %>%
  arrange(chr)

print(per_chr_df)

# If you want overall summaries (keep units separate!)
cat(
  "\nRow-level (conditionally independent (CI) eQTL table):\n",
  "  total CI rows = ", nrow(eqtl_res), "\n",
  "  finemapped CI rows = ", sum(eqtl_res$finemapping_status), "\n",
  "  prop finemapped CI rows = ", round(mean(eqtl_res$finemapping_status), 4), "\n",
  sep = ""
)

###  "Summary independent eQTLs that appeared after finemapping ######
cat( 
  "\nCS-level:",
  "  total CS = ", sum(per_chr_df$n_cs_total), "\n",
  "  CS only in finemapping (no CI overlap) = ", sum(per_chr_df$n_cs_only_in_fm), "\n",
  sep = ""
)


# Row-level (after your all-columns distinct)
total_rows      <- nrow(eqtl_res)
n_eGenes        <- uniqueN(eqtl_res$phenotype_id)
n_rows_fm       <- eqtl_res[, sum(finemapping_status, na.rm = TRUE)]
prop_rows_fm    <- eqtl_res[, mean(finemapping_status, na.rm = TRUE)]
# eGenes with ≥1 fine-mapped variant
n_egenes_total <- uniqueN(eqtl_res$phenotype_id)
egenes_fm <- unique(eqtl_res[finemapping_status == TRUE, phenotype_id])
n_egenes_fm <- length(egenes_fm)
prop_egenes_fm <- if (n_egenes_total) n_egenes_fm / n_egenes_total else NA_real_

cat(
  "**eQTL finemapping summary**\n",
  "TOTAL eQTL rows (after all-columns distinct): ", total_rows, "\n",
  "TOTAL eGenes: ", n_eGenes, "\n",
  "TOTAL finemapped eQTLs: ", n_rows_fm, "\n",
  "Proportion finemapped eQTLs: ", round(prop_rows_fm, 4), "\n",
  "Number of distinct Fine-mapped eGenes: ", n_egenes_fm, "\n",
  "Proportion eGenes fine-mapped: ", round(prop_egenes_fm, 4), "\n",
  sep = ""
)


eqtl_summary<-per_chr_df  %>% summarise(total_ci=sum(n_ci_rows),
                                        total_in_fm=sum(n_cs_only_in_fm),
                                        totalfm=sum(n_cs_only_in_fm+n_fm_rows)) %>% mutate(total=(total_ci+total_in_fm),
                                                                                               n_phenotypes=n_eGenes) %>%
  pivot_longer(cols = everything(), names_to = "unit", values_to = "count") 

eqtl_summary$type=ifelse(eqtl_summary$unit=="n_phenotypes", "eGenes", "cis-eQTLs")
eqtl_summary$type2<-"cis-eQTLs"


########## caQTL Finemapping ##############
caqtl_res  <- data.table()
per_chr_ca <- vector("list", 22)


for (i in 1:22) {
  ## ---- CI caQTLs (row-level; remove only identical rows) ----
  f_ci <- sprintf("/gchm/cd4_caQTL_analysis/variant_to_peak_QTL/run_012625_qc_aware_qsmooth_CPM_MAF5_FDR5_1MB/results/006_caQTLs/filtered_qsmooth_norm_cpm_1mb/cd4_qsmooth_cpm_chromatin_narrowpeaks.independent_cis_QTL_pairs.chr%d.csv", i)
  tmp  <- fread(f_ci) %>%
    filter(pval_perm < 0.05) %>%
    distinct() %>%   # identical rows only
    mutate(
      chr   = paste0("chr", i),
      caqtl = paste(phenotype_id, "-", variant_id, sep = "")  # phenotype_id should be your peak id here
    )
  
  ## ---- Finemap credible sets (do NOT collapse before CS check) ----
  f_fm  <- sprintf("/gcgls/marlis_pj/coloc/SuSiE_finemap_credible_sets/CD4T_chromatin/CD4T_chromatin_chr%d_credible_sets.txt", i)
  fm_raw <- fread(f_fm)
  
  # Normalize variant/peak column names if needed
  if (!"variant_id" %in% names(fm_raw) && "snp" %in% names(fm_raw)) {
    setnames(fm_raw, "snp", "variant_id")
  }
  # If your finemap file uses a different peak column name, add it here:
  # e.g., if it's "peak_id" instead of "peak":
  if (!"peak" %in% names(fm_raw) && "peak_id" %in% names(fm_raw)) {
    setnames(fm_raw, "peak_id", "peak")
  }
  
  fm_raw <- fm_raw %>%
    mutate(
      unique_id = paste(region, peak, cs),             # CS identity
      caqtl     = paste(peak, "-", variant_id, sep = "")
    )
  
  ## ---- Per-credible-set overlap with CI caQTL pairs ----
  cs_status_ca <- fm_raw %>%
    group_by(unique_id) %>%
    summarise(
      n_vars_in_cs     = n(),
      n_overlap_with_CI= sum(caqtl %in% tmp$caqtl),
      any_in_CI        = n_overlap_with_CI > 0,
      .groups = "drop"
    )
  
  # CS that appear only after finemapping (no overlap with CI caQTL pairs)
  cs_only_in_fm_ca <- cs_status_ca %>% filter(!any_in_CI)
  
  ## ---- Pair-level set for row flagging (AFTER CS check to avoid losing CS info) ----
  fm_pairs_ca <- fm_raw %>% distinct(caqtl)
  tmp$finemapping_status <- tmp$caqtl %in% fm_pairs_ca$caqtl
  
  ## ---- Collect ----
  caqtl_res <- rbind(caqtl_res, tmp, fill = TRUE)
  
  per_chr_ca[[i]] <- data.frame(
    chr              = paste0("chr", i),
    # Row-level (CI table)
    n_ci_rows        = nrow(tmp),
    n_fm_rows        = sum(tmp$finemapping_status),
    prop_fm_rows     = ifelse(nrow(tmp) == 0, NA_real_, mean(tmp$finemapping_status)),
    # CS-level
    n_cs_total       = n_distinct(cs_status_ca$unique_id),
    n_cs_only_in_fm  = nrow(cs_only_in_fm_ca)
  )
  
}

per_chr_df_ca <- bind_rows(per_chr_ca) %>% arrange(chr)
print(per_chr_df_ca)

## ===== Row-level (after all-columns distinct) =====
total_rows_caqtl    <- nrow(caqtl_res)
n_capeaks_total     <- uniqueN(caqtl_res$phenotype_id)
n_rows_fm_caqtl     <- caqtl_res[, sum(finemapping_status, na.rm = TRUE)]
prop_rows_fm_caqtl  <- caqtl_res[, mean(finemapping_status, na.rm = TRUE)]

## caPeaks with ≥1 fine-mapped variant (distinct peak ids among finemapped rows)
capeak_fm           <- unique(caqtl_res[finemapping_status == TRUE, phenotype_id])
n_capeak_fm         <- length(capeak_fm)
prop_capeak_fm      <- if (n_capeaks_total) n_capeak_fm / n_capeaks_total else NA_real_

cat(
  "**caQTL finemapping summary**\n",
  "TOTAL caQTL rows (after all-columns distinct): ", total_rows_caqtl, "\n",
  "TOTAL caPeaks: ", n_capeaks_total, "\n",
  "TOTAL finemapped caQTL rows: ", n_rows_fm_caqtl, "\n",
  "Proportion finemapped caQTL rows: ", round(prop_rows_fm_caqtl, 4), "\n",
  "Distinct fine-mapped caPeaks: ", n_capeak_fm, "\n",
  "Proportion caPeaks fine-mapped: ", round(prop_capeak_fm, 4), "\n",
  sep = ""
)

## ===== Summary independent caQTLs that appeared after finemapping =====
cat(
  "\n**caQTL credible-set summary**\n",
  "Total CS (sum over chr): ", sum(per_chr_df_ca$n_cs_total), "\n",
  "CS only in finemapping (no CI overlap): ", sum(per_chr_df_ca$n_cs_only_in_fm), "\n",
  sep = ""
)


caqtl_summary<-per_chr_df_ca  %>% summarise(total_ci=sum(n_ci_rows),
                                            total_in_fm=sum(n_cs_only_in_fm),
                                            totalfm=sum(n_cs_only_in_fm+n_fm_rows)) %>% mutate(total=(total_ci+total_in_fm),
                                                                                               n_phenotypes=n_capeaks_total) %>%
pivot_longer(cols = everything(), names_to = "unit", values_to = "count") 
caqtl_summary$type=ifelse(caqtl_summary$unit=="n_phenotypes", "caPeaks", "cis-caQTLs")
caqtl_summary$type2<-"cis-caQTLs"


summary_all<-rbind(eqtl_summary, caqtl_summary) %>%
  filter(unit %in% c("total", "totalfm", "n_phenotypes")) 

# 1) Lookup with total + totalfm per type2
tot_lookup <- summary_all %>%
  filter(unit %in% c("total", "totalfm")) %>%
  select(type2, unit, count) %>%
  pivot_wider(names_from = unit, values_from = count)  # -> cols: total, totalfm

# 2) Keep only 'total' and 'n_phenotypes' rows, join totals, and set NA for phenotype rows
summary_with_fm <- summary_all %>%
  filter(unit %in% c("total", "n_phenotypes")) %>%
  left_join(tot_lookup, by = "type2") %>%
  mutate(
    total   = if_else(unit == "total", total,   NA_integer_),
    totalfm = if_else(unit == "total", totalfm, NA_integer_)
  ) %>%
  # (optional) nice ordering
  arrange(type2, factor(unit, levels = c("total", "n_phenotypes")))

summary_with_fm
#######

library(dplyr)
library(ggplot2)
library(scales)

df <- summary_with_fm %>% 
  filter(type2 == "cis-eQTLs") %>%
  mutate(
    unit_disp = dplyr::recode(unit, "total" = "eQTLs", "n_phenotypes" = "eGenes")
  )

# fraction label for the eQTL bar
frac_label <- df %>% 
  filter(unit == "total") %>% 
  transmute(
    x_pos = total / 2,
    lab   = scales::percent(totalfm / total, accuracy = 0.1)
  )

p1<-ggplot(df, aes(y = unit_disp, x = count)) +
  # base totals (gray = includes both fm + not fm)
  geom_col(width = 0.8, fill = c("#E0E0E0","#EFC7E6" ), color = "gray20") +
  # overlay just the fine-mapped portion for eQTLs
  geom_col(
    data = df %>% filter(unit == "total"),
    aes(x = totalfm),
    width = 0.8, fill = "#C70E7B", color = "gray20"
  ) +
  # total at the end of each bar
  geom_text(aes(label = comma(count)), hjust = -0.15, color = "gray20", size = 5) +
  # fraction label centered inside the eQTL bar
  geom_text(
    data = frac_label,
    aes(y = "eQTLs", x = x_pos, label = lab),
    color = "white", fontface = "bold", size = 5
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.15))) +
  labs(
    title = "cis-eQTL mapping (fine-mapped coverage)",
    subtitle = "Gray = total; magenta = fine-mapped portion; label = fraction fine-mapped",
    x = "Count", y = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 12)
  )

#--------caqtl plot--------
df <- summary_with_fm %>% 
  filter(type2 == "cis-caQTLs") %>%
  mutate(
    unit_disp = dplyr::recode(unit, "total" = "caQTLs", "n_phenotypes" = "caPeaks")
  )

# fraction label for the eQTL bar
frac_label <- df %>% 
  filter(unit == "total") %>% 
  transmute(
    x_pos = total / 2,
    lab   = scales::percent(totalfm / total, accuracy = 0.1)
  )

p2<-ggplot(df, aes(y = unit_disp, x = count)) +
  # base totals (gray = includes both fm + not fm)
  geom_col(width = 0.8, fill = c("#E0E0E0","#B4D9CC"), color = "gray20") +
  # overlay just the fine-mapped portion for caQTLs
  geom_col(
    data = df %>% filter(unit == "total"),
    aes(x = totalfm),
    width = 0.8, fill = "#9ccb86", color = "gray20"
  ) +
  # total at the end of each bar
  geom_text(aes(label = comma(count)), hjust = -0.15, color = "gray20", size = 5) +
  # fraction label centered inside the caQTLs bar
  geom_text(
    data = frac_label,
    aes(y = "caQTLs", x = x_pos, label = lab),
    color = "white", fontface = "bold", size = 5
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.15))) +
  labs(
    title = "cis-caQTLs mapping (fine-mapped coverage)",
    subtitle = "Gray = total; magenta = fine-mapped portion; label = fraction fine-mapped",
    x = "Count", y = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 12)
  )



plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1, 1))

qtl.plt<-plot_grid(
  ggdraw() + draw_plot(p1, x = 0, y = 0, width = 0.80, height = 1),
  ggdraw() + draw_plot(p2, x = 0, y = 0, width = 1.00, height = 1),
  ncol = 1
)

ggsave("~/cd4_qtl_paper_figures/figure_1/plotting/plots_oct2025/qtl_mapping_barplot.pdf", qtl.plt, width = 4, height = 5)

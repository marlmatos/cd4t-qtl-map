.libPaths(c("/gpfs/commons/home/mmatos/R/x86_64-pc-linux-gnu-library/4.4","/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"))

library(ggplot2)
library(dplyr)
library(readr) 
library(stringr)
library(GenomicRanges)
library(purrr)
library(LaCroixColoR)
library(rcartocolor)
library(ComplexUpset)
library(tidyr)
library(ChIPseeker)
library(stringr)
library(scales) 
library(data.table)
library(forcats)
library(cowplot)
library(patchwork)

options(bitmapType = "cairo")
cat("Reading data") 
out.dir="/gpfs/commons/home/mmatos/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/"
#read the finemapped coloc res
summary_df<-read_delim("/gpfs/commons/home/mmatos/cd4_qtl_paper_figures/figure_1/data/eqtl_caqtl_finemapping_coloc_all.tsv")
#summary_df<-read_delim("/gpfs/commons/groups/lappalainen_lab/mmatos/cd4_aging_project/cd4_qtl_data/colocalization/caqtl_eqtl_coloc_finemapped/eqtl_caqtl_finemap_coloc_summary_chr19.tsv")
summary_df<-summary_df %>% mutate(caqtl_category=ifelse(cs_type_any_var %in% c("in_other_Peak","in_corr_Peak"), "in_other_Peak", cs_type_any_var)) 

#read chromBPNet results
cBPNet_var<-read_delim("~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap/cd4_top_cpbnet_variants_motif_overlap_moreTFS_boolean.tsv") %>% dplyr::select(chr, pos, allele1, allele2, variant_id)
motifs<-read_delim("~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap_hit_caller/cd4_top_cpbnet_variants_motif_overlap_moreTFS_qtl_overlap_boolean_v2.tsv")

cBPNet_var <- cBPNet_var %>%
  left_join(
    motifs %>% dplyr::select(variant_id, category),
    by = "variant_id"
  ) %>%
  mutate(
    motif_overlap = !is.na(category)   # TRUE if in motifs, FALSE otherwise
  )

head(cBPNet_var)

#Please note that the summary_df already contains the joined colocalization results of the finemapped eQTL and caQTL variants, where 
# the column finemapped_cs_coloc already takes into account credible sets from both sets that are statistically the same based on H4 p>0.5.
# in the next commands I first choose the variant from each sets that might match the chrombpnet variant, but select instead the common finedmapped credible set name
# then I deduplicate them.

# get caQTLs credible sets that contain chromBPNet variants
chrombpnet_caqtlvars <- summary_df %>%
  mutate(var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])",
                       "chr\\1_\\2_\\3_\\4",
                       variant_id.x)) %>%
  semi_join(cBPNet_var, by = c("var_id" = "variant_id")) %>%
  distinct(finemapped_cs_coloc)

# NEW: caQTL credible sets that contain chromBPNet variants WITH motif_overlap == TRUE
chrombpnet_caqtlvars_motif <- summary_df %>%
  mutate(var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])",
                       "chr\\1_\\2_\\3_\\4",
                       variant_id.x)) %>%
  inner_join(cBPNet_var %>% dplyr::filter(motif_overlap),
             by = c("var_id" = "variant_id")) %>%
  distinct(finemapped_cs_coloc)

# (opposite) get the chromBPNet variants that are within caQTLs credible sets
bpnet_var2 <- summary_df %>%
  mutate(var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])",
                       "chr\\1_\\2_\\3_\\4",
                       variant_id.x)) %>%
  semi_join(cBPNet_var, by = c("var_id" = "variant_id")) %>%
  distinct(var_id)

# get eQTLs credible sets that contain chromBPNet variants
chrombpnet_eqtlvars <- summary_df %>%
  mutate(var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])",
                       "chr\\1_\\2_\\3_\\4",
                       variant_id.y)) %>%
  semi_join(cBPNet_var, by = c("var_id" = "variant_id")) %>%
  distinct(finemapped_cs_coloc)

# NEW: eQTL credible sets that contain chromBPNet variants WITH motif_overlap == TRUE
chrombpnet_eqtlvars_motif <- summary_df %>%
  mutate(var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])",
                       "chr\\1_\\2_\\3_\\4",
                       variant_id.y)) %>%
  inner_join(cBPNet_var %>% dplyr::filter(motif_overlap),
             by = c("var_id" = "variant_id")) %>%
  distinct(finemapped_cs_coloc)

# (opposite) get the chromBPNet variants that are within eQTLs credible sets
bpnet_var1 <- summary_df %>%
  mutate(var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])",
                       "chr\\1_\\2_\\3_\\4",
                       variant_id.y)) %>%
  semi_join(cBPNet_var, by = c("var_id" = "variant_id")) %>%
  distinct(var_id)

chrombpet_shared_vars <- unique(unlist(c(bpnet_var1, bpnet_var2))) # all the chrombpnet variants found in qtls

# add motif info to chromBPNet-only variants
cBPNet_var_only <- cBPNet_var %>%
  dplyr::filter(!(variant_id %in% chrombpet_shared_vars)) %>%
  dplyr::select(variant_id, motif_overlap) %>%   # <— keep motif here
  mutate(coloc_status = NA_character_,
         caqtl_category = NA_character_)

# all credible sets that contain a chromBPnet var (union)
unique_coloc_vars <- unique(unlist(c(chrombpnet_caqtlvars, chrombpnet_eqtlvars)))

# NEW: all credible sets that contain a chromBPNet var WITH motif overlap (union)
unique_coloc_vars_motif <- unique(unlist(c(chrombpnet_caqtlvars_motif,
                                           chrombpnet_eqtlvars_motif)))

print(paste0("The number of chromBPnet variants that are also QTLs is: ",
             print(length(unique_coloc_vars))))

unique(summary_df$coloc_status)



# build CS table as you already do
credible_sets <- summary_df %>%
  distinct(finemapped_cs_coloc, coloc_status, caqtl_category)
colnames(credible_sets) <- c("variant_id", "coloc_status", "caqtl_category")

credible_sets <- dplyr::bind_rows(credible_sets, cBPNet_var_only) %>%
  dplyr::mutate(
    is_variant_only   = is.na(coloc_status),
    chromBPNet        = dplyr::if_else(
      is_variant_only, TRUE, variant_id %in% unique_coloc_vars
    ),
    chromBPNet_motif  = dplyr::if_else(
      is_variant_only, dplyr::coalesce(motif_overlap, FALSE),
      variant_id %in% unique_coloc_vars_motif
    ),
    caqtl_category    = caqtl_category
  ) %>%
  dplyr::mutate(
    type = dplyr::case_when(
      coloc_status == "coloc"     & chromBPNet ~ "shared_all",
      coloc_status == "coloc"                    ~ "coloc",
      coloc_status == "caQTL_only" & chromBPNet ~ "shared_chrom_caqtl",
      coloc_status == "eQTL_only"  & chromBPNet ~ "shared_chrom_eqtl",
      coloc_status == "caQTL_only"              ~ "caQTL_only",
      coloc_status == "eQTL_only"               ~ "eQTL_only",
      TRUE                                      ~ "chromBPNet_only"
    )
  ) %>%
  dplyr::select(-is_variant_only, -motif_overlap)

rm(unique_coloc_vars,cBPNet_var_only, chrombpet_shared_vars, bpnet_var1, chrombpnet_eqtlvars,bpnet_var2,chrombpnet_caqtlvars)
rm(summary_df)
saveRDS(credible_sets, "~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/upset_credible_sets.rds")
credible_sets<-readRDS("~/cd4_qtl_paper_figures/figure_1/data/upset_credible_sets.rds")

##########################
#--------- UPSET PLOT to visualize variants shared across contexts
#########
# cat("plot sharing across context *VERSION1*") 
# 
# # Step 1: Assign binary membership and convert to logical
# membership_df <- credible_sets %>%
#   mutate(
#     eQTL = type %in% c("shared_all", "coloc", "shared_chrom_eqtl", "eQTL_only"),
#     caQTL = type %in% c("shared_all", "coloc", "shared_chrom_caqtl", "caQTL_only"),
#     chromBPNet = type %in% c("shared_all", "shared_chrom_eqtl", "shared_chrom_caqtl", "chromBPNet_only")
#   ) %>%
#   #select(-c(type, coloc_status))%>%
#   mutate(across(c(eQTL, caQTL, chromBPNet), as.logical))
# 
# # Define your sets
# set_names <- c("chromBPNet", "eQTL", "caQTL")
# 
# colors <- carto_pal( 7, "Temps")
# names(colors) <- set_names
# 
# # UpSet plot with LaCroix coloring
# #regular upset plot with grey bars
# upset_sharing<-upset(
#   membership_df,
#   intersect = set_names,
#   name = "Context",
#   width_ratio = 0.1,
#   wrap = TRUE,
#   n_intersections=9, 
#   base_annotations = list(
#     'Intersection size' = intersection_size(
#       mapping = aes(fill ="gray76")
#     ) +
#       scale_fill_manual(values = "gray76", guide = 'none')),
#   set_sizes=(
#     upset_set_size()
#     + geom_text(aes(label = paste0(round(after_stat(count)/1000), "K")), hjust = 1.1 , stat = 'count')
#     # + expand_limits(y=100000)
#     + theme(axis.text.x=element_blank())),
#   queries=list(
#     upset_query(set='caQTL', fill='#287274', color='#287274'),
#     upset_query(set='eQTL', fill='#C70E7B', color='#C70E7B'),
#     upset_query(set='chromBPNet',fill="#2CB11B", color="#2CB11B"))
# )
# 
# upset_sharing
# ggsave(paste0(out.dir, "upset_sharing_contexts.pdf"), upset_sharing, width = 7, height = 5)
# 
# 
# #########
# cat("plot sharing across context *VERSION2*") 
# #upset plot with bars colored depeding on caQTL categories
# # membership_df <- credible_sets %>%
# #   mutate(
# #     eQTL = type %in% c("shared_all", "coloc", "shared_chrom_eqtl", "eQTL_only"),
# #     caQTL = type %in% c("shared_all", "coloc", "shared_chrom_caqtl", "caQTL_only"),
# #     chromBPNet = type %in% c("shared_all", "shared_chrom_eqtl", "shared_chrom_caqtl", "chromBPNet_only")
# #   ) %>%
# #   select(variant_id, caqtl_category, eQTL, caQTL, chromBPNet) %>%
# #   mutate(across(c(eQTL, caQTL, chromBPNet), as.logical))
# 
# #######
# 
# #assign custom colors to each caQTL category
# 
# label_map <- c(
#   in_caPeak = "Within caQTL peak",
#   in_corr_Peak = "In peak correlated with caPeak",
#   in_other_Peak = "In other open peak",
#   no_Peak_overlap = "No peak overlap"
# )
# membership_df$caqtl_category <- factor(
#   membership_df$caqtl_category,
#   levels = names(label_map),
#   labels = label_map
# )
# category_colors <- carto_pal(5, "Mint")
# names(category_colors) <- c( "Within caQTL peak", "In peak correlated with caPeak", "In other open peak", "No peak overlap")
# 
# upset_sharing2 <- upset(
#   membership_df,
#   intersect = c("chromBPNet", "eQTL", "caQTL"),
#   name = "Context",
#   width_ratio = 0.1,
#   n_intersections = 7,
#   base_annotations = list(
#     'Intersection size' = intersection_size(
#       mapping = aes(fill = caqtl_category)
#     ) +
#       scale_fill_manual(
#         values = category_colors,
#         name = "caQTL category",
#         na.value = "gray76"
#       )
#   ),
#   set_sizes = (
#     upset_set_size() +  # <- remove 'position' argument!
#       geom_text(aes(label = paste0(round(after_stat(count)/1000), "K")), hjust = 1.1, stat = 'count') +
#       theme(axis.text.x = element_blank())
#   ),
#   queries = list(
#     upset_query(set = 'caQTL', fill = '#287274', color = '#287274'),
#     upset_query(set = 'eQTL', fill = '#C70E7B', color = '#C70E7B'),
#     upset_query(set = 'chromBPNet', fill = "#2CB11B", color = "#2CB11B")
#   )
# )
# 
# upset_sharing2
# ggsave(paste0(out.dir, "upset_sharing_contexts_painted_bars.pdf"), upset_sharing2, width = 7, height = 5)


#########
cat("plot sharing across context *VERSION3*") 
# upset plot with a bar graph stacked over it showing the proportions of the caQTL categories matching each group

membership_df <- credible_sets %>%
  mutate(
    eQTL = type %in% c("shared_all", "coloc", "shared_chrom_eqtl", "eQTL_only"),
    caQTL = type %in% c("shared_all", "coloc", "shared_chrom_caqtl", "caQTL_only"),
    chromBPNet = type %in% c("shared_all", "shared_chrom_eqtl", "shared_chrom_caqtl", "chromBPNet_only")
  ) %>%
  select(variant_id, caqtl_category, eQTL, caQTL, chromBPNet, chromBPNet_motif) %>%
  mutate(across(c(eQTL, caQTL, chromBPNet), as.logical))
#######

#assign custom colors to each caQTL category
label_map <- c(
  in_caPeak = "Within caQTL peak",
  in_other_Peak = "In other open peak",
  no_Peak_overlap = "No peak overlap"
)
membership_df$caqtl_category <- factor(
  membership_df$caqtl_category,
  levels = names(label_map),
  labels = label_map
)


category_colors <- carto_pal(5, "Mint")
names(category_colors) <- c("Within caQTL peak", "In other open peak", "No peak overlap")


# Add a bar annotation by caqtl_category
upset_sharing3<-upset(
  membership_df,
  intersect = c("chromBPNet", "eQTL", "caQTL"),
  name = "Context",
  width_ratio = 0.1,
  wrap = TRUE,
  base_annotations = list(
    'Intersection size' = intersection_size()
  ),
  annotations = list(
    'caQTL category'=(
      ggplot(mapping=aes(fill=caqtl_category))
      + geom_bar(stat='count', position='fill')
      + scale_y_continuous(labels=scales::percent_format())
      + scale_fill_manual(values = category_colors,
                          name = "caQTL category",
                          na.value = "gray89"
      ))
      + ylab('caQTL category'),
    'motif'=(
      ggplot(mapping=aes(fill=chromBPNet_motif))
      + geom_bar(stat='count', position='fill')
      + scale_y_continuous(labels=scales::percent_format())
      + scale_fill_manual(values = c("gray89", "red"),
                          name = "motif",
                          na.value = "gray89"
      ))
    + ylab('motif')
    ),
  set_sizes = (
    upset_set_size() +  # <- remove 'position' argument!
      geom_text(aes(label = paste0(round(after_stat(count)/1000), "K")), hjust = 1.1, stat = 'count') +
      theme(axis.text.x = element_blank())
  ),
  queries = list(
    upset_query(set = 'caQTL', fill = '#287274', color = '#287274'),
    upset_query(set = 'eQTL', fill = '#C70E7B', color = '#C70E7B'),
    upset_query(set = 'chromBPNet', fill = "#2CB11B", color = "#2CB11B")
  )
  )
upset_sharing3


#ggsave(paste0(out.dir, "upset_sharing_contexts_painted_bars_top.pdf"), upset_sharing3, width = 7, height = 5)


#####
#version masking "no peak overlap proportions"

## 0) Build explicit intersection labels used by BOTH panels
df_with_intersection <- membership_df %>%
  mutate(
    intersection = apply(select(., chromBPNet, eQTL, caQTL ), 1, function(v) {
      on <- names(v)[as.logical(v)]
      if (length(on) == 0) "None" else paste(on, collapse = " & ")
    })
  )

## 1) Choose EXACT intersections to display (and lock the order)
#    (edit n_keep to match what you want visible)
n_keep <- 7
keep_tbl <- df_with_intersection %>%
  count(intersection, name = "tot") %>%
  arrange(desc(tot)) %>%
  slice_head(n = n_keep)

keep_intersections <- keep_tbl$intersection

intersections_list <- lapply(keep_intersections, function(s) {
  parts <- strsplit(s, " & ", fixed = TRUE)[[1]]
  parts[parts != "None"]
})

## --- add this after you define `keep_intersections` and before building prop_plot ---

levels_inpeak <- c("Within caQTL peak",
                   "In other open peak")
category_colors_inpeak <- category_colors[levels_inpeak]

# totals per intersection (ALL rows, used as denominator in the label)
totals_all <- df_with_intersection %>%
  filter(intersection %in% keep_intersections) %>%
  count(intersection, name = "tot_total")

# proportions among in-peak only (as before)
anno_df <- df_with_intersection %>%
  filter(intersection %in% keep_intersections,
         caqtl_category %in% levels_inpeak) %>%           # drop "No peak overlap"
  count(intersection, caqtl_category, name = "n") %>%
  complete(intersection = keep_intersections,
           caqtl_category = factor(levels_inpeak, levels_inpeak),
           fill = list(n = 0)) %>%
  group_by(intersection) %>%
  mutate(
    tot_inpeak = sum(n),
    prop = ifelse(tot_inpeak > 0, n / tot_inpeak, NA_real_)
  ) %>%
  ungroup() %>%
  left_join(totals_all, by = "intersection") %>%
  mutate(
    intersection   = factor(intersection, levels = keep_intersections),
    caqtl_category = factor(caqtl_category, levels = levels_inpeak)
  )

# label data: “in-peak / total”
label_df <- anno_df %>%
  distinct(intersection, tot_inpeak, tot_total) %>%
  mutate(lab = sprintf("%s / %s", tot_inpeak, tot_total),
         y   = 0.9)

# gray bars for intersections with NO in-peak rows
na_bars <- label_df %>%
  filter(tot_inpeak == 0) %>%
  transmute(intersection, prop = 1)   # full-height gray bar

## --- TOP panel ---
prop_plot <-
  ggplot() +
  # draw gray full bars first for NA intersections
  geom_col(
    data = na_bars,
    aes(x = intersection, y = prop),
    width = 0.9,
    fill = "grey85",
    inherit.aes = FALSE,
    linewidth=0.1, color="gray30"
  ) +
  # draw colored stacks for finite proportions (in-peak only)
  geom_col(
    data = filter(anno_df, is.finite(prop)),
    aes(x = intersection, y = prop, fill = caqtl_category),
    width = 0.9,
    inherit.aes = FALSE
  ) +
  # totals label (in-peak / total)
  geom_text(
    data = label_df,
    aes(x = intersection, y = y, label = lab),
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_x_discrete(limits = keep_intersections, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0,0)) +
  scale_fill_manual(values = category_colors_inpeak,
                    name   = "caQTL category (in-peak)",
                    drop   = TRUE) +
  labs(x = NULL, y = "Proportion (in-peak only)") +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(5, 10, 0, 10)
  )


## 4) BOTTOM panel (ComplexUpset) using the SAME intersections/order
up_base <- upset(
  df_with_intersection,
  intersect            = c("chromBPNet", "eQTL", "caQTL"),
  intersections        = intersections_list,   # lock set & order
  sort_intersections   = FALSE,                # don’t re-sort internally
  wrap                 = TRUE,
  width_ratio          = 0.1,
  base_annotations = list(
    "Intersection size" = intersection_size()
  ),
  annotations = list(),
    # 'motif'=(
    #   ggplot(mapping=aes(fill=chromBPNet_motif))
    #   + geom_bar(stat='count', position='fill')
    #   + scale_y_continuous(labels=scales::percent_format())
    #   + scale_fill_manual(values = c("gray89", "red"),
    #                       name = "motif",
    #                       na.value = "gray89"
    #   ))
    # + ylab('motif')
  #),
  set_sizes = (
    upset_set_size() +
      geom_text(aes(label = paste0(round(after_stat(count)/1000), "K")),
                hjust = 1.1, stat = "count") +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ),
  queries = list(
    upset_query(set = "caQTL",      fill = "#287274", color = "#287274"),
    upset_query(set = "eQTL",       fill = "#C70E7B", color = "#C70E7B"),
    upset_query(set = "chromBPNet", fill = "#2CB11B", color = "#2CB11B")
  )
)

## 5) Stack with patchwork (aligned)
upset_sharing4<-prop_plot / up_base + plot_layout(heights = c(1, 4), guides = "collect")
upset_sharing4

#ggsave(paste0(out.dir, "upset_sharing_contexts_painted_bars_top_cond_inpeak.pdf"), upset_sharing4, width = 8, height = 5)

# 
# ######################################################
# ######################################################
# ######################################################
# ##      GET VARIANTS IN EACH INTERCEPTION 
# #         TO PLOT RELATIVE TSS DISTANCE
# ######################################################
# ######################################################
# ######################################################
# # Create a named list for all 7 intersections
# variants_by_intersection <- list(
#   
#   # 1) chromBPNet only
#   chromBPNet_only = membership_df %>%
#     filter(chromBPNet & !eQTL & !caQTL) %>%
#     pull(variant_id),
#   
#   # 2) eQTL only
#   eQTL_only = membership_df %>%
#     filter(eQTL & !caQTL & !chromBPNet) %>%
#     pull(variant_id),
#   
#   # 3) caQTL only
#   caQTL_only = membership_df %>%
#     filter(caQTL & !eQTL & !chromBPNet) %>%
#     pull(variant_id),
#   
#   # 4) chromBPNet & eQTL
#   chromBPNet_eQTL = membership_df %>%
#     filter(chromBPNet & eQTL & !caQTL) %>%
#     pull(variant_id),
#   
#   # 5) chromBPNet & caQTL
#   chromBPNet_caQTL = membership_df %>%
#     filter(chromBPNet & caQTL & !eQTL) %>%
#     pull(variant_id),
#   
#   # 6) eQTL & caQTL
#   eQTL_caQTL = membership_df %>%
#     filter(eQTL & caQTL & !chromBPNet) %>%
#     pull(variant_id),
#   
#   # 7) all three
#   all_three = membership_df %>%
#     filter(chromBPNet & eQTL & caQTL) %>%
#     pull(variant_id)
# )
# 
# # Check the names and counts
# sapply(variants_by_intersection, length)
# 
# head(summary_df)
# 
# #get variant_ids
# chromBPNet_only_vars<-variants_by_intersection$chromBPNet_only
# 
# eQTL_only_vars <- summary_df %>%
#   filter(finemapped_cs_coloc %in% variants_by_intersection$eQTL_only) %>%
#   select(finemapped_cs_coloc, variant_id.y) %>%
#   mutate(var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])", 
#                        "chr\\1_\\2_\\3_\\4", 
#                        variant_id.y))
# 
# caQTL_only_vars <- summary_df %>%
#   filter(finemapped_cs_coloc %in% variants_by_intersection$caQTL_only) %>%
#   select(finemapped_cs_coloc, variant_id.x)  %>%
#   mutate(var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])", 
#                        "chr\\1_\\2_\\3_\\4", 
#                        variant_id.x))
# 
# chromBPNet_eQTL_vars <- summary_df %>%
#   filter(finemapped_cs_coloc %in% variants_by_intersection$chromBPNet_eQTL) %>%
#   select(finemapped_cs_coloc, variant_id.y) %>%
#   mutate(var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])", 
#                        "chr\\1_\\2_\\3_\\4", 
#                        variant_id.y))
# 
# chromBPNet_caQTL_vars <- summary_df %>%
#   filter(finemapped_cs_coloc %in% variants_by_intersection$chromBPNet_caQTL) %>%
#   select(finemapped_cs_coloc, variant_id.x) %>%
#   mutate(var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])", 
#                        "chr\\1_\\2_\\3_\\4", 
#                        variant_id.x))
# 
# eQTL_caQTL_vars <- summary_df %>%
#   filter(finemapped_cs_coloc %in% variants_by_intersection$eQTL_caQTL) %>%
#   select(finemapped_cs_coloc, variant_id.x, variant_id.y) %>%
#   pivot_longer(cols = c(variant_id.x, variant_id.y),
#                names_to = "source",
#                values_to = "variant_id") %>%
#   filter(!is.na(variant_id)) %>%
#   distinct(finemapped_cs_coloc, variant_id) %>%
#   mutate(var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])", 
#                        "chr\\1_\\2_\\3_\\4", 
#                        variant_id))
# 
# all_three_vars <- summary_df %>%
#   filter(finemapped_cs_coloc %in% variants_by_intersection$all_three) %>%
#   select(finemapped_cs_coloc, variant_id.x, variant_id.y) %>%
#   pivot_longer(cols = c(variant_id.x, variant_id.y),
#                names_to = "source",
#                values_to = "variant_id") %>%
#   filter(!is.na(variant_id)) %>%
#   distinct(finemapped_cs_coloc, variant_id) %>%
#   mutate(var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])", 
#                        "chr\\1_\\2_\\3_\\4", 
#                        variant_id))
# 
# #concatenate all tables
# eQTL_only_vars        <- eQTL_only_vars        %>% select(finemapped_cs_coloc, var_id) %>% distinct()
# caQTL_only_vars       <- caQTL_only_vars       %>% select(finemapped_cs_coloc, var_id) %>% distinct()
# chromBPNet_eQTL_vars  <- chromBPNet_eQTL_vars  %>% select(finemapped_cs_coloc, var_id) %>% distinct()
# chromBPNet_caQTL_vars <- chromBPNet_caQTL_vars %>% select(finemapped_cs_coloc, var_id) %>% distinct()
# eQTL_caQTL_vars       <- eQTL_caQTL_vars       %>% select(finemapped_cs_coloc, var_id) %>% distinct()
# all_three_vars        <- all_three_vars        %>% select(finemapped_cs_coloc, var_id) %>% distinct()
# 
# # chromBPNet_only is just a vector; make it a tibble (no CS/coloc id for these)
# chromBPNet_only_tbl <- tibble(
#   finemapped_cs_coloc = NA_character_,
#   var_id = chromBPNet_only_vars
# )
# 
# # concatenate + tag source group
# all_vars <- bind_rows(
#   mutate(eQTL_only_vars,        group = "eQTL_only"),
#   mutate(caQTL_only_vars,       group = "caQTL_only"),
#   mutate(chromBPNet_eQTL_vars,  group = "chromBPNet_eQTL"),
#   mutate(chromBPNet_caQTL_vars, group = "chromBPNet_caQTL"),
#   mutate(eQTL_caQTL_vars,       group = "eQTL_caQTL"),
#   mutate(all_three_vars,        group = "all_three"),
#   mutate(chromBPNet_only_tbl,   group = "chromBPNet_only")
# ) %>%
#   relocate(group, finemapped_cs_coloc, var_id) %>%
#   distinct()
# ######################################################
# ######################################################
# ######################################################


# coord_pat <- "(?:[0-9]+|X|Y|M):\\d+-\\d+"
# 
# all_vars <- all_vars %>%
#   mutate(
#     # Peak label = everything before the first genomic interval
#     peak_label = if_else(
#       is.na(finemapped_cs_coloc),
#       NA_character_,
#       str_replace(finemapped_cs_coloc, paste0("_(?:", coord_pat, ").*"), "")
#     ),
#     # Gene when two intervals (shared): token between them
#     gene_between = str_match(
#       finemapped_cs_coloc,
#       paste0(".*?", coord_pat, "_[^_]+_([^_]+)_", coord_pat)
#     )[, 2],
#     # Gene when a single interval and leading token is the gene (e.g. PFKP_10:...)
#     gene_alt = str_match(
#       finemapped_cs_coloc,
#       paste0("^([^_]+)_(?:", coord_pat, ")")
#     )[, 2],
#     gene_label = coalesce(gene_between, gene_alt),
#     
#     # === Enforce your mapping ===
#     feature1 = case_when(
#       group %in% c("eQTL_only", "chromBPNet_eQTL", "eQTL_caQTL", "all_three") ~ gene_label,
#       TRUE ~ NA_character_
#     ),
#     feature2 = case_when(
#       group %in% c("caQTL_only", "chromBPNet_caQTL", "eQTL_caQTL", "all_three") ~ peak_label,
#       TRUE ~ NA_character_
#     )
#   ) %>%
#   select(group, finemapped_cs_coloc, var_id, feature1, feature2)

# write_csv(all_vars, "~/cd4_qtl_paper_figures/figure_1/data/bpnet_qtl_upset_interception_comb_finemap_vas.csv")

# cat("Finished!")

all_vars<-read_csv("~/cd4_qtl_paper_figures/figure_1/data/bpnet_qtl_upset_interception_comb_finemap_vas.csv")

# Normalizer: supports both "chr10_3126444_C_T" and "10:5913362[b38]AC,A"
normalize_variant_id <- function(x) {
  x <- trimws(x)
  out <- rep(NA_character_, length(x))

  # Case 1: chr10_3126444_C_T   (allow optional 'chr', any case)
  m1 <- str_match(x, "^(?:chr)?([0-9XYM]{1,2})_(\\d+)_([ACGTacgt]+)_([ACGTacgt]+)$")
  keep1 <- !is.na(m1[,1])
  out[keep1] <- paste0("chr", m1[keep1,2], "_", m1[keep1,3], "_",
                       toupper(m1[keep1,4]), "_", toupper(m1[keep1,5]))

  # Case 2: 10:5913362[b38]AC,A (allow optional 'chr', optional [b..])
  m2 <- str_match(x, "^(?:chr)?([0-9XYM]{1,2}):(\\d+)(?:\\[b\\d+\\])?([ACGTacgt]+),([ACGTacgt]+)$")
  keep2 <- !is.na(m2[,1]) & is.na(out)
  out[keep2] <- paste0("chr", m2[keep2,2], "_", m2[keep2,3], "_",
                       toupper(m2[keep2,4]), "_", toupper(m2[keep2,5]))

  out
}

# 2) Standardize and parse
all_vars <- all_vars %>%
  mutate(
    var_id_std = normalize_variant_id(var_id)
  ) %>% mutate(
    finemapped_cs_coloc=ifelse(is.na(finemapped_cs_coloc),var_id_std,  finemapped_cs_coloc)
  )

# #read the gencode gene annotation
# gencode <- fread("~/resources/genome/hg38_gencode_raw/gencode.v44.primary_assembly.annotation.genes.bed",
#                  sep="\t", data.table=FALSE) %>%
#   #filter(V8 == "protein_coding",
#    #      grepl("^chr(\\d+|X|Y)$", V1)) %>%                 # drop alts/patches
#   transmute(
#     chr        = V1,
#     strand     = V6,
#     gene_id    = V7,
#     gene       = V4,                                       # gene symbol
#     tss_1based = if_else(V6 == "+", V2 + 1L, V3)           # BED start is 0-based
#   ) %>%
#   distinct(chr, gene, .keep_all = TRUE)                    # ensure 1 row per chr×gene
# 
# # 2) Groups of interest
# groups_interest <- c("eQTL_only", "chromBPNet_eQTL", "eQTL_caQTL", "all_three")
# exp <- all_vars %>% filter(group %in% groups_interest)
# 
# # 3) Parse chr/pos from var_id
# exp_parsed <- exp %>%
#   mutate(
#     chr_chr  = str_match(var_id_std, "^(chr[0-9XYM]+)_")[,2],
#     pos_chr  = str_match(var_id_std, "^chr[0-9XYM]+_([0-9]+)_")[,2],
#     chr_b38  = str_match(var_id_std, "^([0-9XYM]+):")[,2],
#     pos_b38  = str_match(var_id_std, "^[0-9XYM]+:([0-9]+)\\[b38\\]")[,2],
#     chr = coalesce(chr_chr, if_else(!is.na(chr_b38), paste0("chr", chr_b38), NA_character_)),
#     pos = coalesce(as.integer(pos_chr), as.integer(pos_b38))
#   ) %>%
#   select(-chr_chr, -pos_chr, -chr_b38, -pos_b38) %>%
#   filter(!is.na(chr), !is.na(pos), !is.na(feature1))  # feature1 is the gene
# 
# 
# # 4) Join to gene TSS (by gene symbol + chromosome) and compute distances
# exp_with_tss <- exp_parsed %>%
#   left_join(gencode, by = c("feature1" = "gene", "chr" = "chr")) %>%
#   mutate(
#     distance_to_TSS = abs(pos - tss_1based),
#     signed_distance = case_when(
#       strand == "+" ~ pos - tss_1based,     # downstream = positive
#       strand == "-" ~ tss_1based - pos,     # downstream (5'→3') still positive
#       TRUE           ~ NA_integer_
#     )
#   ) 
# 
# exp_summary <- exp_with_tss %>%
#   distinct(group, finemapped_cs_coloc, feature1, var_id, .keep_all = TRUE) %>%
#   group_by(group, finemapped_cs_coloc, feature1) %>%
#   summarise(mean_tss = mean(signed_distance, na.rm = TRUE),   # use signed_distance!
#             median_tss = median(signed_distance, na.rm = TRUE),
#             .groups = "drop") %>%
#   mutate(
#     abs_dist = abs(mean_tss),   
#     TSS_bin = cut(
#       abs_dist,
#       breaks = c(0, 1e3, 1e4, 5e4, 1e5, 5e5, Inf),
#       labels = c("<1 kb", "1–10 kb", "10–50 kb", "50–100 kb", "100–500 kb", ">500 kb"),
#       include.lowest = TRUE, right = FALSE
#     ),
#     # lock factor levels so they stay even if empty
#     TSS_bin = factor(
#       TSS_bin,
#       levels = c("<1 kb", "1–10 kb", "10–50 kb", "50–100 kb", "100–500 kb", ">500 kb")
#     )
#   )
# 
# 
# colors <- carto_pal(7, "PinkYl")
# 
# bplot_tss <- ggplot(exp_summary, aes(x = group, fill = TSS_bin)) +
#   geom_bar(position = "fill") +
#   scale_fill_manual(values = colors, name = "TSS distance") +
#   labs(x = "Group", y = "Prop. variants") +
#   theme_minimal()
# 
# bplot_tss
# ggsave(paste0(out.dir, "upset_sharing_contexts_bplot_tss.pdf"), bplot_tss, width = 8, height = 5)

## Map your groups → the same intersection labels used by ComplexUpset
# group_to_intersection <- c(
#   "eQTL_only"        = "eQTL",
#   "chromBPNet_eQTL"  = "chromBPNet & eQTL",
#   "eQTL_caQTL"       = "eQTL & caQTL",
#   "all_three"        = "chromBPNet & eQTL & caQTL"
# )
# 
# # palette (4 groups you’re plotting)
# group_levels <- c("eQTL_only","chromBPNet_eQTL","eQTL_caQTL","all_three")
# 
# 
# # Map groups -> intersections used in your UpSet
# exp_summary2 <- exp_summary %>%
#   mutate(
#     intersection = recode(group, !!!group_to_intersection),
#     intersection = factor(intersection, levels = keep_intersections)
#   )
# 
# # Count per (intersection × TSS_bin); keep empty combinations
# counts <- exp_summary2 %>%
#   mutate(TSS_bin = fct_na_value_to_level(TSS_bin, NA)) %>%  # <-- this line
#   count(intersection, TSS_bin, name = "n") %>%
#   group_by(intersection) %>%
#   mutate(total_n  = sum(n, na.rm = TRUE)) %>%
#   ungroup()
# 
# # Finite (non-NA) stacks: proportions per TSS_bin within intersection
# anno_df2 <- counts %>%
#   filter(!is.na(TSS_bin)) %>%
#   group_by(intersection, TSS_bin, total_n) %>%
#   summarise(n = sum(n), .groups = "drop") %>%
#   mutate(prop = if_else(total_n > 0, n / total_n, 0))
# 
# # Gray background bars: full height = 1 (so overlay shows in-range fraction)
# na_bars2 <- counts %>%
#   distinct(intersection) %>%
#   mutate(prop = 1)
# 
# 
# # Colors for your TSS bins (must match TSS_bin levels you set earlier)
# tss_levels <- c("<1 kb","1–10 kb","10–50 kb","50–100 kb","100–500 kb",">500 kb")
# category_colors_inrange <- setNames(carto_pal(length(tss_levels), "BurgYl"), tss_levels)
# 
# barplot_tss <-
#   ggplot() +
#   # gray full bars
#   geom_col(
#     data = na_bars2,
#     aes(x = intersection, y = prop),
#     width = 0.9,
#     fill = "grey85",
#     inherit.aes = FALSE
#   ) +
#   # colored overlay (finite/in-range only)
#   geom_col(
#     data = anno_df2,
#     aes(x = intersection, y = prop, fill = TSS_bin),
#     width = 0.9,
#     inherit.aes = FALSE
#   ) +
#   scale_x_discrete(limits = keep_intersections, drop = FALSE) +
#   scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0,0)) +
#   scale_fill_manual(
#     values = category_colors_inrange,
#     name   = "Distance to TSS",
#     drop   = FALSE
#   ) +
#   labs(x = NULL, y = "Proportion (in-range bins)") +
#   theme_minimal() +
#   theme(
#     legend.position = "right",
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     panel.grid.major.x = element_blank(),
#     plot.margin = margin(5, 10, 0, 10)
#   )
# 
# barplot_tss
# 
# ## Stack with your existing panels
# combined <- barplot_tss / prop_plot / up_base +
#   plot_layout(heights = c(1, 1, 4), guides = "collect")
# 
# combined

#ggsave(paste0(out.dir, "upset_sharing_contexts_painted_bars_top_cond_inpeak_tss.pdf"), combined, width = 8, height = 5)

#############
# Convert "chr10_3126444_C_T"  -->  "10:3126444[b38]C,T"

to_b38_variant_id <- function(x, build = "b38") {
  m <- str_match(
    x,
    regex("^(?:chr)?([0-9]{1,2}|X|Y|M|MT)[_:-](\\d+)[_:-]([^_:-]+)[_:-]([^_:-]+)$",
          ignore_case = TRUE)
  )
  chr <- toupper(m[, 2])
  pos <- m[, 3]
  ref <- toupper(m[, 4])    # accepts multi-base alleles
  alt <- toupper(m[, 5])
  
  # normalize M -> MT
  chr[chr == "M"] <- "MT"
  
  ifelse(is.na(chr), NA_character_, paste0(chr, ":", pos, "[", build, "]", ref, ",", alt))
}


# Example on your tibble:
all_vars <- all_vars %>%
  dplyr::mutate(variant_id.plink = to_b38_variant_id(var_id_std))

#read plink calculated MAFs
mafs<-fread("/gpfs/commons/home/mmatos/cd4_QTL_analysis/01_genotype_snps_covar/02_genotype_covariates/analysis/002.v4_calculating_pcs_MAF5/CD4_all_chr_ashkenazi.363.AF1.QC.BA.king2.hwe.scrnaseq.annot.freq.frq")
#hist(mafs$MAF, main = "Allele Frequency Distributions of Merged Data", xlab = "Minor Allele Frequency", xlim = c(0.0, 0.5), col = "plum2")

all_vars<-all_vars%>% left_join(mafs, by = c("variant_id.plink" = "SNP")) 
summary(is.na(all_vars$MAF))

group_to_intersection <- c("caQTL_only"= "caQTL",
"eQTL_caQTL" = "eQTL & caQTL",
"chromBPNet_only"= "chromBPNet",
"all_three" = "chromBPNet & eQTL & caQTL",
"chromBPNet_caQTL"= "chromBPNet & caQTL", 
"eQTL_only"= "eQTL",
"chromBPNet_eQTL"= "chromBPNet & eQTL"  )


# Map groups -> intersections used in your UpSet
all_vars <- all_vars %>%
  mutate(
    intersection = recode(group, !!!group_to_intersection),
    intersection = factor(intersection, levels = keep_intersections)
  )


# 1) Deduplicate: one row per (intersection, var_id_std)
dedup <- all_vars %>%
  filter(!is.na(MAF)) %>%
  group_by(intersection, var_id_std) %>%
  summarise(MAF = first(MAF), .groups = "drop")   # or slice_max(pip, n=1) if you prefer

# 2) Bin globally (use cut_number to avoid duplicate quantile breaks)
barplot_stacked <- dedup %>%
  mutate(
    MAF_bin = cut_number(MAF, n = 5),   # same count as your 20% quantiles
    MAF_bin = fct_rev(MAF_bin)
  ) %>%
  ggplot(aes(x = intersection, fill = MAF_bin)) +
  geom_bar(position = "fill", linewidth = 0.1, color = "gray30") +
  scale_fill_brewer(palette = "GnBu", name = "MAF bin", direction = -1, drop = FALSE) +
  labs(x = "Group", y = "Prop. variants in MAF bins") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank())

barplot_stacked


upset_sharing4

#ggsave(paste0(out.dir, "upset_sharing_contexts_interceptions_mafs.pdf"), upset_sharing4, width = 8, height = 8)


###### chromatin regulatory annotations 
cre_annot<-readRDS("~/cd4_qtl_paper_figures/figure_1/data/rel_tss_distance_CRE_peak_annotation.rds")
# make GRanges
cre_annot_gr <- GRanges(
  seqnames = cre_annot$chr,
  ranges   = IRanges(as.integer(cre_annot$start), as.integer(cre_annot$end)),
  peak_name = cre_annot$peak_name,
  annotation = cre_annot$annotation
)

library(dplyr)
library(stringr)
library(GenomicRanges)
library(S4Vectors)

# 0) Add a per-row unique id to the source df
all_vars <- all_vars %>% mutate(row_id = row_number())

# 1) Parse standardized IDs
parsed <- str_match(all_vars$var_id_std, "^(chr[0-9XYM]+)_(\\d+)_([ACGT]+)_([ACGT]+)$")
ok <- !is.na(parsed[,1])

# 2) Build variant GRanges only for parsed rows, **carry row_id**
var_in_interceptions <- GRanges(
  seqnames = parsed[ok, 2],
  ranges   = IRanges(as.integer(parsed[ok, 3]), as.integer(parsed[ok, 3])),
  var_id   = all_vars$var_id_std[ok],
  row_id   = all_vars$row_id[ok]
)

# 3) Overlaps (first hit per query; change to collapse if you want all)
ov  <- findOverlaps(var_in_interceptions, cre_annot_gr, ignore.strand = TRUE)
idx <- selectHits(ov, select = "first")

q_idx <-queryHits(ov)
s_idx <- subjectHits(ov)

# 4) Annotation vector aligned to var_in_interceptions
annotation_over <- rep(NA_character_, length(var_in_interceptions))
if (length(idx)) {
  annotation_over[q_idx] <- as.character(mcols(cre_annot_gr)$annotation[s_idx])
}
mcols(var_in_interceptions)$annotation_cre <- annotation_over

# 5) Build a 1:1 mapping by row_id (one row per original row)
map_df <- tibble(
  row_id        = mcols(var_in_interceptions)$row_id,
  annotation_cre = mcols(var_in_interceptions)$annotation_cre
)

# If you want to be extra-safe, enforce one row per row_id:
map_df <- map_df %>% group_by(row_id) %>%
  summarise(annotation_cre = dplyr::first(annotation_cre), .groups = "drop")

# 6) Join back by row_id (guaranteed many-to-one)
all_vars_annot <- all_vars %>%
  left_join(map_df, by = "row_id", relationship = "many-to-one")

# Sanity checks
stopifnot(nrow(all_vars_annot) == nrow(all_vars))
table(is.na(all_vars_annot$annotation_cre))

all_vars_annot <- all_vars_annot %>%
  mutate(
    anno_variant = case_when(
      !is.na(annotation_cre) & str_detect(annotation_cre, "(?i)promoter") ~ "promoter",
      !is.na(annotation_cre) & str_detect(annotation_cre, "(?i)distal")   ~ "distal",
      TRUE ~ NA_character_
    )
  ) 

# 1) Collapse to one annotation per credible set using the hierarchy
cs_annot <- all_vars_annot %>%
  group_by(finemapped_cs_coloc) %>%
  summarise(
    anno_cs = case_when(
      any(anno_variant == "promoter", na.rm = TRUE) ~ "promoter",
      any(anno_variant == "distal",   na.rm = TRUE) ~ "distal",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    anno_cs = forcats::fct_explicit_na(
      factor(anno_cs, levels = c("promoter", "distal")),
      na_level = "NA"
    )
  )

cs_annot_majority <- all_vars_annot %>%
  filter(!is.na(anno_variant)) %>%                             # ignore NA votes
  count(finemapped_cs_coloc, anno_variant, name = "n") %>%     # tally per CS × annotation
  group_by(finemapped_cs_coloc) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%                   # pick the majority
  ungroup() %>%
  transmute(
    finemapped_cs_coloc,
    anno_cs_majority = factor(anno_variant, levels = c("promoter", "distal"))
  ) %>%
  right_join(
    all_vars_annot %>% distinct(finemapped_cs_coloc),
    by = "finemapped_cs_coloc"
  ) %>%
  mutate(anno_cs_majority = fct_explicit_na(anno_cs_majority, na_level = "NA"))

# 2) Attach back to variants
all_vars_annot <- all_vars_annot %>%
  left_join(cs_annot, by = "finemapped_cs_coloc")

all_vars_annot <- all_vars_annot %>%
  left_join(cs_annot_majority, by = "finemapped_cs_coloc")


all_vars_annot<-all_vars_annot %>%
  mutate(
    intersection = recode(group, !!!group_to_intersection),
    intersection = factor(intersection, levels = keep_intersections))

cre_annot_barplot<-all_vars_annot %>%
  group_by(intersection, anno_cs) %>%
  summarise(count = n(), .groups = "drop") %>%
  ggplot(aes(x = intersection, y = count, fill = anno_cs)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    y = "Proportion of variants",
    fill = "CRE annotation"
  ) +
  theme_minimal() +
  theme(#axis.text.x = element_blank(),
        axis.title.x =element_blank() )

cre_annot_barplot2<-all_vars_annot %>%
  group_by(intersection, anno_cs_majority) %>%
  summarise(count = n(), .groups = "drop") %>%
  ggplot(aes(x = intersection, y = count, fill = anno_cs_majority)) +
  geom_col(position = "fill", linewidth = 0.1, color = "gray30") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_brewer(palette="BrBG", direction = -1) +
  labs(
    y = "Proportion of variants",
    fill = "CRE annotation"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x =element_blank() )

plot_grid(cre_annot_barplot, cre_annot_barplot2, ncol=1)

###########

upset_sharing4 <- barplot_stacked /cre_annot_barplot2 / prop_plot / up_base +
  plot_layout(
    heights = c(1, 1, 1, 4),
    guides  = "collect"
  )
upset_sharing4

ggsave(paste0(out.dir, "upset_sharing_contexts_interceptions_mafs.pdf"), upset_sharing4, width = 8, height = 8)




# ##########################
# #--------- Plot categories caQTL
# cat("Plot categories caQTL") 
# 
# colors <- carto_pal(5, "Mint")
# names(colors) <- c("in_caPeak", "in_corr_Peak", "in_other_Peak", "no_Peak_overlap")
# label_map <- c(
#   in_caPeak = "Within caQTL peak",
#   in_corr_Peak = "In peak correlated with caPeak",
#   in_other_Peak = "In other open peak",
#   no_Peak_overlap = "No peak overlap"
# )
# 
# # Plot with custom fill scale
# 
# stack_plot<-summary_df %>%
#   filter(!is.na(finemapped_cs_caqtl)) %>%
#   group_by(caqtl_category) %>%
#   summarise(n = n_distinct(finemapped_cs)) %>%
#   mutate(unit = "caQTLs") %>% 
#   ggplot(aes(x = unit, y = n, fill = caqtl_category)) +
#   geom_bar(stat = "identity", position = "stack") +
#   geom_text(
#     aes(label = n),
#     position = position_stack(vjust = 0.5),
#     color = "black",
#     size = 3
#   ) +
#   scale_fill_manual(
#     values = colors,
#     labels = label_map,
#     name = "Category"
#   ) +
#   labs(
#     title = "caQTL Categories",
#     x = NULL,
#     y = "Count"
#   ) +
#   theme_minimal()
# 
# ggsave(paste0(out.dir,"peaks_vs_cs_plot_stacked.pdf"), stack_plot, width = 5, height = 5)
# 
# ##########################
# #--------- Plot relative distance to caPeak
# cat("Plot relative distance to caPeak") 
# #Read caQTL finemapping
# caqtl_finemapping <- as.data.frame(read_delim("~/cd4_qtl_paper_figures/figure_1/data/CD4T_combined_finemapping_annotated.csv"))
# 
# peak_count_path <- "/gpfs/commons/home/mmatos/ATAC-seq_analysis/diff_accesibility_ana/results/peak_counts/RAW_cd4_atac_peakcounts_ranges_scrna_union.csv"
# 
# # Load peak coordinates as GRanges once
# ranges.table <- read.csv2(peak_count_path, sep = ",") %>%
#   filter(Chr %in% paste0("chr", 1:22))
# rownames(ranges.table) <- ranges.table$X
# gr <- GRanges(seqnames = ranges.table$Chr,
#               ranges = IRanges(ranges.table$Start, ranges.table$End),
#               mcols = data.frame(peakID = rownames(ranges.table)))
# names(gr) <- gr$mcols.peakID
# 
# 
# gr$peakID <- names(gr)
# gr
# 
# caqtl_finemapping <- caqtl_finemapping %>%
#   mutate(
#     var_chr = str_extract(variant_id, "^[^:]+"),
#     var_pos = as.numeric(str_extract(variant_id, "(?<=:)\\d+"))
#   )
# 
# # Join data by peakID
# # Merge to get matched peaks from GRanges gr
# matched_peaks <- caqtl_finemapping %>%
#   filter(!is.na(peak)) %>%
#   left_join(as.data.frame(gr), by = c("peak" = "peakID"))
# 
# # Compute distances directly between variant and matched peak center
# matched_peaks <- matched_peaks %>%
#   mutate(
#     peak_center = start + (end - start) / 2,
#     distance = var_pos - peak_center
#   )
# 
# binned_matched <- matched_peaks %>%
#   mutate(bin = cut(distance, breaks = seq(-50000, 50000, by = 1000))) %>%
#   count(bin, caqtl_category) %>%
#   mutate(bin_center = as.numeric(sub("\\((.+),.*", "\\1", bin)))
# 
# # Normalize by number of variants per category
# n_variants_matched <- matched_peaks %>%
#   filter(!is.na(caqtl_category)) %>%
#   count(caqtl_category, name = "n_variants")
# 
# binned_norm_matched <- binned_matched %>%
#   left_join(n_variants_matched, by = "caqtl_category") %>%
#   mutate(norm_peak_count = n / n_variants)
# 
# colors <- carto_pal(7, "BluYl")
# colors <- colors[2:6]
# 
# names(colors) <- c("in_caPeak", "in_corr_Peak", "in_other_Peak", "no_Peak_overlap")
# 
# var_distance<-ggplot(binned_norm_matched, aes(x = bin_center, y = norm_peak_count, color = caqtl_category)) +
#   geom_line(linewidth = 1) +
#   scale_x_continuous(labels = function(x) paste0(x / 1000, "kb")) +
#   labs(
#     x = "Relative distance of caQTL to caPeak (kb)",
#     y = "Normalized count per variant",
#     color = "caQTL category"
#   ) +  ylim(0, 0.05)+
#   theme_minimal() +
#   scale_color_manual(
#     values = colors,
#     name = "Category"
#   )
# var_distance
# 
# ggsave(paste0(out.dir,"var_distance_caqtl_cat_50k.pdf"), var_distance, width = 4, height = 3)




########################################
### statistical testing ####
########################################
membership_df2 <- membership_df %>%
  mutate(bpnet_specific = chromBPNet & !eQTL & !caQTL)

#caQTL-eQTL and caQTL-eQTL-ChromBPNet intersections, there were higher proportions of C2 caQTLs than for the other categories
#Within caQTL∩eQTL loci, is chromBPNet membership associated with a higher probability of being C2?
c2_label <- "in_other_Peak" 

df2 <- membership_df2 %>%
  filter(caQTL, eQTL) %>%
  mutate(group = ifelse(chromBPNet, "triple", "double_only"),
         is_C2 = (caqtl_category == c2_label))

tab <- table(df2$group, df2$is_C2)
tab
res<-fisher.test(tab)


#############################
#chrombpnet specific
########################################
cbpnet_specific<-all_vars_annot[all_vars_annot$intersection=="chromBPNet",]
tb <- table(cbpnet_specific$anno_cs_majority, useNA = "ifany")

# % including NA
pct_all <- round(100 * prop.table(tb), 2)

sprintf(
  "The percentage of ChromBPNet-specific variants that do not fall in any annotated caPeak region is %.2f%%",
  pct_all[[3]]
)

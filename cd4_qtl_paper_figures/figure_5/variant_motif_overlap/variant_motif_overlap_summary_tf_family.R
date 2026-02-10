library(ggplot2)
library(readr)
library(dplyr)
library(cowplot)
library(GenomicRanges)
library(scales)
library(tidyr)
library(ggseqlogo)   # geom_logo(), theme_logo()
library(purrr)       # map()
library(stringr)     # file-name helpers
library(tibble)
library(forcats)
library(RColorBrewer)
library(forcats)

options(bitmapType = "cairo")
source("~/cd4_qtl_paper_figures/utils/color_pallete_helper.R")

# load hit counts per motif 
unique_hits <- read_tsv("~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap/motifs_hits_filtered_c90_40seq_w_variant_hits_boolean_1b.tsv")
vars <- read_tsv("~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap/cd4_top_cpbnet_variants_motif_overlap_moreTFS_boolean_1b.tsv")
variant_motif_map<- read_tsv("~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap/motif_variant_overlap_all_1b.tsv")

unique_hits %>% ggplot(aes(hit_coefficient)) + geom_histogram()
unique_hits %>% filter(hit_coefficient>5) %>% ggplot(aes(hit_coefficient)) + geom_histogram()
nrow(unique_hits %>% filter(hit_coefficient>5) )

# Step 1: Compute order of TFs by total disrupted motif hits
tf_order <- unique_hits %>%
  filter(motifs_disrupted == TRUE) %>%
  group_by(family) %>%
  summarise(total_hits = sum(!is.na(match_confidence))) %>%
  arrange(desc(total_hits)) %>%
  pull(family)

# Step 2: Plot with fill by match_confidence
p4 <- unique_hits %>%
  filter(motifs_disrupted == TRUE) %>%
  group_by(family) %>%
  summarise(n_hits = n(), .groups = "drop") %>%
  mutate(TF_family = factor(family, levels = rev(tf_order))) %>%
  ggplot(aes(y = TF_family, x = n_hits, fill = TF_family)) +
  geom_col(color="black") +
  scale_fill_viridis_d(option = "D", guide = "none") +
  labs(
    x = "Count" ,
    title = "Variant Overlaps")+
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.y = element_blank(),
        axis.title.y= element_blank())
# ) +
# scale_x_continuous(
#   breaks = seq(0, 100000, by = 20000),
#   labels = label_number(scale_cut = cut_short_scale())
# )

p4

#read variant predicted effect
var_prediction<-read_tsv("~/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_scores/average_cd4_tcells_AJ_common_variants.mean.variant_scores.tsv") %>%
  filter(variant_id %in% vars$variant_id) %>% dplyr::select(c(variant_id, logfc.mean)) %>% unique()


##boxplot
plot_df <- var_prediction %>%
  inner_join(variant_motif_map, by = "variant_id") %>%
  filter(!is.na(logfc.mean))  %>% mutate(TF_family = factor(TF_family, levels = rev(tf_order))) 
n_colors <- length(levels(plot_df$TF_family))
my_colors <- colorRampPalette(brewer.pal(8, "YlGnBu"))(n_colors)

fc_boxplot <- ggplot(plot_df, aes(x = logfc.mean, y = TF_family, fill = TF_family)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "gray30") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", linewidth = 0.4) +  # ← added line
  scale_fill_manual(values = my_colors) +
  labs(
    x = "Allelic logFC (Alt/Ref)"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 12),
    legend.position = "none"
  )

fc_boxplot

importance_boxplot <- ggplot(plot_df, aes(x = hit_coefficient_global, y = TF_family, fill = TF_family)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "gray30") +
  scale_fill_manual(values = my_colors) +
  labs(
    x = "Global Hit Importance"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 12),
    legend.position = "none"
  )

importance_boxplot

#counts

pt.varcounts<-variant_motif_map %>% group_by(TF_family)  %>% summarise(n_vars=n()) %>%
  mutate(TF_family = factor(TF_family, levels = rev(tf_order))) %>%
  ggplot(aes(y=TF_family, x=n_vars, fill = TF_family)) +
  geom_col(color = "black", linewidth=0.4)+
  geom_text(aes(label=n_vars, hjust = -0.3)) +
  scale_fill_manual(values = my_colors) +
  labs(
    x = "# Variants"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 12),
    legend.position = "none"
  )
pt.varcounts
#Scatter plot, n disrupted motifs vs total motifs

#Get all hits per TF
df_all <- unique_hits %>%
  group_by(family) %>%
  summarise(all_hits = n(), .groups = "drop")

# Get disrupted hits per TF
df_disrupted <- unique_hits %>%
  filter(motifs_disrupted == TRUE) %>%
  group_by(family) %>%
  summarise(disrupted_hits = n(), 
            hit_importance.mean=mean(hit_coefficient_global),
            .groups = "drop")

#Merge the two
df_scatter <- left_join(df_all, df_disrupted, by = "family") %>%
  mutate(disrupted_hits = replace_na(disrupted_hits, 0))  %>% mutate(TF_family = factor(family, levels = rev(tf_order)))



var_prediction_summary <- var_prediction %>%
  inner_join(variant_motif_map, by = "variant_id") %>%              # <- was join=
  group_by(TF_family) %>%
  summarise(
    median_allelic_fc = median(logfc.mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(df_scatter, by = "TF_family")        

# Scatter plot
p_scatter <- ggplot(df_scatter, aes(x = all_hits, y = disrupted_hits, label = TF_family, fill = TF_family)) +
  geom_point(shape = 21, color = "black", size = 6) +  # shape 21 allows fill
  geom_text(check_overlap = TRUE, size = 3, vjust = -1.2, hjust=0.5) +
  scale_fill_manual(values = my_colors) +
  scale_x_continuous(
    breaks = seq(0, 100000, by = 20000),
    labels = label_number(scale_cut = cut_short_scale())
  ) +
  labs(
    x = "Total motif hits",
    y = "Motif Hits Overlapping Variants"  ) +
  theme_classic() +
  theme(legend.position = "none")

p_scatter

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
  "PUTATIVE_RUNX2"          = "RUNX",
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

label_plot <- df_disrupted %>%
  filter(family %in% tf_order) %>%
  mutate(TF_family = factor(family, levels = rev(tf_order))) %>%
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



summary_df<-read_delim("/gpfs/commons/home/mmatos/cd4_qtl_paper_figures/figure_1/data/eqtl_caqtl_finemapping_coloc_all.tsv")

### what fraction of chromBPNet variants are any type of QTLs?
##### since a single ChromBPNet variant can show up in multiple credible sets and multiple QTL types we are going global QTL appearance
#--- summary_df contains the all finemapped variants for both eQTLs and caQTLs
# ChromBPNet variants universe
chrombpnet_vars <- vars$variant_id 
## ---- ca side: all variants that are in any caQTL credible set ----
ca_all <- summary_df %>%
  mutate(
    var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])",
                  "chr\\1_\\2_\\3_\\4",
                  variant_id.x)
  ) %>%
  semi_join(vars, by = c("var_id" = "variant_id")) %>%
  pull(var_id) %>%
  unique()

## ---- eq side: all variants that are in any eQTL credible set ----
eq_all <- summary_df %>%
  mutate(
    var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])",
                  "chr\\1_\\2_\\3_\\4",
                  variant_id.y)
  ) %>%
  semi_join(vars, by = c("var_id" = "variant_id")) %>%
  pull(var_id) %>%
  unique()

rm(summary_df)

qtl_any <- union(ca_all, eq_all)

# find 

variant_status <- tibble(
  variant_id = chrombpnet_vars
) %>%
  mutate(
    is_qtl = variant_id %in% qtl_any,
    category = ifelse(is_qtl, "QTL_shared", "ChromBPNet_specific")
  )


var_prediction_summary <- var_prediction %>%
  left_join(variant_status, by = "variant_id") %>%
  inner_join(variant_motif_map, by = "variant_id") %>%
  group_by(TF_family) %>%
  summarise(
    n_variants        = dplyr::n(),
    n_qtl             = sum(is_qtl, na.rm = TRUE),   # or is_qtl
    frac_qtl          = n_qtl / n_variants,
    median_allelic_fc = median(logfc.mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(df_scatter, by = "TF_family") 


sum_df_long <- var_prediction_summary %>%
  transmute(
    TF_family,
    n_qtl,
    n_nonqtl = n_variants - n_qtl
  ) %>%
  tidyr::pivot_longer(
    cols = c(n_qtl, n_nonqtl),
    names_to = "var_type",
    values_to = "n"
  ) %>%
  mutate(
    bucket = recode(var_type,
                    n_qtl    = "QTL_shared",
                    n_nonqtl = "ChromBPNet_spec"),
    TF_family = factor(TF_family, levels = rev(tf_order))
  )

motif_qtl.plt <- ggplot(sum_df_long, aes(x = TF_family, y = n, fill = bucket)) +
  geom_col(position = "fill", width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = variant_cmap3, drop = FALSE) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "QTL overlap", fill = "Variant Set")


write_tsv(variant_status, "~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/cd4_top_cpbnet_variants_motif_overlap_moreTFS_qtl_overlap_boolean_v2.tsv")

#########--- check for motif overlap on peak annotations (rpomoter vs. distal)

#read atac called peaks for chromBPnet  
bp.peaks<-fread("~/cd4_chrombpnet/data/inputs/peaks/merged_cd4_samples_peaks_no_blacklist.sorted.bed")
bp.peaksGR <- GRanges(bp.peaks$V1, IRanges(bp.peaks$V2, bp.peaks$V3), peak_name = bp.peaks$V4)

#read the gencode gene annotation
gencode <- fread("~/resources/genome/hg38_gencode_raw/gencode.v44.primary_assembly.annotation.genes.bed",
                 sep="\t", data.table=FALSE) %>%
  #filter(V8 == "protein_coding",
  #      grepl("^chr(\\d+|X|Y)$", V1)) %>%                 # drop alts/patches
  transmute(
    chr        = V1,
    strand     = V6,
    gene_id    = V7,
    gene       = V4,                                       # gene symbol
    tss_1based = if_else(V6 == "+", V2 + 1L, V3)           # BED start is 0-based
  ) %>%
  distinct(chr, gene, .keep_all = TRUE)                    # ensure 1 row per chr×gene


###### data driven annotation of regulatory elements ####
### promoter: peak needs to overlap at least one TSS of a gene
### distal regulatory element -- no TSS overlap


# 1) Build 1-bp TSS GRanges
tss_gr <- with(gencode,
               GRanges(
                 seqnames = chr,
                 ranges   = IRanges(start = tss_1based, end = tss_1based),
                 strand   = strand
               )
)

# 5) Overlap: promoter if overlaps ANY TSS (ignore strand for promoter call)
#tss_win <- resize(tss_gr, width = 201, fix = "center")  # ±100bp
ov_tss  <- countOverlaps(bp.peaksGR, tss_gr, ignore.strand = TRUE)

# Write annotations back onto the GRanges
mcols(bp.peaksGR)$promoter   <- ov_tss > 0
mcols(bp.peaksGR)$annotation <- ifelse(mcols(bp.peaksGR)$promoter, "promoter", "distal")

## now annotate motif_vars
vars_GR<-readRDS("~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap/cd4_top_cpbnet_variants_granges_motif_overlap_moreTFS_boolean.rds")

# 1) Overlap mapping (variant -> peak)
hits <- findOverlaps(vars_GR, bp.peaksGR, ignore.strand = TRUE)

## 3) Transfer peak metadata to variants
# Prepare vectors aligned to variants
ann_vec  <- rep(NA_character_, length(vars_GR))
name_vec <- rep(NA_character_, length(vars_GR))

ann_vec[queryHits(hits)]  <- mcols(bp.peaksGR)$annotation[subjectHits(hits)]
name_vec[queryHits(hits)] <- mcols(bp.peaksGR)$peak_name[subjectHits(hits)]

mcols(vars_GR)$peak_annotation <- ann_vec
mcols(vars_GR)$peak_name       <- name_vec


var_tss <- tibble(
  var_id= vars_GR$variant_id,
  peak_annotation=vars_GR$peak_annotation
)

##join with back to the chrombpnet var df contaitning the motifs ids
var_prediction_summary2 <- var_prediction %>%
  left_join(variant_status, by = "variant_id") %>%
  left_join(var_tss, join_by(variant_id==var_id)) %>%
  inner_join(variant_motif_map, by = "variant_id")


sum_df2 <- var_prediction_summary2 %>%
  filter(!is.na(TF_family), !is.na(peak_annotation)) %>%
  group_by(TF_family, peak_annotation) %>%
  summarise(n_vars = n_distinct(variant_id), .groups = "drop") %>%
  mutate(
    TF_family = factor(TF_family, levels = rev(tf_order))
  )



# 3) Plot: 100% stacked by peak annotation, faceted by bucket
motif_tss.plt <- ggplot(sum_df2, aes(x = TF_family, y = n_vars, fill = peak_annotation)) +
  geom_col(position = "fill", width = 0.75) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Proportion of motif-overlapping variants",
    fill = "Genomic annotation"
  )

motif_tss.plt
########



grid<-plot_grid(label_plot, p,fc_boxplot,importance_boxplot,pt.varcounts, motif_qtl.plt, motif_tss.plt, rel_widths = c(1.1, 1,1,1, 1.2, 1.2, 1.2), nrow = 1,
                align = "h",
                axis = "tb",
                label_size = 10,
                labels =c("TF Family", "De Novo Pattern", "CA Allelic Effect", 'Hit Importance', "Ovelapping Variants", "QTL overlap", "Annotation"))
grid

ggsave("~/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/variant_motif_overlap/plots/variant_motif_overlap_summary_tf_family_coeff5_b1_v3.pdf", grid, width =12, height = 5)




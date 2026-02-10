# ============================================================
# Setup
# ============================================================
.libPaths(c(
  "/gpfs/commons/home/mmatos/R/x86_64-pc-linux-gnu-library/4.4",
  "/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"
))

setwd("~/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/cbpnet_examples")

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(data.table)
  library(tidyr)
  library(readr)
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
  library(BSgenome)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(gridExtra)
  library(readr)
})

source("~/cd4_qtl_paper_figures/utils/track_plots_helpers.R")
source("~/cd4_qtl_paper_figures/utils/cbpnet_extract_allelic_contrib_plus_atrr_bigwigs_v2.R")

# Ensure exporter initialized once
stopifnot(exists(".init_variant_bw_exporter"))
.init_variant_bw_exporter()

# ============================================================
# Paths / Inputs
# ============================================================
paths <- list(
  atrr        = "/gpfs/commons/home/mmatos/cd4_chrombpnet/chrombpnet_model_b7/variant_contribution_scores/averaged_cd4_tcells_AJ_common_variants.shap.counts.h5",
  data_h5     = "/gpfs/commons/home/mmatos/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_variant_prediction_scores.h5",
  tsv         = "/gpfs/commons/home/mmatos/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_scores/average_cd4_tcells_AJ_common_variants.mean.variant_scores.tsv",
  chrom_sizes = "/gpfs/commons/home/mmatos/resources/genome/hg38.chrom.sizes",
  bw_path     = "/gpfs/commons/groups/lappalainen_lab/mmatos/cd4_aging_project/data/ATAC-seq/atac_preprocessing/nextflow_tmp/a4/4c07c68cd77cef99a8b320e8e66238/cd4_atac.bigWig",
  peaks_bed   = "/gpfs/commons/groups/lappalainen_lab/mmatos/cd4_aging_project/data/ATAC-seq/atac_preprocessing/merged_library/peak_calling/MACS3/BAMPE/peaks_102024/cd4_atac_padded_summits_peaks.bed",
  filtered_peaks="/gpfs/commons/home/mmatos/cd4_caQTL_analysis/variant_to_peak_QTL/run_012625_qc_aware_qsmooth_CPM_MAF5_FDR5_1MB/results/001_peaks/peak_names.txt",
  dynamic_int = "~/cd4_CellRegMap/002_interaction_analysis/results/results_01312025/res_34_5factors/2025_03_21_mofa_interaction_res34_significant_summary_fdr1.csv",
  var_pred    = "~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap_hit_caller/unique_variants_overlapping_motifs.tsv",
  motifs_tsv  = "~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap_hit_caller/unique_hits_filtered_motif_instances.tsv",
  summary_df  = "/gpfs/commons/home/mmatos/cd4_qtl_paper_figures/figure_1/data/eqtl_caqtl_finemapping_coloc_all.tsv",
  coloc_csv_tpl = "/gpfs/commons/groups/lappalainen_lab/sghatan/marlis_pj/coloc/coloc_results/ca_eqtl_coloc/%s_coloc_results.csv",
  metadata = "/gpfs/commons/home/mmatos/ATAC-seq_analysis/diff_accesibility_ana/results/metadata/cd4_atac_metadata_scrna_union.csv",
  dir = '/gpfs/commons/home/mmatos/cd4_caQTL_analysis/variant_to_peak_QTL/run_012625_qc_aware_qsmooth_CPM_MAF5_FDR5_1MB/results/004_genotypes/plink/per_chr_plink_files/chr',
  bigwig_dir = "/gpfs/commons/groups/lappalainen_lab/mmatos/cd4_aging_project/data/ATAC-seq/atac_preprocessing/STAR/bigwig"
)
# ============================================================
# Helpers
# ============================================================

parse_variant_id <- function(variant_id) {
  x <- trimws(variant_id)
  
  mA <- stringr::str_match(x, "^chr([0-9XYMT]+)_([0-9]+)_[ACGT]_[ACGT]$")
  if (!is.na(mA[1, 1])) return(list(chr = paste0("chr", mA[1, 2]), pos = as.integer(mA[1, 3])))
  
  mB <- stringr::str_match(x, "^([0-9XYMT]+):([0-9]+)\\[b\\d+\\][ACGT],[ACGT]$")
  if (!is.na(mB[1, 1])) return(list(chr = paste0("chr", mB[1, 2]), pos = as.integer(mB[1, 3])))
  
  stop('Unrecognized variant_id format. Allowed: "chr6_113306741_G_A" or "6:113306741[b38]G,A".')
}

make_ld_bins <- function(r2) {
  cut(
    r2,
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    labels = c("0.0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"),
    include.lowest = TRUE,
    right = FALSE
  )
}

ld_colors <- c(
  "0.0–0.2" = "gray",
  "0.2–0.4" = "skyblue",
  "0.4–0.6" = "green",
  "0.6–0.8" = "orange",
  "0.8–1.0" = "red"
)

add_ld_from_lead <- function(markers_df, LD_mat, lead_pos) {
  LD_mat <- as.matrix(LD_mat)
  lead_index <- which.min(abs(markers_df$pos_hg38 - lead_pos))
  ld_cor <- LD_mat[, lead_index]
  markers_df$ld_r_gwas  <- ld_cor
  markers_df$ld_r2_gwas <- ld_cor^2
  markers_df$ld_bin_gwas <- make_ld_bins(markers_df$ld_r2_gwas)
  markers_df
}

plot_qtl_scatter <- function(markers_df, minStart, maxEnd, label_txt, region,
                             ylab_expr = bquote(-log[10](p))) {
  
  markers_df <- markers_df %>%
    dplyr::mutate(
      p = 2 * stats::pnorm(-abs(z_pos0)),   # two-sided p from Z
      y = -log10(p)
    )
  
  y_finite <- markers_df$y[is.finite(markers_df$y)]
  ymax <- if (length(y_finite)) max(y_finite, na.rm = TRUE) else 1
  
  step <- dplyr::case_when(
    ymax <= 5  ~ 1,
    ymax <= 20 ~ 2,
    ymax <= 50 ~ 5,
    TRUE       ~ 10
  )
  maxlimit <- ceiling(ymax / step) * step
  
  ggplot(markers_df, aes(x = pos_hg38, y = y)) +
    geom_point(
      aes(fill = ld_bin_gwas),
      shape = 21, size = 2, alpha = 0.8, color = "gray20", stroke = 0.2
    ) +
    scale_fill_manual(values = ld_colors, name = expression(LD~(r^2))) +
    scale_x_continuous(expand = c(0, 0), limits = c(minStart, maxEnd)) +
    coord_cartesian(ylim = c(0, maxlimit)) +
    annotate(
      "text",
      x = minStart,
      y = maxlimit,
      label = label_txt,
      hjust = 0, vjust = 1,
      size = 4
    ) +
    labs(y = ylab_expr, x = "") +
    BPCells:::trackplot_theme() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 11, angle = 90, vjust = 0.5),
      legend.position = "none",
      plot.title = element_text(size = 8, hjust = 0.5)
    ) +
    highlight_region(region, color = "skyblue", alpha = 0.2)
}



# ============================================================
# Load reference data
# ============================================================
transcripts <- read_gencode_genes("./references/", release = "44")

dynamic_int <- read.csv(paths$dynamic_int, row.names = 1)

# ============================================================
# Selecting the example to plot
# ============================================================
# var_prediction <- fread(paths$var_pred)
# 
# motifs_tbl <- fread(paths$motifs_tsv) %>%
#    select(start, end, variant_hit, family, strand, chr, variant_loc) %>%
#    mutate(
#      variant_loc2      = paste0(chr, "_", variant_loc),
#      motifs_disrupted  = variant_hit,
#      TF_family         = family
#    )
# 
# var_prediction <- var_prediction %>%
#    left_join(motifs_tbl, join_by(var_chr_pos == variant_loc2))
# 
# summary_df <- read_delim(paths$summary_df)
# 
# chrombpnet_vars <- var_prediction$variant_id
# 
# #get caqtls that contain chrombpnet variants
# ca_all <- summary_df %>%
#    mutate(var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])", "chr\\1_\\2_\\3_\\4", variant_id.x)) %>%
#    filter(var_id %in% chrombpnet_vars)
# #get eqtl that contain chrombpnet variants
# eq_all <- summary_df %>%
#    mutate(var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])", "chr\\1_\\2_\\3_\\4", variant_id.y)) %>%
#    filter(var_id %in% chrombpnet_vars)
# 
# qtl_any <- bind_rows(ca_all, eq_all) %>% filter(gene %in% dynamic_int$gene)
# rm(summary_df, ca_all, eq_all)
# 
# qtl_any_dyn <- qtl_any %>%
#  left_join(var_prediction, join_by(var_id == variant_id)) %>%
#  select(
#    var_id, finemapped_cs_coloc, chr.x.x, gene, peak, variant_id.x, variant_id.y, PP.H4.abf, coloc_status,
#    caQTL_variant, eQTL_variant,  pos, allele1, allele2, logfc.mean, jsd.mean, var_chr_pos,
#    is_qtl, category, peak_annotation, variant_loc,  motifs_disrupted ,TF_family)
# 
# fwrite(qtl_any_dyn, "~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/molqtls_chrombpnet_dynamic_list.csv" )
# rm(qtl_any)
qtl_any_dyn<-fread("~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/molqtls_chrombpnet_dynamic_list.csv") 

test<-qtl_any_dyn %>% filter(gene=='GIMAP4')
unique(test$var_id)
 # ============================================================
 # END
 # ============================================================
# Peaks BED -> GRanges
peak_names <- read_tsv(paths$filtered_peaks, col_names = FALSE)
bed_peaks <- read_tsv(paths$peaks_bed, col_select = c(1:4), col_names = FALSE) %>% filter(X4 %in% peak_names$X1)
setnames(bed_peaks, c("peak_chr", "peak_start", "peak_end", "peak_name"))

bed_peaks_gr <- GRanges(
  seqnames = bed_peaks$peak_chr,
  ranges   = IRanges(start = bed_peaks$peak_start, end = bed_peaks$peak_end),
  name     = gsub("merged_cd4_samples_", "", bed_peaks$peak_name)
)

# Motifs GRanges for genome annotation track
motifs <- fread(paths$motifs_tsv) %>%
  select(start, end, variant_hit, family, strand, chr, variant_loc) %>%
  mutate(
    variant_loc2      = paste0(chr, "_", variant_loc),
    motifs_disrupted  = variant_hit,
    TF_family         = family
  )
motifs$motifs_disrupted <- factor(motifs$motifs_disrupted, levels = c(TRUE, FALSE))

# ============================================================
# Example settings
# ============================================================
gene    <- "GIMAP4"
variant <- "chr7_150567557_C_T"
vp <- parse_variant_id(variant)
var_contig <- vp$chr
var_pos    <- vp$pos

# Coloc results 
coloc_csv <- sprintf(paths$coloc_csv_tpl, var_contig) 
coloc_results <- fread(coloc_csv)


# Windows
region_window_broad  <- 400000L
region_window_broad2 <- 50000L
region_window_broad2_5 <-20000L
region_window_broad3 <- 500L
region_window        <- 50L

# Regions (strings) + GRanges
region_broad      <- paste0(var_contig, ":", max(1, var_pos - region_window_broad),  "-", var_pos + region_window_broad)
region_zoom2      <- paste0(var_contig, ":", max(1, var_pos - region_window_broad2), "-", var_pos + region_window_broad2)
region_zoom2_5      <- paste0(var_contig, ":", max(1, var_pos - region_window_broad2_5), "-", var_pos + region_window_broad2_5)
region_zoom3      <- paste0(var_contig, ":", max(1, var_pos - region_window_broad3), "-", var_pos + region_window_broad3)
region_zoom       <- paste0(var_contig, ":", max(1, var_pos - region_window),        "-", var_pos + region_window)
region_variant    <- paste0(var_contig, ":", var_pos, "-", var_pos + 1)  # 1bp window (start <= end)

region_broad_gr <- str_to_gr(region_broad)
region_zoom2_gr <- str_to_gr(region_zoom2)
region_zoom2_5_gr <- str_to_gr(region_zoom2_5)


minStart <- max(1, var_pos - region_window_broad)
maxEnd   <- var_pos + region_window_broad

# ============================================================
# caQTL finemapping scatter
# ============================================================
x <- 1158
region_ca <- gsub("^7:", "", coloc_results$region.x[x])
peak_ca   <- coloc_results$peak[x]
caqtl_variant <- coloc_results$caQTL_variant[x]

prep_ca   <- prepare_caqtl_data(var_contig, peak_ca, region_ca, caqtl_variant)
markers_ca <- prep_ca$markers
LD_ca      <- prep_ca$LD

markers_ca <- markers_ca %>%
  mutate(
    pos_hg38 = as.integer(str_extract(marker, "(?<=:)\\d+")),
    z_pos0   =z_caqtl
  ) %>%
  add_ld_from_lead(LD_ca, lead_pos = var_pos)

p.caqtl <- plot_qtl_scatter(
  markers_df = markers_ca,
  minStart   = minStart,
  maxEnd     = maxEnd,
  region=region_zoom2,
  label_txt  = paste0("caQTL ", caqtl_variant, " ", peak_ca)
)

# ============================================================
# eQTL finemapping scatter
# ============================================================
region_eq <- gsub("^7:", "", coloc_results$region.x[x])
gene_eq   <- coloc_results$gene[x]
eqtl_variant <- coloc_results$eQTL_variant[x]

prep_eq   <- prepare_eqtl_data(var_contig, gene_eq, region_eq, variant)
markers_eqtl <- prep_eq$markers
LD_eqtl      <- prep_eq$LD

markers_eqtl <- markers_eqtl %>%
  mutate(
    pos_hg38 = as.integer(str_extract(marker, "(?<=:)\\d+")),
    z_pos0   = z_eqtl
  ) %>%
  add_ld_from_lead(LD_eqtl, lead_pos = var_pos)

p.eqtl <- plot_qtl_scatter(
  markers_df = markers_eqtl,
  minStart   = minStart,
  maxEnd     = maxEnd,
  region=region_zoom2,
  label_txt  = paste0("eQTL ", eqtl_variant, " ", gene_eq)
)

# ============================================================
# Observed ATAC coverage + peaks + genes
# ============================================================
atac_sub <- import.bw(paths$bw_path, which = region_broad_gr)

atac_cov <- trackplot_bw(
  bw          = atac_sub,
  region      = region_broad,
  facet_label = "ATAC",
  plot_as     = "area",
  tile_width  = 1,
  track_label = "ATAC",
  color       = "darkgreen"
) +
  highlight_region(region_zoom, color = "red", alpha = 0.2) +
  highlight_region(region_zoom2_5, color = "yellow", alpha = 0.2) +
  theme(legend.position = "none", plot.title = element_text(size = 8, hjust = 0.5)) 
  

p.peak_delim <- trackplot_peak_delim(bed_peaks_gr, region_broad, "darkgreen") +
  highlight_region(region_zoom, color = "red", alpha = 0.2) +
  highlight_region(region_zoom2_5, color = "yellow", alpha = 0.2) +
  theme(legend.position = "none", plot.title = element_text(size = 8, hjust = 0.5))

track_genes <- trackplot_gene_custom(transcripts, str_to_gr(region_broad)) +
  guides(color = "none") +
  highlight_region(region_zoom, color = "red", alpha = 0.2) +
  highlight_region(region_zoom2_5, color = "yellow", alpha = 0.2) +
  theme(legend.position = "none", plot.title = element_text(size = 8, hjust = 0.5))

# ============================================================
# Observed ATAC coverage + peaks + genes (ZOOMED)
# ============================================================
metadata<-read_csv(paths$metadata)

track_data <- prepareTrackData_grouped_genotype(
  variant_id = caqtl_variant,
  metadata   = metadata,
  base_dir   = paths$dir,
  bigwig_dir = paths$bigwig_dir,
  flanking = c(20000, 20000)
)

p.peaks <- plotGroupedCoverage_allelic_fill(
  track_data, region_zoom2_5_gr,
  genotype_levels = c("Homozygous Alt","Heterozygous", "Homozygous Ref" ),
  rasterize = TRUE, raster_dpi = 400, alpha = 0.5) 


p.peak_delim2_5 <- trackplot_peak_delim(bed_peaks_gr, region_zoom2_5, "darkgreen") 

track_genes2_5 <- trackplot_gene_custom(transcripts, str_to_gr(region_zoom2_5)) +
  guides(color = "none")
# ============================================================
# Export / import ChromBPNet bigWigs for the variant
# ============================================================
dir.create("bw_out", showWarnings = FALSE, recursive = TRUE)

res <- export_variant_bigwigs_safe(
  pred_h5          = paths$data_h5,
  attr_h5          = paths$atrr,
  tsv_path         = paths$tsv,
  variant_id       = variant,
  chrom_sizes_path = paths$chrom_sizes,
  outdir           = "bw_out/",
  write_channels   = FALSE
)

ref_bw      <- import.bw(res$pred_REF)
alt_bw      <- import.bw(res$pred_ALT)
atrr_ref_bw <- import.bw(res$contribSum_REF)
atrr_alt_bw <- import.bw(res$contribSum_ALT)

# ============================================================
# Pred + contribution tracks
# ============================================================
pred <- trackplot_bw_allelic(
  bw_ref      = ref_bw,
  bw_alt      = alt_bw,
  region      = region_zoom3,
  track_label = "ATAC signal",
  facet_label = "ATAC signal",
  score_cmap  = c(ref = "darkgreen", alt = "darkviolet"),
  ymin_zero   = TRUE,
  score_shift = 1
) +
  highlight_region(region_zoom, color = "skyblue", alpha = 0.2) +
  theme(legend.position = "none", plot.title = element_text(size = 8, hjust = 0.5)) 

all_scores <- c(as.numeric(mcols(atrr_ref_bw)$score), as.numeric(mcols(atrr_alt_bw)$score))
all_scores <- all_scores[is.finite(all_scores)]

ymax <- max(all_scores) * 1.05
ymin <- min(all_scores) * 1.05
ylim_shared <- c(ymin, ymax)

ymax_accuracy <- 10^as.integer(log10(0.01 * abs(ymax)))
ymin_accuracy <- 10^as.integer(log10(0.01 * abs(ymin)))

range_label <- sprintf(
  "[%s-%s]",
  scales::label_comma(accuracy = ymin_accuracy, big.mark = " ")(ymin),
  scales::label_comma(accuracy = ymax_accuracy, big.mark = " ")(ymax)
)

Ref_contribs <- trackplot_contribs(
  bw          = atrr_ref_bw,
  region      = region_zoom,
  genome      = BSgenome.Hsapiens.UCSC.hg38,
  ylim        = ylim_shared,
  range_label = range_label,
  track_label = "REF"
) +
  highlight_relative_region(region_window + 0.5, region_window + 1.5, color = "gold", alpha = 0.25)

Alt_contribs <- trackplot_contribs(
  bw          = atrr_alt_bw,
  region      = region_zoom,
  genome      = BSgenome.Hsapiens.UCSC.hg38,
  ylim        = ylim_shared,
  range_label = range_label,
  track_label = "ALT"
) +
  highlight_relative_region(region_window + 0.5, region_window + 1.5, color = "gold", alpha = 0.25)

track_hits <- trackplot_genome_annotation(
  loci        = motifs,
  region      = region_zoom,
  color_by    = "motifs_disrupted",
  show_strand = TRUE,
  colors      = c("red", "gray"),
  label_size  = 3,
  label_by    = "TF_family",
  track_label = "TF Motifs"
)


# ============================================================
# Combine trackplot
# ============================================================
title_txt <- paste0(variant, ": ", str_to_pretty(region_zoom))

track_plot <- trackplot_combine2(
  list(
    BPCells:::wrap_trackplot(p.caqtl,       unit(0.20, "null")) + indiv_theme,
    BPCells:::wrap_trackplot(p.eqtl,        unit(0.20, "null")) + indiv_theme,
    BPCells:::wrap_trackplot(atac_cov,      unit(0.18, "null")) + indiv_theme,
    BPCells:::wrap_trackplot(p.peak_delim,  unit(0.10, "null")) + indiv_theme,
    BPCells:::wrap_trackplot(track_genes,   unit(0.10, "null")) + indiv_theme
    
  ),
  title = title_txt
)

track_plot2 <- trackplot_combine2(
  list(
    BPCells:::wrap_trackplot(pred,          unit(0.22, "null")) + indiv_theme,
    BPCells:::wrap_trackplot(Ref_contribs,  unit(0.17, "null")) + indiv_theme,
    BPCells:::wrap_trackplot(Alt_contribs,  unit(0.17, "null")) + indiv_theme,
    BPCells:::wrap_trackplot(track_hits,    unit(0.10, "null")) + indiv_theme
    
  ),
  title = title_txt
)

track_plot3 <- trackplot_combine2(
  list(
    BPCells:::wrap_trackplot(p.peaks,      unit(0.22, "null")) + indiv_theme,
    BPCells:::wrap_trackplot(p.peak_delim2_5,  unit(0.10, "null")) + indiv_theme,
    BPCells:::wrap_trackplot(track_genes2_5,   unit(0.10, "null")) + indiv_theme
    
  ),
  title = title_txt
)

pt1<-cowplot::plot_grid(track_plot2,track_plot3, align = "h")
pt<-cowplot::plot_grid(track_plot, pt1, nrow = 2)
# ============================================================
# Save
# ============================================================
ggsave(
  filename = "~/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/cbpnet_examples/plots/GIMAP4_coloc_qtl_chrom_dyna2.pdf",
  plot     = pt,
  width    = 10, height = 6, units = "in"
)
pt1<-cowplot::plot_grid(track_plot2,track_plot3, align = "v", ncol=1)
ggsave(
  filename = "~/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/cbpnet_examples/plots/GIMAP4_coloc_qtl_chrom_dyna_zoomplots.pdf",
  plot     = pt1,
  width    = 7, height = 6, units = "in"
)
ggsave(
  filename = "~/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/cbpnet_examples/plots/GIMAP4_coloc_qtl_chrom_dyna_coloc.pdf",
  plot     = track_plot,
  width    = 6, height = 6, units = "in"
)

# ============================================================
# Dynamic eQTL
# ============================================================

library(dplyr)
library(LaCroixColoR)
library(ggplot2)
library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)
library(reticulate)
library(SeuratDisk)
library(cowplot)
options(bitmapType = "cairo")
source("~/cd4_qtl_paper_figures/utils/dynamic_mofa_plot_beta_helper.R")
source("~/cd4_qtl_paper_figures/utils/qtl_violion_helper.R")



gene_name="GIMAP4"
index <- which(dynamic_int$gene == gene_name)
gene <- dynamic_int$gene[index]
chr  <- dynamic_int$chr[index]
variant_id <- dynamic_int$snpID[index]
pos <- as.character(sub(".*:(\\d+)\\[b38\\].*", "\\1", variant_id))


pseudocells=("/gpfs/commons/home/mmatos/cd4_CellRegMap/001_preprocessing/results_01312025/prep_phenotype_vector/Data/pseudobulk_metacells/Pseudobulk_per_donor_all_cells_all_conditions_Leiden_res34.h5")
mofa = "/gpfs/commons/home/mmatos/cd4_CellRegMap/001_preprocessing/results_01312025/mofa_factors/mofa_trained/cd4_aging_filt_sce_mofa_expectations_model_factors_res34.csv"

sce = readH5AD(pseudocells)
df = read.csv(mofa, row.names = 1)
head(df,2)

rownames(df) = colnames(sce)
head(df,2)

#Add mofa factors as the PCA
reducedDim(sce, "MOFA") <- df 
head(reducedDim(sce))

sce.seurat <- as.Seurat(sce, counts = "X", data = "X")

pt1<-plot_gene_mofa(gene, chr, pos, 
                    x_factor = "MOFA1", 
                    y_factor = "MOFA4",
                    nclus=10,
                    alpha = 0.9,
                    output_dir = "~/cd4_qtl_paper_figures/figure_1/plotting/plots_oct2025/dynamic/")
cowplot::plot_grid(pt1)


pt3<-FeaturePlot(sce.seurat, reduction = "MOFA", dims = c(1, 4), features = c("CCR7"), pt.size = 1)
pt4<-FeaturePlot(sce.seurat, reduction = "MOFA", dims = c(1, 4), features = c("FUT7"), pt.size = 1)

pt5<-plot_eqtl_violin_from_string(
  gene           = gene_name,
  variant_string = variant_id,
  celltype       = "allcells"
)


state<-cowplot::plot_grid(pt3, pt4, ncol=1)
 
g <- cowplot::plot_grid(state, pt1, pt5, rel_widths = c(1,1.5, 1), ncol=3)

ggsave("~/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/cbpnet_examples/plots/GIMAP4_dynamic_qtl.pdf", g, width = 12, height = 5, units = "in")

.libPaths(c("/gpfs/commons/home/mmatos/R/x86_64-pc-linux-gnu-library/4.4","/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"))
setwd("~/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/cbpnet_examples")
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(rtracklayer)
  library(BSgenome)
  library(data.table)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(gridExtra)   
})
source("~/cd4_qtl_paper_figures/utils/track_plots_helpers.R")
source("~/cd4_qtl_paper_figures/utils/cbpnet_extract_allelic_contrib_plus_atrr_bigwigs_v2.R")

#Ensure exporter initialized once
if (!exists(".init_variant_bw_exporter")) {
  stop(".init_variant_bw_exporter() not found. Did you source the Python helper?")
}
.init_variant_bw_exporter()

# ----------------------- Variant PARSER ----------------------------------------
parse_variant_id <- function(variant_id) {
  x <- trimws(variant_id)
  mA <- stringr::str_match(x, "^chr([0-9XYMT]+)_([0-9]+)_[ACGT]_[ACGT]$")
  if (!is.na(mA[1,1])) return(list(chr = paste0("chr", mA[1,2]), pos = as.integer(mA[1,3])))
  mB <- stringr::str_match(x, "^([0-9XYMT]+):([0-9]+)\\[b\\d+\\][ACGT],[ACGT]$")
  if (!is.na(mB[1,1])) return(list(chr = paste0("chr", mB[1,2]), pos = as.integer(mB[1,3])))
  stop('Unrecognized variant_id format. Allowed only: "chr6_113306741_G_A" or "6:113306741[b38]G,A".')
}

#get variants to plot
vars_motifs<-fread("~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap_hit_caller/unique_variants_overlapping_motifs.tsv")
head(vars_motifs)

vars_motifs_top16 <- vars_motifs %>%
  filter(abs_IPS.mean>0.1 & abs_logfc.mean >1.0) %>% 
  group_by(TF_annotation) %>%
  slice_max(hit_coefficient_global, n = 16, with_ties = FALSE) %>%
  ungroup() %>% arrange(desc(abs_logfc.mean))


vars_motifs_top16 <- vars_motifs %>%
  filter(jsd.mean > 0.1 & abs_logfc.mean >0.05) %>% 
  group_by(TF_annotation) %>%
  slice_max(hit_coefficient_global, n = 16, with_ties = FALSE) %>%
  ungroup() %>% arrange(desc(abs_logfc.mean))
# writeLines(vars_motifs_top16$variant_id, "~/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/cbpnet_examples/variant_plot.txt")
# # write one TSV per TF_annotation
# vars_motifs_top16 %>%
#   group_by(TF_annotation) %>%
#   group_walk(~ {
#     tf <- .y$TF_annotation[1]
#     tf_sane <- str_replace_all(tf, "[^A-Za-z0-9._-]", "_")
#     out_path <- paste0("top16_", tf_sane, ".tsv")
#     write_tsv(.x, out_path)
#   })

#variants <- as.character(vars_motifs_top16$variant_id)



######################## INPUT VCARIANT NAME OR AS ARGUMENT #######
variant="chr7_142701355_G_C"
# variant=args[1]
###########################

for (variant in vars_motifs_top16$variant_id) {
  variant="chr7_142701355_G_C"
  vp <- parse_variant_id(variant)
  var_contig <- vp$chr
  var_pos    <- vp$pos
  region_window_broad <- 100000L
  region_broad <- paste0(var_contig, ":", max(1, var_pos - region_window_broad), "-", var_pos + region_window_broad)
  region_broad.ls <- region_gr <- str_to_gr(region_broad)
  #region_ft <- list(chr = var_contig, start = max(1, var_pos - region_window_broad) , end = var_pos + region_window_broad)
  
  #1st zoom
  region_window_broad2 <-200L
  region_broad_zoom2 <- paste0(var_contig, ":", max(1, var_pos - region_window_broad2), "-", var_pos + region_window_broad2)
  
  #2st zoom
  region_window_broad2 <- 50L
  region_broad_zoom2 <- paste0(var_contig, ":", max(1, var_pos - region_window_broad2), "-", var_pos + region_window_broad2)
  
  region_window <- 50L
  region_zoom <- paste0(var_contig, ":", max(1, var_pos - region_window), "-", var_pos + region_window)
  
  region_variant<- paste0(var_contig, ":", max(1, var_pos+1), "-", var_pos)
  
  # ----------------------- RESOURCES -------------------------------------------
  # atrr: averaged variant contribution H5 across all folds
  # data_h5: averaged variant prediction H5 across all folds
  # tsv: output tsv from variant-scorer, or the input variant file, we need it to get the index of the variant and parse the predictions
  # chrom_sizes: hg38 chromsizes file
  # bw_path: BigWig from ATAC-seq data
  
  atrr         <- "/gpfs/commons/home/mmatos/cd4_chrombpnet/chrombpnet_model_b7/variant_contribution_scores/averaged_cd4_tcells_AJ_common_variants.shap.counts.h5"
  data_h5      <- "/gpfs/commons/home/mmatos/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_variant_prediction_scores.h5"
  tsv          <- "/gpfs/commons/home/mmatos/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_scores/average_cd4_tcells_AJ_common_variants.mean.variant_scores.tsv"
  chrom_sizes  <- "/gpfs/commons/home/mmatos/resources/genome/hg38.chrom.sizes"
  motifs_path  <- "~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap/motifs_hits_filtered_c90_40seq_w_variant_hits_granges_1b.rds"
  bw_path <- "/gpfs/commons/groups/lappalainen_lab/mmatos/cd4_aging_project/data/ATAC-seq/atac_preprocessing/nextflow_tmp/a4/4c07c68cd77cef99a8b320e8e66238/cd4_atac.bigWig"
  
  filtered_peaks<-'/gpfs/commons/home/mmatos/cd4_chrombpnet/data/inputs/peaks/merged_cd4_samples_peaks_no_blacklist.sorted.bed' #peakcall used for ChromBPNet
  bed_peaks <- read_tsv(filtered_peaks, col_select = c(1:4), col_names = FALSE) 
  colnames(bed_peaks) <- c("peak_chr", "peak_start", "peak_end", "peak_name")
  
  #motifs_with_tf_annotated_filt_c90_40seq_unique_variant_hits_granges.rds
  
  motifs       <- readRDS(motifs_path)
  
  #load observed ATAC coverage for the specified region
  atac_sub <- import.bw(bw_path, which = region_broad.ls)
  
  #----plot coverage ----------
  atac_cov <- trackplot_bw(
    bw          = atac_sub,
    region      = region_broad,
    facet_label = "ATAC",
    plot_as     = "area",
    tile_width  = 1,
    track_label = "ATAC",
    color       = "aquamarine4") +
    highlight_region(region_broad_zoom2, color = "skyblue", alpha = 0.2) +
    theme(legend.position = "none",
          plot.title = element_text(size = 8, hjust = 0.5)) 
  
  
  #make granges of narrowpeaks
  bed_peaks_gr <- GRanges(
    seqnames = bed_peaks$peak_chr,
    ranges = IRanges(start = bed_peaks$peak_start, end = bed_peaks$peak_end),
    name    = gsub("merged_cd4_samples_", "", bed_peaks$peak_name)
  )
  
  #plot peak delimiters ---------
  color="aquamarine4"
  p.peak_delim<-trackplot_peak_delim(bed_peaks_gr, region_broad, color) + #bed_peaks_gr: rgranges object of peaks 
    highlight_region(region_broad_zoom2, color = "skyblue", alpha = 0.2) +
    theme(legend.position = "none",
          plot.title = element_text(size = 8, hjust = 0.5)) 
  
  
  # ----------------------- EXPORT / IMPORT BIGWIGS ------------------------------
  dir.create("bw_out", showWarnings = FALSE, recursive = TRUE)
  
  # Export your variant with proper paths
  res <- export_variant_bigwigs_safe(
    pred_h5 = data_h5,
    attr_h5 = atrr,
    tsv_path = tsv,
    variant_id = variant ,
    chrom_sizes_path = chrom_sizes,  
    outdir = 'bw_out/', 
    write_channels = FALSE
  )
  
  ref_bw      <- rtracklayer::import.bw(res$pred_REF)
  alt_bw      <- rtracklayer::import.bw(res$pred_ALT)
  atrr_ref_bw <- rtracklayer::import.bw(res$contribSum_REF)
  atrr_alt_bw <- rtracklayer::import.bw(res$contribSum_ALT)
  
  # ----------------------- BUILD TRACKS -----------------------------------------
  pred <- trackplot_bw_allelic(
    bw_ref      = ref_bw, #BigWig ref allele (predictions)
    bw_alt      = alt_bw,#BigWig alt allele (predictions)
    region      = region_broad_zoom2,
    track_label = "ATAC signal",
    facet_label = "ATAC signal",
    score_cmap  = c(ref = "darkgreen", alt = "darkviolet"),
    ymin_zero   = TRUE,
    score_shift = 1
  ) +
    highlight_region(region_zoom, color = "skyblue", alpha = 0.2) +
    theme(legend.position = "none",
          plot.title = element_text(size = 8, hjust = 0.5)) 
  
  all_scores <- c(
    as.numeric(mcols(atrr_ref_bw)$score),
    as.numeric(mcols(atrr_alt_bw)$score)
  )
  all_scores <- all_scores[is.finite(all_scores)]
  
  ymax <- max(all_scores) * 1.05
  ymin <- min(all_scores) * 1.05
  
  ymax_accuracy <- 10^as.integer(log10(0.01 * abs(ymax)))
  ymin_accuracy <- 10^as.integer(log10(0.01 * abs(ymin)))
  
  range_label <- sprintf(
    "[%s-%s]",
    scales::label_comma(accuracy = ymin_accuracy, big.mark = " ")(ymin),
    scales::label_comma(accuracy = ymax_accuracy, big.mark = " ")(ymax)
  )
  
  ylim_shared <- c(ymin, ymax)
  
  Ref_contribs <- trackplot_contribs(
    bw          = atrr_ref_bw,#BigWig ref allele
    region      = region_zoom,
    genome      = BSgenome.Hsapiens.UCSC.hg38,
    ylim        = ylim_shared,
    range_label = range_label,
    track_label = "REF"
  ) +
    highlight_relative_region(region_window+0.5, region_window+1.5, color = "gold", alpha = 0.25)
  
  
  Alt_contribs <- trackplot_contribs(
    bw          = atrr_alt_bw, #BigWig alt allele
    region      = region_zoom,
    genome      = BSgenome.Hsapiens.UCSC.hg38,
    ylim        = ylim_shared,
    range_label = range_label,
    track_label = "ALT"
    
  ) +
    highlight_relative_region(region_window+0.5, region_window+1.5, color = "gold", alpha = 0.25)
  
  motifs$motifs_disrupted<-factor(motifs$motifs_disrupted, levels = c(TRUE, FALSE)) #setting factor levels so color is consistent across variants
  
  track_hits <- trackplot_genome_annotation(
    loci        = motifs, #GRanges
    region      = region_zoom,
    color_by    = "motifs_disrupted", #this a true/false column in my granges
    show_strand = TRUE,
    colors      = c("red", "gray"), 
    label_size  = 3,
    label_by    = "TF_family",
    track_label = "TF Motifs"
  )
  
  transcripts <- read_gencode_genes("./references/", release="44")
  head(transcripts)
  
  track_genes <- trackplot_gene_custom(transcripts, region_broad.ls) +
    ggplot2::guides(color = "none")
  track_genes
  
  # subset_g    <- get_subset_genes_for_region(region_broad)
  # track_genes <- trackplot_gene(subset_g, region_broad) + ggplot2::guides(color = "none")
  
  title_txt <- paste0(variant, ": ", str_to_pretty(region_zoom))
  track_plot <- trackplot_combine2(list(
    BPCells:::wrap_trackplot(atac_cov,      unit(0.17, "null")) + indiv_theme,
    BPCells:::wrap_trackplot(track_genes,      unit(0.17, "null")) + indiv_theme,
    BPCells:::wrap_trackplot(p.peak_delim,  unit(0.10, "null")) + indiv_theme,
    BPCells:::wrap_trackplot(pred,          unit(0.22, "null")) + indiv_theme,
    BPCells:::wrap_trackplot(Ref_contribs,  unit(0.17, "null")) + indiv_theme,
    BPCells:::wrap_trackplot(Alt_contribs,  unit(0.17, "null")) + indiv_theme,
    BPCells:::wrap_trackplot(track_hits,    unit(0.10, "null")) + indiv_theme
  ), title = title_txt)
  
  track_plot
  
  out_fig_path=paste0("~/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/cbpnet_examples/fig_out/", variant, ".svg")
  
  ggsave(out_fig_path, device = svglite::svglite, bg = "transparent", units = "px", height = 103, width = 173, dpi=300)
  
  ggsave(
    filename = out_fig_path,        # e.g., "track_plot.svg"
    plot     = track_plot + coord_cartesian(clip = "off"),
    device   = svglite::svglite,
    bg       = "transparent",
    width    = 173, height = 103,   # A5-ish landscape
    units    = "mm",
    limitsize = FALSE
  )
  
  ggsave(out_fig_path, track_plot, height = 6, width = 10, dpi=300)
}


# ggplot(vars_motifs) + geom_point(aes(-1*abs_IPS.mean, abs_logfc.mean))
# ggplot(vars_motifs) + geom_point(aes(abs_IPS.mean, abs_logfc.mean, color=hit_coefficient_global))
# ggplot(vars_motifs) + geom_point(aes(abs_IPS.mean, abs_logfc.mean, color=-log(abs_IPS.mean.pval)))
# 
# 
# ggplot(vars_motifs_top16) + geom_point(aes(hit_coefficient_global, abs_IPS.mean))
# ggplot(vars_motifs_top16) + geom_point(aes(abs_logfc.mean, hit_coefficient_global, ))
# 
# vars_motifs %>% filter(abs_IPS.mean<0.015 & abs_logfc.mean >0.5) %>% filter(hit_coefficient_global) arrange(desc(hit_coefficient_global))
#  
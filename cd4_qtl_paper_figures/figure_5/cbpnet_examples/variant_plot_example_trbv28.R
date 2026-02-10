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
# ----------------------- RESOURCES -------------------------------------------
# atrr: averaged variant contribution H5 across all folds
# data_h5: averaged variant prediction H5 across all folds
# tsv: output tsv from variant-scorer, or the input variant file, we need it to get the index of the variant and parse the predictions
# chrom_sizes: hg38 chromsizes file
# bw_path: BigWig from ATAC-seq data
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

transcripts <- read_gencode_genes("./references/", release="44")
head(transcripts)

atrr         <- "/gpfs/commons/home/mmatos/cd4_chrombpnet/chrombpnet_model_b7/variant_contribution_scores/averaged_cd4_tcells_AJ_common_variants.shap.counts.h5"
data_h5      <- "/gpfs/commons/home/mmatos/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_variant_prediction_scores.h5"
tsv          <- "/gpfs/commons/home/mmatos/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_scores/average_cd4_tcells_AJ_common_variants.mean.variant_scores.tsv"
chrom_sizes  <- "/gpfs/commons/home/mmatos/resources/genome/hg38.chrom.sizes"
motifs_path  <- "~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap/motifs_hits_filtered_c90_40seq_w_variant_hits_granges_1b.rds"
bw_path <- "/gpfs/commons/groups/lappalainen_lab/mmatos/cd4_aging_project/data/ATAC-seq/atac_preprocessing/nextflow_tmp/a4/4c07c68cd77cef99a8b320e8e66238/cd4_atac.bigWig"
filtered_peaks<-'/gpfs/commons/home/mmatos/cd4_chrombpnet/data/inputs/peaks/merged_cd4_samples_peaks_no_blacklist.sorted.bed' #peakcall used for ChromBPNet
dynamic_int<-read.csv("~/cd4_CellRegMap/002_interaction_analysis/results/results_01312025/res_34_5factors/2025_03_21_mofa_interaction_res34_significant_summary_fdr1.csv", 
                      row.names = 1)
var_prediction<-fread("~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap_hit_caller/unique_variants_overlapping_motifs.tsv")

motifs<-fread("~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap_hit_caller/unique_hits_filtered_motif_instances.tsv") %>% 
  select(c(start, end, variant_hit, family, strand, chr, variant_loc)) %>% 
  mutate(variant_loc2=paste0(chr, "_", variant_loc),
         motifs_disrupted=variant_hit,
         TF_family=family)

bed_peaks <- read_tsv(filtered_peaks, col_select = c(1:4), col_names = FALSE) 
colnames(bed_peaks) <- c("peak_chr", "peak_start", "peak_end", "peak_name")

var_prediction<-var_prediction %>% left_join(motifs, join_by(var_chr_pos==variant_loc2))

summary_df <- read_delim(
  "/gpfs/commons/home/mmatos/cd4_qtl_paper_figures/figure_1/data/eqtl_caqtl_finemapping_coloc_all.tsv"
)

chrombpnet_vars <- var_prediction$variant_id

## caQTL side
ca_all <- summary_df %>%
  mutate(
    var_id = gsub(
      "([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])",
      "chr\\1_\\2_\\3_\\4",
      variant_id.x
    )
  ) %>%
  filter(var_id %in% chrombpnet_vars) 

## eQTL side
eq_all <- summary_df %>%
  mutate(
    var_id = gsub(
      "([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])",
      "chr\\1_\\2_\\3_\\4",
      variant_id.y
    )
  ) %>%
  filter(var_id %in% chrombpnet_vars)

qtl_any <- bind_rows(ca_all, eq_all)  

rm(summary_df)
rm(ca_all)
rm(eq_all)

head(qtl_any)

qtl_any_dyn<-qtl_any %>% filter(gene %in% dynamic_int$gene)

qtl_any_dyn<-qtl_any_dyn %>% 
  left_join(var_prediction, join_by("var_id"=="variant_id")) %>% 
  select(var_id, finemapped_cs_coloc, gene, peak, variant_id.x, variant_id.y, PP.H4.abf, coloc_status,caQTL_variant, eQTL_variant, caQTL_variant,
         chr, pos,  allele1,  allele2, logfc.mean, jsd.mean, var_chr_pos, is_qtl, category,          
         peak_annotation, seqnames, variant_loc,  motifs_collapsed )

#############################################
# Plot example
#######################################
gene="TRBV28"
variant="chr7_142701355_G_C"
vp <- parse_variant_id(variant)
var_contig <- vp$chr
var_pos    <- vp$pos

region_window_broad <- 400000L
region_broad <- paste0(var_contig, ":", max(1, var_pos - region_window_broad), "-", min(var_pos + region_window_broad))
#region_ft <- list(chr = var_contig, start = max(1, var_pos - region_window_broad) , end = var_pos + region_window_broad)
region_broad.lsgwas <- str_to_gr(region_broad)

#1st zoom
region_window_broad2 <-10000L
region_broad_zoom2 <- paste0(var_contig, ":", max(1, var_pos - region_window_broad2), "-", var_pos + region_window_broad2)
region_broad.ls <- str_to_gr(region_broad_zoom2)

#2st zoom
region_window_broad3 <- 500L
region_broad_zoom3 <- paste0(var_contig, ":", max(1, var_pos - region_window_broad3), "-", var_pos + region_window_broad3)

region_window <- 50L
region_zoom <- paste0(var_contig, ":", max(1, var_pos - region_window), "-", var_pos + region_window)

region_variant<- paste0(var_contig, ":", max(1, var_pos+1), "-", var_pos)
###################
#caQTL finemapping
###### get caQTL finemap 
file_name_ca  = "/gpfs/commons/groups/lappalainen_lab/sghatan/marlis_pj/coloc/coloc_results/ca_eqtl_coloc/chr7_coloc_results.csv"
coloc_file_ca= gsub("_coloc_results","",gsub(".csv","",basename(file_name_ca)))

# Load colocalization results
coloc_results_ca = fread(file_name_ca) 


x=1031 #index of the variant in caQTL-GWAS coloc results
region=coloc_results_ca$region.x[x]
region=gsub("7:", "", region)
peak=coloc_results_ca$peak[x]

prep_ca <- prepare_caqtl_data(var_contig, peak, region, variant)
markers_ca <- prep_ca$markers
LD_ca <- prep_ca$LD
caqtl_variant = coloc_results_ca$caQTL_variant[x]

markers_ca <- markers_ca %>%
  mutate(
    pos_hg38 = as.integer(str_extract(marker, "(?<=:)\\d+")),
    z_pos0= as.numeric(scale(z_caqtl)) - min(as.numeric(scale(z_caqtl)), na.rm = TRUE) 
  ) 

position2_ca <- markers_ca$pos_hg38[markers_ca$marker == caqtl_variant]

minStart = max(1, var_pos - region_window_broad)
maxEnd = var_pos + region_window_broad
maxlimit = (ceiling(max((-1)*markers_ca$z_pos0 ) / 5) * 10) + 5


ld_colors <- c(
  "0.0–0.2" = "gray",
  "0.2–0.4" = "skyblue",
  "0.4–0.6" = "green",
  "0.6–0.8" = "orange",
  "0.8–1.0" = "red"
)
# ld_colors <- c(
#   "0.0–0.2" = "darkblue",
#   "0.2–0.4" = "skyblue",
#   "0.4–0.6" = "yellow",
#   "0.6–0.8" = "orange",
#   "0.8–1.0" = "red"
# )

# Identify index of lead SNP (use the one closest to snp_pos)
lead_index_ca <- which.min(abs(markers_ca$pos_hg38 - var_pos))
lead_snp_ca <- markers_ca$marker[lead_index_ca]
LD_ca1<-as.matrix(LD_ca)
# Extract correlation with lead SNP
ld_cor_ca <- LD_ca1[, lead_index_ca]
markers_ca$ld_r_gwas <- ld_cor_ca  # add it to markers
markers_ca$ld_r2_gwas <- markers_ca$ld_r_gwas^2

markers_ca$ld_bin_gwas <- cut(
  markers_ca$ld_r2_gwas,
  breaks = c(0, 0.2,  0.4, 0.6, 0.8, 1),
  labels = c("0.0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"),
  include.lowest = TRUE,
  right = FALSE
)
p.caqtl <- ggplot() +
  geom_point(
    data = markers_ca,
    aes(x = pos_hg38, y = z_pos0, fill = ld_bin_gwas),
    shape = 21, size = 2, alpha = 0.6, color = "gray20", stroke = 0.2
  ) +
  scale_fill_manual(values = ld_colors, name = expression(LD~(r^2))) +
  scale_x_continuous(expand = c(0, 0), limits = c(minStart, maxEnd)) +
  annotate(
    "text",
    x = minStart,
    y = max(markers_eqtl$z_pos0, na.rm = TRUE),
    label = paste0("caQTL ", caqtl_variant, " ", peak),
    hjust = 0, vjust = 1,
    size = 4
  ) +
  ylim(0, maxlimit) +
  labs(y = bquote(-log[10](p)), x = "") +
  BPCells:::trackplot_theme() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 11, angle = 90, vjust = 0.5)
  )


###################
#eQTL finemapping
###### get eQTL finemap 
file_name_ca  = "/gpfs/commons/groups/lappalainen_lab/sghatan/marlis_pj/coloc/coloc_results/ca_eqtl_coloc/chr7_coloc_results.csv"
coloc_file_ca= gsub("_coloc_results","",gsub(".csv","",basename(file_name_ca)))

# Load colocalization results
coloc_results_ca = fread(file_name_ca) 


x=1031 #index of the variant in caQTL-GWAS coloc results
region=coloc_results_ca$region.x[x]
region=gsub("7:", "", region)
gene=coloc_results_ca$gene[x]

prep_eqtl <- prepare_eqtl_data(var_contig, gene, region, variant)
markers_eqtl <- prep_eqtl$markers
LD_eqtl <- prep_eqtl$LD
eqtl_variant = coloc_results_ca$eQTL_variant[x]

markers_eqtl <- markers_eqtl %>%
  mutate(
    pos_hg38 = as.integer(str_extract(marker, "(?<=:)\\d+")),
    z_pos0= as.numeric(scale(z_eqtl)) - min(as.numeric(scale(z_eqtl)), na.rm = TRUE) 
  ) 

position2_e <- markers_eqtl$pos_hg38[markers_eqtl$marker == eqtl_variant]

minStart = max(1, var_pos - region_window_broad)
maxEnd = var_pos + region_window_broad
maxlimit = (ceiling(max((-1)*markers_eqtl$z_pos0 ) / 5) * 10) + 5

ld_colors <- c(
  "0.0–0.2" = "gray",
  "0.2–0.4" = "skyblue",
  "0.4–0.6" = "green",
  "0.6–0.8" = "orange",
  "0.8–1.0" = "red"
)
# ld_colors <- c(
#   "0.0–0.2" = "darkblue",
#   "0.2–0.4" = "skyblue",
#   "0.4–0.6" = "yellow",
#   "0.6–0.8" = "orange",
#   "0.8–1.0" = "red"
# )

# Identify index of lead SNP (use the one closest to snp_pos)
lead_index_eqtl <- which.min(abs(markers_eqtl$pos_hg38 - var_pos))
lead_snp_eqtl <- markers_eqtl$marker[lead_index_eqtl]
LD_e1<-as.matrix(LD_eqtl)
# Extract correlation with lead SNP
ld_cor_eqtl <- LD_e1[, lead_index_eqtl]
markers_eqtl$ld_r_gwas <- ld_cor_eqtl  # add it to markers
markers_eqtl$ld_r2_gwas <- markers_eqtl$ld_r_gwas^2

markers_eqtl$ld_bin_gwas <- cut(
  markers_eqtl$ld_r2_gwas,
  breaks = c(0, 0.2,  0.4, 0.6, 0.8, 1),
  labels = c("0.0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"),
  include.lowest = TRUE,
  right = FALSE
)
p.eqtl <- ggplot() +
  geom_point(
    data = markers_eqtl,
    aes(x = pos_hg38, y = z_pos0, fill = ld_bin_gwas),
    shape = 21, size = 2, alpha = 0.6, color = "gray20", stroke = 0.2
  ) +
  scale_fill_manual(values = ld_colors, name = expression(LD~(r^2))) +
  scale_x_continuous(expand = c(0, 0), limits = c(minStart, maxEnd)) +
  annotate(
    "text",
    x = minStart,
    y = max(markers_eqtl$z_pos0, na.rm = TRUE),
    label = paste0("eQTL ", eqtl_variant, " ", gene),
    hjust = 0, vjust = 1,
    size = 4
  ) +
  ylim(0, maxlimit) +
  labs(y = bquote(-log[10](p)), x = "") +
  BPCells:::trackplot_theme() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 11, angle = 90, vjust = 0.5)
  )


####################################
#
####################################
#motifs_with_tf_annotated_filt_c90_40seq_unique_variant_hits_granges.rds
motifs       <- readRDS(motifs_path)

#load observed ATAC coverage for the specified region
atac_sub <- import.bw(bw_path, which = region_broad.ls)

#----plot coverage ----------
atac_cov <- trackplot_bw(
  bw          = atac_sub,
  region      = region_broad_zoom2,
  facet_label = "ATAC",
  plot_as     = "area",
  tile_width  = 1,
  track_label = "ATAC",
  color       = "aquamarine4") +
  highlight_region(region_zoom, color = "skyblue", alpha = 0.2) +
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
p.peak_delim<-trackplot_peak_delim(bed_peaks_gr, region_broad_zoom2, color) + #bed_peaks_gr: rgranges object of peaks 
  highlight_region(region_zoom, color = "skyblue", alpha = 0.2) +
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
  region      = region_broad_zoom3,
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

motifs$variant_hit<-factor(motifs$variant_hit, levels = c(TRUE, FALSE)) #setting factor levels so color is consistent across variants

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



track_genes <- trackplot_gene_custom(transcripts, region_broad.ls) +
  ggplot2::guides(color = "none") +
  highlight_region(region_zoom, color = "skyblue", alpha = 0.2) +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, hjust = 0.5)) 

track_genes


title_txt <- paste0(variant, ": ", str_to_pretty(region_zoom))
track_plot <- trackplot_combine2(list(
  BPCells:::wrap_trackplot(p.caqtl,      unit(0.22, "null")) + indiv_theme ,
  BPCells:::wrap_trackplot(p.eqtl,      unit(0.22, "null")) + indiv_theme,
  BPCells:::wrap_trackplot(atac_cov,      unit(0.22, "null")) + indiv_theme,
  BPCells:::wrap_trackplot(p.peak_delim,  unit(0.10, "null")) + indiv_theme,
BPCells:::wrap_trackplot(track_genes,      unit(0.10, "null")) + indiv_theme,
  BPCells:::wrap_trackplot(pred,          unit(0.22, "null")) + indiv_theme,
  BPCells:::wrap_trackplot(Ref_contribs,  unit(0.17, "null")) + indiv_theme,
  BPCells:::wrap_trackplot(Alt_contribs,  unit(0.17, "null")) + indiv_theme,
  BPCells:::wrap_trackplot(track_hits,    unit(0.10, "null")) + indiv_theme
), title = title_txt)

track_plot

ggsave(
  filename = "~/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/cbpnet_examples/plots/TRBV28_coloc_qtl_chrom_dyna.pdf",
  plot     = track_plot,
  width    = 6, height = 6, units = "in"
)

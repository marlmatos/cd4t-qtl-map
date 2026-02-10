##############
# helper scripts for all trackplot examples and chromBPNet related figures
#  adapted from https://github.com/GreenleafLab/HDMA
#########

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(scales)
library(glue)
library(purrr)
library(stringr)
library(ggrepel)
library(ggseqlogo)
library(BPCells)
library(patchwork)
library(ggrastr)
library(BPCells)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(genio)
options(bitmapType = "cairo")

### Track plot helpers

prepareTrackData_grouped_genotype <- function(variant_id, metadata, base_dir, bigwig_dir, flanking = c(5000, 5000)) {
  chr_num <- as.numeric(strsplit(variant_id, ":")[[1]][1])
  
  # Read genotype files
  bim <- read_bim(paste0(base_dir, chr_num, '_ashkenazi.362.AF1.QC.BA.king2.hwe.annot.bim'))
  fam <- read_fam(paste0(base_dir, chr_num, '_ashkenazi.362.AF1.QC.BA.king2.hwe.annot.fam'))
  gen <- read_bed(
    paste0(base_dir, chr_num, '_ashkenazi.362.AF1.QC.BA.king2.hwe.annot.bed'),
    bim$id,
    gsub('-', '_', fam$id)
  )
  
  # Get bigWig paths
  bigwig_files <- list.files(path = bigwig_dir, pattern = "\\.bigWig$", full.names = TRUE)
  track_data_list <- setNames(as.list(bigwig_files), gsub("\\.bigWig$", "", basename(bigwig_files)))
  
  # Extract genotypes for the variant
  snp <- as.data.frame(gen[variant_id, ])
  colnames(snp) <- c(variant_id)
  snp$ind <- rownames(snp)
  
  # Merge genotype info with sample metadata
  snp_meta <- metadata %>%
    filter(WGS_sampleID %in% snp$ind) %>%
    merge(snp, ., by.x = 'ind', by.y = 'WGS_sampleID')
  
  # Join with bigWig paths
  track_df <- data.frame(
    ATAC_Sample_Name = names(track_data_list),
    bigwig_path = unlist(track_data_list)
  )
  
  snp_meta_tracks <- snp_meta %>%
    left_join(track_df, by = "ATAC_Sample_Name")
  
  # Build track_data
  track_data <- snp_meta_tracks %>%
    dplyr::select(ATAC_Sample_Name, bigwig_path, !!sym(variant_id)) %>%
    dplyr::rename(
      sample_id = ATAC_Sample_Name,
      bigWig = bigwig_path,
      colour_group = !!sym(variant_id)
    ) %>%
    mutate(
      track_id = "All Samples",
      colour_group = case_when(
        colour_group == 0 ~ "Homozygous Ref",
        colour_group == 1 ~ "Heterozygous",
        colour_group == 2 ~ "Homozygous Alt"
      ),
      colour_group = factor(colour_group, levels = c("Homozygous Ref", "Heterozygous", "Homozygous Alt")),
      scaling_factor = 1
    ) %>%
    filter(!is.na(colour_group))
  
  return(track_data)
}



#############
plotGroupedCoverage_allelic <- function(track_data,
                                gene_range,
                                fill_palette = c(
                                  "Homozygous Ref" = "darkgreen",
                                  "Heterozygous"   = "gray90",
                                  "Homozygous Alt" = "darkviolet"
                                ),
                                alpha = 0.6,
                                ymin_zero = TRUE) {
  
  # Internal helper to read and tidy one bigWig
  readCoverageFromBigWig <- function(bigwig_path, gene_range) {
    sel <- rtracklayer::BigWigSelection(gene_range)
    coverage_ranges <- rtracklayer::import.bw(bigwig_path, selection = sel)
    GenomeInfoDb::seqlevels(coverage_ranges) <- S4Vectors::as.vector.Rle(GenomicRanges::seqnames(gene_range), mode = "character")
    coverage_rle <- GenomicRanges::coverage(coverage_ranges, weight = GenomicRanges::score(coverage_ranges))[[1]]
    region_start <- start(gene_range)
    region_end <- end(gene_range)
    coverage_rle <- coverage_rle[region_start:region_end]
    data.frame(position = seq(region_start, region_end),
               coverage = as.numeric(coverage_rle))
  }
  
  # Read and annotate all bigWigs
  all_cov <- purrr::pmap_dfr(
    list(bigwig = track_data$bigWig,
         sample_id = track_data$sample_id,
         colour_group = track_data$colour_group),
    function(bigwig, sample_id, colour_group) {
      readCoverageFromBigWig(bigwig, gene_range) %>%
        dplyr::mutate(sample_id = sample_id, colour_group = colour_group)
    }
  )
  
  # Normalize by library size
  all_cov <- dplyr::left_join(all_cov, track_data[, c("sample_id", "scaling_factor")], by = "sample_id") %>%
    dplyr::mutate(coverage = coverage / scaling_factor)
  
  # Mean coverage per group
  grouped_cov <- all_cov %>%
    dplyr::group_by(position, colour_group) %>%
    dplyr::summarise(mean_coverage = mean(coverage), .groups = "drop")
  
  # Compute Y range
  ymax <- max(grouped_cov$mean_coverage)*1.05
  ymax_accuracy <- 10^as.integer(log10(0.01 * ymax))
  
  if (!ymin_zero) {
    ymin <- min(grouped_cov$mean_coverage)
    ymin_accuracy <- 10^as.integer(log10(0.01 * abs(ymin)))
    range_label <- glue::glue("[{scales::label_comma(accuracy = ymin_accuracy)(ymin)}-{scales::label_comma(accuracy = ymax_accuracy)(ymax)}]")
  } else {
    ymin <- 0
    range_label <- sprintf("[0-%s]", scales::label_comma(accuracy = ymax_accuracy)(ymax))
  }
  
  # Plot
  ggplot(grouped_cov, aes(x = position, y = mean_coverage, color = colour_group)) +
    geom_line(alpha = alpha, size=1) +
    scale_color_manual(values = fill_palette) +
    scale_x_continuous(limits = c(start(gene_range), end(gene_range)), expand = c(0, 0), labels = scales::label_comma()) +
    scale_y_continuous(limits = c(ymin, ymax), expand = c(0, 0)) +
    annotate("text", x = start(gene_range), y = ymax, label = range_label, vjust = 1.5, hjust = -0.1, size = 3.5) +
    labs(x = "Genomic Position (bp)", y = "Mean Normalized Coverage", color = "Genotype") +
    BPCells:::trackplot_theme()
}


##filled in

plotGroupedCoverage_allelic_fill <- function(track_data,
                                             gene_range,
                                             fill_palette = c(
                                               "Homozygous Ref" = "darkgreen",
                                               "Heterozygous"   = "gray90",
                                               "Homozygous Alt" = "darkviolet"
                                             ),
                                             genotype_levels = NULL,   # bottom -> top
                                             alpha = 0.6,
                                             ymin_zero = TRUE,
                                             rasterize = TRUE,
                                             raster_dpi = 400) {
  
  readCoverageFromBigWig <- function(bigwig_path, gene_range) {
    sel <- rtracklayer::BigWigSelection(gene_range)
    coverage_ranges <- rtracklayer::import.bw(bigwig_path, selection = sel)
    GenomeInfoDb::seqlevels(coverage_ranges) <-
      S4Vectors::as.vector.Rle(GenomicRanges::seqnames(gene_range), mode = "character")
    coverage_rle <- GenomicRanges::coverage(
      coverage_ranges,
      weight = GenomicRanges::score(coverage_ranges)
    )[[1]]
    rs <- start(gene_range); re <- end(gene_range)
    coverage_rle <- coverage_rle[rs:re]
    data.frame(position = seq(rs, re), coverage = as.numeric(coverage_rle))
  }
  
  all_cov <- purrr::pmap_dfr(
    list(bigwig = track_data$bigWig,
         sample_id = track_data$sample_id,
         colour_group = track_data$colour_group),
    function(bigwig, sample_id, colour_group) {
      readCoverageFromBigWig(bigwig, gene_range) %>%
        dplyr::mutate(sample_id = sample_id, colour_group = colour_group)
    }
  ) %>%
    dplyr::left_join(track_data[, c("sample_id", "scaling_factor")], by = "sample_id") %>%
    dplyr::mutate(coverage = coverage / scaling_factor)
  
  grouped_cov <- all_cov %>%
    dplyr::group_by(position, colour_group) %>%
    dplyr::summarise(mean_coverage = mean(coverage, na.rm = TRUE), .groups = "drop")
  
  # --- Order groups so small area is on top ---
  if (is.null(genotype_levels)) {
    lvls <- grouped_cov %>%
      dplyr::group_by(colour_group) %>%
      dplyr::summarise(total_area = sum(mean_coverage, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(total_area)) %>%  # largest first (drawn first -> bottom)
      dplyr::pull(colour_group)
  } else {
    lvls <- genotype_levels
  }
  
  grouped_cov <- grouped_cov %>%
    dplyr::mutate(colour_group = factor(colour_group, levels = lvls))
  
  ymax <- max(grouped_cov$mean_coverage, na.rm = TRUE) * 1.05
  ymax_accuracy <- 10^as.integer(log10(0.01 * ymax))
  
  if (!ymin_zero) {
    ymin <- min(grouped_cov$mean_coverage, na.rm = TRUE)
    ymin_accuracy <- 10^as.integer(log10(0.01 * abs(ymin)))
    range_label <- glue::glue(
      "[{scales::label_comma(accuracy = ymin_accuracy)(ymin)}-{scales::label_comma(accuracy = ymax_accuracy)(ymax)}]"
    )
  } else {
    ymin <- 0
    range_label <- sprintf("[0-%s]", scales::label_comma(accuracy = ymax_accuracy)(ymax))
  }
  
  # 1) build the heavy layer (optionally rasterised)
  lyr <- ggplot2::geom_area(
    ggplot2::aes(x = position, y = mean_coverage, fill = colour_group),
    alpha = alpha,
    position = "identity",
    na.rm = TRUE
  )
  
  if (rasterize) {
    # requires: library(ggrastr)
    lyr <- ggrastr::rasterise(lyr, dpi = raster_dpi)
  }
  
  # 2) assemble the plot
  ggplot2::ggplot(grouped_cov) +
    lyr +
    ggplot2::scale_fill_manual(values = fill_palette[levels(grouped_cov$colour_group)]) +
    ggplot2::scale_x_continuous(
      limits = c(start(gene_range), end(gene_range)),
      expand = c(0, 0),
      labels = scales::label_comma()
    ) +
    ggplot2::scale_y_continuous(limits = c(ymin, ymax), expand = c(0, 0)) +
    ggplot2::annotate(
      "text",
      x = start(gene_range), y = ymax,
      label = range_label,
      vjust = 1.5, hjust = -0.1,
      size = 3.5
    ) +
    ggplot2::labs(
      x = "Genomic Position (bp)",
      y = "Mean Normalized Coverage",
      fill = "Genotype"
    ) +
    BPCells:::trackplot_theme()
}




#######
prepare_manhattan_gwas_eqtl_data <- function(x,
                                             coloc_results,
                                             gwas_file) {  
  print("making plot")
  
  eqtl_variant = coloc_results$variant_id_GRC38[x]
  eqtl_name = "All_CD4T_cells"
  disease_name = str_extract(file_name, "[^_/]+_[^_/]+(?=_preprocessed)")
  region = str_split_fixed(coloc_results$region_GRC38[x], ":",2)[,2]
  chr = strsplit(eqtl_variant,":", 4)[[1]][1]
  gene_name = coloc_results$gene[x]
  
  # Load GWAS sumstats
  gwas_stats = fread(paste0('/gcgnl/finemapping_autoimmune/data/preprocessed_v1/',gwas_file,'.tsv')) %>%
    mutate(variant_id_GRC37 = paste(chromosome, position,sep=":"))
  
  # Load eqtl
  eqtl_lbf = fread(paste0("/gcgl/lappalainen_lab/sghatan/marlis_pj/coloc/SuSiE_finemap_lbf_GRC37/",eqtl_name,"/",eqtl_name,"_chr",chr,"_lbfs_GRC37.txt")) %>%
    mutate(variant_id_GRC37 = paste(gsub("chr","",chr), pos_GRC37, sep = ":")) %>%
    filter(gene == gene_name)
  
  # Load LD matrix
  ld_dir = paste0('/gcgl/sghatan/marlis_pj/coloc/Ashkenazi_LD_matrices/chr',chr,'/',region,'/',region,'_imputed')
  LD <- fread(paste0(ld_dir, ".ld"))
  BIM <- fread(paste0(ld_dir, ".bim"))
  
  # assign headers
  setnames(BIM, c("chr", "rsid", "dk", "pos", "a0", "a1"))
  BIM[,chr_pos_GRC38 := paste(chr, pos, sep = ":")]
  
  # assign SNP labels to LD matrix
  setnames(LD, BIM$rsid)
  LD[, variant_id_GRC38 := BIM$rsid]
  
  # Remove columns from the LD matrix that are all NAs
  na_columns = colSums(is.na(LD)) == nrow(LD)
  nacol_names = names(na_columns[na_columns == TRUE])
  LD = LD[, !names(LD) %in% nacol_names, with = FALSE]
  # Filter rows with complete cases
  LD = LD[complete.cases(LD)]
  
  # Add GRC37 variant id to bim file
  BIM = BIM %>% 
    left_join(eqtl_lbf[,c("variant_id_GRC38", "variant_id_GRC37")], by = c("rsid" = "variant_id_GRC38")) %>%
    distinct(variant_id_GRC37, .keep_all=T) %>% drop_na()
  
  # Find variants across all three lists
  intermediate_result <- Reduce(intersect, list(
    eqtl_lbf$variant_id_GRC37,
    gwas_stats$variant_id_GRC37,
    BIM$variant_id_GRC37
  ))
  
  # Filter BIM file to this list
  BIM = BIM[variant_id_GRC37 %in% intermediate_result,]
  
  # Filter snps in LD matrix for those intersecting
  LD <- LD[LD$variant_id_GRC38 %in% BIM$rsid,]
  
  # Make sure LD only contains columns in LD$SNP and "SNP"
  required_cols <- unique(c(LD$variant_id_GRC38,"variant_id_GRC38"))
  LD <- LD[, ..required_cols]
  # Remove duplicated SNPs
  LD_clean <- unique(LD, by = "variant_id_GRC38")
  
  eqtl_region = eqtl_lbf[variant_id_GRC37 %in% intermediate_result,] %>% 
    distinct(variant_id_GRC37, .keep_all = T) %>%
    arrange(variant_id_GRC37) %>%
    mutate(z = beta/slope_se)
  
  gwas_region = gwas_stats[variant_id_GRC37 %in% intermediate_result,] %>%
    arrange(variant_id_GRC37) %>%
    mutate(z = beta/se)
  
  identical(gwas_region$variant_id_GRC37, eqtl_region$variant_id_GRC37)
  
  markers <- data.frame(marker = eqtl_region$variant_id_GRC38,
                        chr = chr,
                        pos = gwas_region$position,
                        z_1 = gwas_region$z,
                        z_2 = eqtl_region$z)
  
  LD = LD_clean %>% dplyr::select(-variant_id_GRC38)
  
  markers <- markers[match(colnames(LD), markers$marker), ]
  
  identical(colnames(LD), markers$marker)
  
  return(list(markers = markers, LD = LD))
}

######

prepare_manhattan_gwas_caqtl_data <- function(x,
                                             coloc_results,
                                             gwas_file) {  
  print("making plot")
  
  caqtl_variant = coloc_results$variant_id_GRC38[x]
  caqtl_name = "CD4T_chromatin"
  disease_name = str_extract(file_name, "[^_/]+_[^_/]+(?=_preprocessed)")
  region = str_split_fixed(coloc_results$region_GRC38[x], ":",2)[,2]
  chr = strsplit(caqtl_variant,":", 4)[[1]][1]
  peak_name = coloc_results$peak[x]
  
  gwas_file_name <- str_remove(gwas_file, "^CD4T_chromatin_")  
  # Load GWAS sumstats
  gwas_stats = fread(paste0('/gcgnl/finemapping_autoimmune/data/preprocessed_v1/',gwas_file_name, ".tsv")) %>%
    mutate(variant_id_GRC37 = paste(chromosome, position,sep=":"))
  
  # Load eqtl
  caqtl_lbf = fread(paste0("/gcgl/sghatan/marlis_pj/coloc/SuSiE_finemap_lbf_GRC37/",caqtl_name,"/",caqtl_name,"_chr",chr,"_lbfs_GRC37.txt")) %>%
    mutate(variant_id_GRC37 = paste(gsub("chr","",chr), pos_GRC37, sep = ":")) %>%
    filter(peak == peak_name)
  
  fread("/gcgl/sghatan/marlis_pj/coloc/SuSiE_finemap_lbf_GRC37/CD4T_chromatin/CD4T_chromatin_chr6_lbfs_GRC37.txt")
  # Load LD matrix
  ld_dir = paste0('/gcgl/sghatan/marlis_pj/coloc/Ashkenazi_LD_matrices/',"chr", chr,'/',region,'/',region,'_imputed')
  LD <- fread(paste0(ld_dir, ".ld"))
  BIM <- fread(paste0(ld_dir, ".bim"))
  
  # assign headers
  setnames(BIM, c("chr", "rsid", "dk", "pos", "a0", "a1"))
  BIM[,chr_pos_GRC38 := paste(chr, pos, sep = ":")]
  
  # assign SNP labels to LD matrix
  setnames(LD, BIM$rsid)
  LD[, variant_id_GRC38 := BIM$rsid]
  
  # Remove columns from the LD matrix that are all NAs
  na_columns = colSums(is.na(LD)) == nrow(LD)
  nacol_names = names(na_columns[na_columns == TRUE])
  LD = LD[, !names(LD) %in% nacol_names, with = FALSE]
  # Filter rows with complete cases
  LD = LD[complete.cases(LD)]
  
  # Add GRC37 variant id to bim file
  BIM = BIM %>% 
    left_join(caqtl_lbf[,c("variant_id_GRC38", "variant_id_GRC37")], by = c("rsid" = "variant_id_GRC38")) %>%
    distinct(variant_id_GRC37, .keep_all=T) %>% drop_na()
  
  # Find variants across all three lists
  intermediate_result <- Reduce(intersect, list(
    caqtl_lbf$variant_id_GRC37,
    gwas_stats$variant_id_GRC37,
    BIM$variant_id_GRC37
  ))
  
  # Filter BIM file to this list
  BIM = BIM[variant_id_GRC37 %in% intermediate_result,]
  
  # Filter snps in LD matrix for those intersecting
  LD <- LD[LD$variant_id_GRC38 %in% BIM$rsid,]
  
  # Make sure LD only contains columns in LD$SNP and "SNP"
  required_cols <- unique(c(LD$variant_id_GRC38,"variant_id_GRC38"))
  LD <- LD[, ..required_cols]
  # Remove duplicated SNPs
  LD_clean <- unique(LD, by = "variant_id_GRC38")
  
  caqtl_region = caqtl_lbf[variant_id_GRC37 %in% intermediate_result,] %>% 
    distinct(variant_id_GRC37, .keep_all = T) %>%
    arrange(variant_id_GRC37) %>%
    mutate(z = beta/slope_se)
  
  gwas_region = gwas_stats[variant_id_GRC37 %in% intermediate_result,] %>%
    arrange(variant_id_GRC37) %>%
    mutate(z = beta/se)
  
  identical(gwas_region$variant_id_GRC37, caqtl_region$variant_id_GRC37)
  
  markers <- data.frame(marker = caqtl_region$variant_id_GRC38,
                        chr = chr,
                        pos = gwas_region$position,
                        z_1 = gwas_region$z,
                        z_2 = caqtl_region$z)
  
  LD = LD_clean %>% dplyr::select(-variant_id_GRC38)
  
  markers <- markers[match(colnames(LD), markers$marker), ]
  
  identical(colnames(LD), markers$marker)
  
  return(list(markers = markers, LD = LD))
}


#####
prepare_caqtl_data <- function(chr, peak_name, region, variant_id = NULL) {  
  print("preparing caQTL data")
  
  caqtl_name = "CD4T_chromatin"
  # Clean chromosome format (remove "chr" if present)
  chr <- gsub("chr", "", chr)
  
  # Load caQTL
  caqtl_lbf = fread(paste0("/gcgl/sghatan/marlis_pj/coloc/SuSiE_finemap_lbf_GRC37/",caqtl_name,"/",caqtl_name,"_chr",chr,"_lbfs_GRC37.txt")) %>%
    mutate(variant_id_GRC37 = paste(gsub("chr","",chr), pos_GRC37, sep = ":")) %>% filter(peak==peak_name)
  
  # Load LD matrix
  ld_dir = paste0('/gcgl/sghatan/marlis_pj/coloc/Ashkenazi_LD_matrices/chr',chr,'/',region,'/',region,'_imputed')
  LD <- fread(paste0(ld_dir, ".ld"))
  BIM <- fread(paste0(ld_dir, ".bim"))
  
  # assign headers
  setnames(BIM, c("chr", "rsid", "dk", "pos", "a0", "a1"))
  BIM[,chr_pos_GRC38 := paste(chr, pos, sep = ":")]
  
  # assign SNP labels to LD matrix
  setnames(LD, BIM$rsid)
  LD[, variant_id_GRC38 := BIM$rsid]
  
  # Remove columns from the LD matrix that are all NAs
  na_columns = colSums(is.na(LD)) == nrow(LD)
  nacol_names = names(na_columns[na_columns == TRUE])
  LD = LD[, !names(LD) %in% nacol_names, with = FALSE]
  # Filter rows with complete cases
  LD = LD[complete.cases(LD)]
  
  # Add GRC37 variant id to bim file
  BIM = BIM %>% 
    left_join(caqtl_lbf[,c("variant_id_GRC38", "variant_id_GRC37")], by = c("rsid" = "variant_id_GRC38")) %>%
    distinct(variant_id_GRC37, .keep_all=T) %>% 
    drop_na()
  
  # Find variants present in both caQTL and LD data
  intermediate_result <- intersect(
    caqtl_lbf$variant_id_GRC37,
    BIM$variant_id_GRC37
  )
  
  # Filter BIM file to this list
  BIM = BIM[variant_id_GRC37 %in% intermediate_result,]
  
  # Filter snps in LD matrix for those intersecting
  LD <- LD[LD$variant_id_GRC38 %in% BIM$rsid,]
  
  # Make sure LD only contains columns in LD$SNP and "SNP"
  required_cols <- unique(c(LD$variant_id_GRC38,"variant_id_GRC38"))
  LD <- LD[, ..required_cols]
  # Remove duplicated SNPs
  LD_clean <- unique(LD, by = "variant_id_GRC38")
  
  # Prepare caQTL data
  caqtl_region = caqtl_lbf[variant_id_GRC37 %in% intermediate_result,] %>% 
    distinct(variant_id_GRC37, .keep_all = T) %>%
    arrange(variant_id_GRC37) %>%
    mutate(z = beta/slope_se)
  
  # Create markers data frame with caQTL data only
  markers <- data.frame(
    marker = caqtl_region$variant_id_GRC38,
    chr = rep(paste0("chr", chr), nrow(caqtl_region)),
    pos = caqtl_region$pos_GRC37,  # Use GRC37 position (this exists in caQTL data)
    z_caqtl = caqtl_region$z,
    beta_caqtl = caqtl_region$beta,
    se_caqtl = caqtl_region$slope_se,
    pval_caqtl = caqtl_region$pval_nominal
  )
  
  # Clean LD matrix
  LD = LD_clean %>% dplyr::select(-variant_id_GRC38)
  
  # Match markers to LD matrix order
  markers <- markers[match(colnames(LD), markers$marker), ]
  
  # Verify consistency
  if(!identical(colnames(LD), markers$marker)) {
    warning("LD matrix columns don't match marker order")
  }
  
  return(list(
    markers = markers, 
    LD = LD,
    peak_name = peak_name,
    variant_id = variant_id,
    chr = chr,
    region = region
  ))
}


##### eqtl
prepare_eqtl_data <- function(chr, gene_name, region, variant_id = NULL) {  
  print("preparing eQTL data")
  
  caqtl_name = "All_CD4T_cells"
  # Clean chromosome format (remove "chr" if present)
  chr <- gsub("chr", "", chr)
  
  # Load caQTL
  eqtl_lbf = fread(paste0("/gcgl/sghatan/marlis_pj/coloc/SuSiE_finemap_lbf_GRC37/",caqtl_name,"/",caqtl_name,"_chr",chr,"_lbfs_GRC37.txt")) %>%
    mutate(variant_id_GRC37 = paste(gsub("chr","",chr), pos_GRC37, sep = ":")) %>% filter(gene==gene_name)
  
  # Load LD matrix
  ld_dir = paste0('/gcgl/sghatan/marlis_pj/coloc/Ashkenazi_LD_matrices/chr',chr,'/',region,'/',region,'_imputed')
  LD <- fread(paste0(ld_dir, ".ld"))
  BIM <- fread(paste0(ld_dir, ".bim"))
  
  # assign headers
  setnames(BIM, c("chr", "rsid", "dk", "pos", "a0", "a1"))
  BIM[,chr_pos_GRC38 := paste(chr, pos, sep = ":")]
  
  # assign SNP labels to LD matrix
  setnames(LD, BIM$rsid)
  LD[, variant_id_GRC38 := BIM$rsid]
  
  # Remove columns from the LD matrix that are all NAs
  na_columns = colSums(is.na(LD)) == nrow(LD)
  nacol_names = names(na_columns[na_columns == TRUE])
  LD = LD[, !names(LD) %in% nacol_names, with = FALSE]
  # Filter rows with complete cases
  LD = LD[complete.cases(LD)]
  
  # Add GRC37 variant id to bim file
  BIM = BIM %>% 
    left_join(eqtl_lbf[,c("variant_id_GRC38", "variant_id_GRC37")], by = c("rsid" = "variant_id_GRC38")) %>%
    distinct(variant_id_GRC37, .keep_all=T) %>% 
    drop_na()
  
  # Find variants present in both caQTL and LD data
  intermediate_result <- intersect(
    eqtl_lbf$variant_id_GRC37,
    BIM$variant_id_GRC37
  )
  
  # Filter BIM file to this list
  BIM = BIM[variant_id_GRC37 %in% intermediate_result,]
  
  # Filter snps in LD matrix for those intersecting
  LD <- LD[LD$variant_id_GRC38 %in% BIM$rsid,]
  
  # Make sure LD only contains columns in LD$SNP and "SNP"
  required_cols <- unique(c(LD$variant_id_GRC38,"variant_id_GRC38"))
  LD <- LD[, ..required_cols]
  # Remove duplicated SNPs
  LD_clean <- unique(LD, by = "variant_id_GRC38")
  
  # Prepare caQTL data
  eqtl_region = eqtl_lbf[variant_id_GRC37 %in% intermediate_result,] %>% 
    distinct(variant_id_GRC37, .keep_all = T) %>%
    arrange(variant_id_GRC37) %>%
    mutate(z = beta/slope_se)
  
  # Create markers data frame with caQTL data only
  markers <- data.frame(
    marker = eqtl_region$variant_id_GRC38,
    chr = rep(paste0("chr", chr), nrow(eqtl_region)),
    pos = eqtl_region$pos_GRC37,  # Use GRC37 position (this exists in caQTL data)
    z_eqtl = eqtl_region$z,
    beta_caqtl = eqtl_region$beta,
    se_caqtl = eqtl_region$slope_se,
    pval_caqtl = eqtl_region$pval_nominal
  )
  
  # Clean LD matrix
  LD = LD_clean %>% dplyr::select(-variant_id_GRC38)
  
  # Match markers to LD matrix order
  markers <- markers[match(colnames(LD), markers$marker), ]
  
  # Verify consistency
  if(!identical(colnames(LD), markers$marker)) {
    warning("LD matrix columns don't match marker order")
  }
  
  return(list(
    markers = markers, 
    LD = LD,
    gene_name = gene_name,
    variant_id = variant_id,
    chr = chr,
    region = region
  ))
}

####
# Annotation helpers -----------------------------------------------------------

#' Highlight a small genomic region on a trackplot
#' 
#' @param region character, specifying region to highlight
#'
#' @examples 
#' trackplot_coverage(
#'     region = region,
#'     fragments= bp_obj$frags,
#'     groups = bp_obj$cluster_metadata$Cluster) +
#'     highlight_region(region_small)
highlight_region <- function(region, alpha = 0.2, color = "red") {
  
  region_gr <- str_to_gr(region)
  
  ggplot2::annotate("rect",
                    alpha = alpha,
                    xmin = start(region_gr), xmax = end(region_gr),
                    ymin = -Inf, ymax = Inf,
                    fill = color)
  
}

highlight_relative_region <- function(start, end, alpha = 0.2, color = "red") {
  
  ggplot2::annotate("rect",
                    alpha = alpha,
                    xmin = start, xmax = end,
                    ymin = -Inf, ymax = Inf,
                    fill = color)
  
}

####
#' Convert strings representing genomic regions to GRanges
#' 
#' @param regions character, vector of regions in the form chr1:10000-20000
#' @param sep character, two delimiters which separate chrom, start, end in the
#' strings provided to \code{regions}
#' 
#' @example 
#' str_to_gr("chr15:52785497-52791921")
str_to_gr <- function(regions, sep = c(":", "-")) {
  
  ranges.df <- data.frame(ranges = regions)
  separate(
    data = ranges.df,
    col = "ranges",
    sep = paste0(sep[[1]], "|", sep[[2]]),
    into = c("chr", "start", "end")
  ) %>% 
    GRanges()
  
}

#' Make a prettier version of the coordinate string, using commas to separate 1kbs
#' so that it's more readable
#' 
#' @param region character, region in the form of chr1:10000-20000
#' 
#' @value "chr1:10,000-20,000"
str_to_pretty <- function(region) {
  
  region_gr <- str_to_gr(region)
  paste0(seqlevels(region_gr), ":", scales::comma(start(region_gr)), "-", scales::comma(end(region_gr)))
  
}


####
plot_manhattan <- function(markers, snp_pos, z_col = "z_1", flank = 2e5) {
  minStart <- snp_pos - flank
  maxEnd <- snp_pos + flank
  y_vals <- markers[[z_col]]
  maxlimit <- (ceiling(max(y_vals, na.rm = TRUE) / 5) * 5) + 0.5
  
  ggplot(markers, aes(x = pos, y = .data[[z_col]])) +
    geom_point(color = "grey50", size = 0.8, alpha = 0.5) +
    scale_x_continuous(limits = c(minStart, maxEnd), expand = c(0, 0)) +
    ylim(0, maxlimit) +
    labs(y = bquote(beta / SE), x = NULL) +
    BPCells:::trackplot_theme() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 11, angle = 90, vjust = 0.5)
    )
}


indiv_theme = theme(plot.margin = unit(c(0, 0, 0, 0), "pt"),
                    legend.position = "none")
#' Combine track plots
#' 
#' Combines multiple track plots of the same region into a single grid.
#' Uses the `patchwork` package to perform the alignment.
#'
#' @param tracks List of tracks in order from top to bottom, generally ggplots as output from
#'    the other `trackplot_*()` functions.
#' @param side_plot Optional plot to align to the right (e.g. RNA expression per cluster). Will be aligned to the first
#'    `trackplot_coverage()` output if present, or else the first generic ggplot in the alignment. Should be in horizontal orientation and 
#'    in the same cluster ordering as the coverage plots.
#' @param title Text for overarching title of the plot
#' @param side_plot_width Fraction of width that should be used for the side plot relative to the main track area
#' @return A plot object with aligned genome plots. Each aligned row has
#'    the text label, y-axis, and plot body. The relative height of each row is given
#'    by heights. A shared title and x-axis are put at the top.
#' @seealso `trackplot_coverage()`, `trackplot_gene()`, `trackplot_loop()`, `trackplot_scalebar()`
#' @export
trackplot_combine2 <- function(tracks, side_plot = NULL, title = NULL, side_plot_width = 0.3) {
  for (plot in tracks) {
    BPCells:::assert_is(plot, "ggplot")
  }
  if (!is.null(side_plot)) {
    BPCells:::assert_is(side_plot, "ggplot")
  }
  
  # Calculate layout information on the plots
  heights <- list()
  collapse_upper_margin <- rep.int(TRUE, length(tracks))
  side_plot_row <- NULL
  areas <- NULL
  last_region <- NULL
  for (i in seq_along(tracks)) {
    if (is(tracks[[i]], "trackplot")) {
      if (tracks[[i]]$trackplot$takes_sideplot && is.null(side_plot_row)) {
        side_plot_row <- i
      }
      # If we switch regions, don't collapse margins into the track above
      if (!is.null(tracks[[i]]$trackplot$region)) {
        if (!is.null(last_region) && last_region != tracks[[i]]$trackplot$region && i > 1) collapse_upper_margin[i] <- FALSE 
        last_region <- tracks[[i]]$trackplot$region
      }
      
      # Preserve top and bottom margins if `keep_vertical_margin`
      if (tracks[[i]]$trackplot$keep_vertical_margin) {
        collapse_upper_margin[i] <- FALSE
        if (i < length(tracks)) collapse_upper_margin[i+1] <- FALSE
      }
    } else {
      if (is.null(side_plot_row)) side_plot_row <- i
    }
    heights <- c(heights, list(get_trackplot_height(tracks[[i]])))
    areas <- c(areas, list(patchwork::area(i, 1)))
  }
  heights <- do.call(grid::unit.c, heights)
  if (!is.null(side_plot) && is.null(side_plot_row)) {
    rlang::warn("Did not find a row to place the side_plot: no trackplot_coverage() or base ggplot tracks found. Defaulting to first row")
    side_plot_row <- 1L
  }
  
  # Collapse margins as needed among plots
  for (i in seq_along(tracks)) {
    plot.margin <- c(TRUE, TRUE, TRUE, TRUE) # Top, right, bottom, left
    if (!is.null(side_plot)) plot.margin[2] <- FALSE # Side plot should be flush on right side
    if (i < length(tracks) && collapse_upper_margin[i+1]) plot.margin[3] <- FALSE # Plot below should be flush
    if (collapse_upper_margin[i]) plot.margin[1] <- FALSE
    
    if (!plot.margin[3]) {
      tracks[[i]] <- tracks[[i]] +
        ggplot2::guides(x="none") +
        ggplot2::labs(x=NULL)
    }
    
    # Independent of showing the axis, we'll remove the bottom margin if the next row has the side_plot, since the
    # axis tick labels will add in some natural margin already
    # TODO: raise issue in BPCells
    if (!is.null(side_plot) && !is.null(side_plot_row) && (i+1 == side_plot_row)) plot.margin[3] <- FALSE  # ADDED CHECK HERE
    
    
    tracks[[i]] <- tracks[[i]] + ggplot2::theme(plot.margin=ggplot2::unit(5.5*plot.margin, "pt"))
  }
  
  # Reduce cut-off y-axis labels. Put plots with y axis labels later in the plot list, as they will take layer priority with patchwork
  has_y_axis <- vapply(tracks, function(t) is(t, "ggplot") && !is.null(t$labels$y), logical(1))
  tracks <- c(tracks[!has_y_axis], tracks[has_y_axis])
  areas <- c(areas[!has_y_axis], areas[has_y_axis])
  
  if (is.null(side_plot)) {
    widths <- c(1)
  } else {
    # Decide whether to put legends below/above side plot by adding up the height of all relatively-sized tracks
    height_above <- sum(as.vector(heights)[seq_along(heights) < side_plot_row & grid::unitType(heights) == "null"])
    height_below <- sum(as.vector(heights)[seq_along(heights) > side_plot_row & grid::unitType(heights) == "null"])
    if (height_above < height_below) {
      guide_position <- patchwork::area(side_plot_row+1L, 2, length(tracks))
    } else {
      guide_position <- patchwork::area(1L, 2, side_plot_row-1L)
    }
    
    widths <- c(1, side_plot_width)
    areas <- c(areas, list(patchwork::area(side_plot_row, 2), guide_position))
    # Make adjustments to the side plot style to fit in with tracks
    side_plot <- side_plot + 
      ggplot2::scale_x_discrete(limits=rev, position="top") +
      ggplot2::scale_y_continuous(position="right") +
      ggplot2::coord_flip() +
      ggplot2::labs(x=side_plot$labels$y, y=NULL) +
      ggplot2::theme(
        plot.margin=ggplot2::unit(c(0,0,0,0), "pt"),
        axis.ticks.length.x.top=grid::unit(-2.75, "pt")
      ) 
    guide_area <- patchwork::guide_area() + ggplot2::theme(plot.margin=ggplot2::unit(c(0,0,0,0), "pt"))
    tracks <- c(tracks, list(side_plot, guide_area))
  }
  
  patch <- patchwork::wrap_plots(tracks) +
    patchwork::plot_layout(ncol = 1, byrow = FALSE, heights = heights, widths=widths, guides = "collect", design=do.call(c, areas))
  
  if (!is.null(side_plot)) {
    # If a side plot is present, switch legend layout to use the horizontal space better
    patch <- patch * ggplot2::theme(
      legend.direction="horizontal", 
      legend.title.position = "top"
    )
  }
  if (!is.null(title)) {
    patch <- patch + patchwork::plot_annotation(title = title, theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }
  return(patch)
}


#' Make a trackplot annotating peaks, that's compatible with BPCells plots.
#' 
#' In these plots, each hit in the region is represented by a rectangle, and plot
#' alongside its name.
#'
#' @param hits GRanges object obtained from \code{rtracklayer::import.bed(hits_bed_path)},
#' corresponding to named motif hits.
#' @param region character, genomic coordinates in the form chr:start-end
#' @param track_label character, the y-axis label for the track
#' @param facet_label character, the facet label for the track
#' 
#' @value
#' BPCells trackplot as returned by \code{BPCells:::wrap_trackplot}
trackplot_peak_delim <- function(hits_gr, region, color, track_label = "Peaks",
                           facet_label = NULL, rel_height = 0.6,
                           label_size = 4) {
  
  region_gr <- str_to_gr(region)
  
  hits_filt <- hits_gr[overlapsAny(hits_gr, region_gr) ]
  
  # build the annotation table representing the coordinates of the rectangle for each hit
  hits_filt_anno <- tibble(
    xmin  = start(hits_filt),
    xmax  = end(hits_filt),
    name  = hits_filt$name,
    ymin  = -0.2, ymax = 0.2) %>% 
    mutate(x = xmin + (xmax - xmin)/2)
  
  plot <- ggplot(hits_filt_anno) + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = color) +
    geom_text_repel(
      aes(label = name, x = x, y = ymin - 1),
      color = "black", fontface = "bold", size = 3, hjust = 0,
      force = 0.5, min.segment.length = 1
    ) +
    scale_y_continuous(limits = c(-3, 1)) +
    scale_x_continuous(
      limits = c(start(region_gr), end(region_gr)),
      expand = c(0, 0),
      labels = scales::label_comma(big.mark = " ")
    ) +
    labs(x = "Genomic Position (bp)", y = NULL) +
    BPCells:::trackplot_theme() +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none"
    )
  
  # Wrap as BPCells track
  trackplot <- BPCells:::wrap_trackplot(plot, ggplot2::unit(rel_height, "null"), takes_sideplot = FALSE)
  
  if (!is.null(facet_label)) trackplot <- BPCells:::set_trackplot_label(trackplot, labels = facet_label)
  
  return(trackplot)
  
}


#' Make a trackplot annotating motif hits, that's compatible with BPCells plots.
#' 
#' In these plots, each hit in the region is represented by a rectangle, and plot
#' alongside its name.
#'
#' @param hits GRanges object obtained from \code{rtracklayer::import.bed(hits_bed_path)},
#' corresponding to named motif hits.
#' @param region character, genomic coordinates in the form chr:start-end
#' @param track_label character, the y-axis label for the track
#' @param facet_label character, the facet label for the track
#' 
#' @value
#' BPCells trackplot as returned by \code{BPCells:::wrap_trackplot}
trackplot_hits <- function(hits, region, color, track_label = "Hits",
                           facet_label = NULL, rel_height = 0.1,
                           label_size = 4) {

  region_gr <- str_to_gr(region)

  hits_filt <- hits[overlapsAny(hits, region_gr) ]
  
  # build the annotation table representing the coordinates of the rectangle for each hit
  hits_filt_anno <- tibble(
    xmin  = start(hits_filt),
    xmax  = end(hits_filt),
    alpha = hits_filt$score/1000,
    name  = hits_filt$name,
    ymin  = -0.2, ymax = 0.2) %>% 
    mutate(x = xmin + (xmax - xmin)/2)
  
  # plot hits
  plot <- ggplot2::ggplot(hits_filt_anno) +
    ggplot2::geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha = alpha), fill = color) +
    ggplot2::scale_y_continuous(limits = c(-3, 1)) +
    ggrepel::geom_text_repel(aes(label = name, x = x, y = ymin - 1), color = "black", fontface = "bold", size = label_size, hjust = 0,
                             force = 0.5,
                             min.segment.length = 1) +
    ggplot2::scale_alpha(range = c(0.5, 0.8), limits = c(0.7, 1)) +
    BPCells:::trackplot_theme() +
    ggplot2::scale_x_continuous(limits = c(start(region_gr), end(region_gr)),
                                expand = c(0, 0), labels = scales::label_comma(big.mark=" ")) +
    ggplot2::labs(x = "Genomic Position (bp)", y = NULL) +
    ggplot2::theme(axis.ticks.y = element_blank(),
                   axis.text.y = element_blank(),
                   legend.position = "none")
  
  # make this one a bit shorter
  trackplot <- BPCells:::wrap_trackplot(plot, ggplot2::unit(rel_height, "null"), takes_sideplot = FALSE)
  
  if (!is.null(facet_label)) trackplot <- BPCells:::set_trackplot_label(trackplot, labels = facet_label)
  
  return(trackplot)
  
}

#' Adapting trackplot_contribs2 to a second format of contribution scores as BED files as returned by 
#' variant scoring code.
#' 
#' @param gr GRanges object with at least four metadata columns: A, C, G, T which
#' hold the per nucleotide contribution scores per position
trackplot_contribs2 <- function(gr, region,
                                track_label = "Contributions",
                                facet_label = NULL,
                                clip_quantile = 0.999,
                                rel_height = 0.6) {
  
  region_gr <- str_to_gr(region)
  
  if (width(region_gr) > 500) warning("@ not recommended to plot basepair-level contribs for 500 bp")
  else message("@ plotting basepair-level contribs for width ", width(region_gr))
  
  contrib_filt <- gr[gr %over% region_gr]
  
  # rows are nucleotides and cols are positions
  contribs_ohe <- as.matrix(mcols(contrib_filt)[, c("A", "C", "G", "T")]) %>% 
    t()
  
  # make track using ggseqlogo
  plot <- ggseqlogo::ggseqlogo(contribs_ohe, method = "custom", seq_type = "dna") +
    ggplot2::geom_hline(yintercept = 0, color = "gray90") +
    # TODO: right now, the seqlogo has positions 1:N, not start:end, so we can't
    # rescale the x-axis and add the ticks easily. Need to map the tick labels to
    # genomic coordinates.
    # But, still must remove padding! Important for making sure things are aligned
    # to other plots, even if we can't scale the x axis
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = "Genomic Position (bp)", y = track_label) +
    ggplot2::guides(y="none", fill="none") +
    BPCells:::trackplot_theme() +
    
    ggplot2::theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
  
  # make this plot a bit shorter
  trackplot <- BPCells:::wrap_trackplot(plot, ggplot2::unit(rel_height, "null"), takes_sideplot = FALSE)
  
  if (!is.null(facet_label)) trackplot <- BPCells:::set_trackplot_label(trackplot, labels = facet_label)
  
  return(trackplot)
  
  
}


#' Make a trackplot for basepair-level neural network attributions, given
#' an external bw file, that's compatible with BPCells plots.
#' 
#' In these plots, at each position, the nucleotide in the reference genome
#' is plot at a height corresponding to its contribution score. This is implemented
#' using \code{ggseqlogo}.
#' 
#' TODO: we can't get genomic coordinates here because ggseqlogo treats x-axis
#' as positions starting at 1. So we need a way to map 1:N to the actual
#' genomic range for the labels.
#'
#' @param bw GRanges object obtained from \code{rtracklayer::import.bw(bigwig_path)},
#' corresponding to base-resolution contribution scores.
#' @param region character, genomic coordinates in the form chr:start-end
#' @param genome BSGenomes genome object e.g. BSgenome.Hsapiens.UCSC.hg38, used to
#' fetch sequence data.
#' @param track_label character, the y-axis label for the track
#' 
#' @value
#' BPCells trackplot as returned by \code{BPCells:::wrap_trackplot}
trackplot_contribs <- function(bw, region, genome,
                               track_label = "Contributions",
                               facet_label = NULL,
                               clip_quantile = 0.999,
                               rel_height = 0.6,
                               ylim = NULL,          # NEW: override y limits
                               range_label = NULL) { # NEW: optional [min-max] label
  
  # region as GRanges
  region_gr <- str_to_gr(region)
  if (length(region_gr) != 1L) stop("`region` must resolve to a single interval.")
  L <- GenomicRanges::width(region_gr)
  
  if (L > 500) {
    warning("@ not recommended to plot basepair-level contribs for >500 bp")
  } else {
    message("@ plotting basepair-level contribs for width ", L)
  }
  
  # subset GRanges to region
  if (!inherits(bw, "GRanges")) stop("`bw` must be a GRanges with numeric 'score'.")
  contrib_filt <- bw[bw %over% region_gr]
  
  # --------- build per-base score vector aligned to region ----------
  start0 <- GenomicRanges::start(region_gr)
  scores_per_base <- numeric(L)  # default zeros
  
  if (length(contrib_filt) > 0) {
    pos_idx <- as.integer(GenomicRanges::start(contrib_filt)) - start0 + 1L
    keep <- !is.na(pos_idx) & pos_idx >= 1L & pos_idx <= L & is.finite(contrib_filt$score)
    if (any(keep)) {
      pos_idx <- pos_idx[keep]
      vals    <- as.numeric(contrib_filt$score)[keep]
      
      # clip here if clip_quantile is used
      if (is.null(ylim)) {
        # internal clipping based on this track only
        ymax_local <- stats::quantile(vals, clip_quantile, na.rm = TRUE)
        vals <- pmin(vals, ymax_local)
      }
      
      # aggregate if multiple records hit same base (max by default)
      agg <- tapply(vals, pos_idx, max, na.rm = TRUE)
      scores_per_base[as.integer(names(agg))] <- as.numeric(agg)
    }
  }
  
  # --------- override y-limits if provided ----------
  # If ylim is given, clip scores to ylim[2] for consistent across tracks
  if (!is.null(ylim)) {
    scores_per_base <- pmin(scores_per_base, ylim[2])
    scores_per_base <- pmax(scores_per_base, ylim[1])
  }
  
  # extract reference sequence
  region_seq_str <- as.character(Biostrings::getSeq(genome, region_gr))
  bases <- toupper(strsplit(region_seq_str, "", fixed = TRUE)[[1]])
  if (length(bases) != L) stop("@ sequence and region length are not equal.")
  
  # one-hot matrix: positions x {A,C,G,T} with contrib heights
  base_map <- c("A", "C", "G", "T")
  col_idx  <- match(bases, base_map)  # A=1,C=2,G=3,T=4, others NA
  
  mat_ohe <- matrix(0, nrow = L, ncol = 4, dimnames = list(NULL, base_map))
  valid <- !is.na(col_idx)
  if (any(valid)) {
    mat_ohe[cbind(seq_len(L)[valid], col_idx[valid])] <- scores_per_base[valid]
  }
  contribs_ohe <- t(mat_ohe)  # ggseqlogo expects rows=letters, cols=positions
  
  # --------- build ggplot ---------
  plot <- ggseqlogo::ggseqlogo(contribs_ohe, method = "custom", seq_type = "dna") +
    ggplot2::geom_hline(yintercept = 0, color = "gray90") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = "Genomic Position (bp)", y = track_label) +
    ggplot2::guides(y = "none", fill = "none") +
    BPCells:::trackplot_theme() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  
  # attach y-scale if provided
  if (!is.null(ylim)) {
    plot <- plot + ggplot2::scale_y_continuous(limits = ylim, expand = c(0, 0))
  }
  
  # optional label
  if (!is.null(range_label) && !is.null(ylim)) {
    plot <- plot +
      ggplot2::annotate(
        "text",
        x = 1,
        y = ylim[2],
        label = range_label,
        vjust = 1.5,
        hjust = -0.1,
        size = 11 * .8 / ggplot2::.pt
      )
  }
  
  # wrap as track
  trackplot <- BPCells:::wrap_trackplot(
    plot,
    ggplot2::unit(rel_height, "null"),
    takes_sideplot = FALSE
  )
  
  if (!is.null(facet_label)) {
    trackplot <- BPCells:::set_trackplot_label(trackplot, labels = facet_label)
  }
  
  trackplot
}


  
#' Make a trackplot for basepair-level neural network attributions from a BigWig,
#' rendered as a DNA sequence logo with letter heights = contribution scores.
#'
#' Minimal tweak: allow overriding exactly one base (the variant) so that
#' the displayed sequence reflects REF/ALT (or a literal base) while using
#' the provided bw scores as-is.
#'
#' @param bw GRanges from rtracklayer::import.bw(), base-resolution contrib scores.
#' @param region character "chr:start-end"
#' @param genome BSgenome (e.g., BSgenome.Hsapiens.UCSC.hg38)
#' @param track_label character y-axis label
#' @param facet_label optional label for BPCells
#' @param clip_quantile numeric (0–1) to clip tall spikes
#' @param rel_height numeric track relative height
#' @param var_pos integer(1) genomic coordinate of the variant (1-based, hg38 etc.)
#' @param ref,alt optional single letters "A","C","G","T"
#' @param allele_choice character: "ref","alt", or one of "A","C","G","T".
#'        If NULL (default), no override is performed.
trackplot_contribs_allele <- function(bw, region, genome,
                                      track_label = "Contributions",
                                      facet_label = NULL,
                                      clip_quantile = 0.999,
                                      rel_height = 0.6,
                                      var_pos = NULL,
                                      ref = NULL,
                                      alt = NULL,
                                      allele_choice = NULL) {
  
  region_gr <- str_to_gr(region)
  
  if (width(region_gr) > 500) warning("@ not recommended to plot basepair-level contribs for >500 bp")
  else message("@ plotting basepair-level contribs for width ", width(region_gr))
  
  # subset scores for region and order by position to be safe
  contrib_filt <- bw[bw %over% region_gr]
  contrib_filt <- sort(contrib_filt)
  
  # clip extreme values (purely visual)
  ymax <- stats::quantile(contrib_filt$score, clip_quantile)
  ymax_accuracy <- 10^as.integer(log10(0.01 * ymax))
  range_label <- sprintf("[0-%s]", scales::label_comma(accuracy = ymax_accuracy, big.mark=" ")(ymax))
  contrib_filt$score <- pmin(contrib_filt$score, ymax)
  
  # --- fetch reference sequence for the region ---
  dna <- genome[[as.character(seqnames(region_gr))]][start(region_gr):end(region_gr)]
  seq_vec <- strsplit(as.character(dna), "", fixed = TRUE)[[1]]
  stopifnot(length(seq_vec) == width(region_gr))
  L <- length(seq_vec)
  
  # --- optionally override exactly one base (the variant) ---
  if (!is.null(var_pos) && !is.null(allele_choice)) {
    idx <- as.integer(var_pos) - start(region_gr) + 1L
    if (idx >= 1L && idx <= L) {
      allele_choice <- toupper(allele_choice)
      base_to_use <- NULL
      if (allele_choice %in% c("A","C","G","T")) {
        base_to_use <- allele_choice
      } else if (allele_choice == "REF" && !is.null(ref)) {
        base_to_use <- toupper(ref)
      } else if (allele_choice == "ALT" && !is.null(alt)) {
        base_to_use <- toupper(alt)
      }
      if (!is.null(base_to_use) && base_to_use %in% c("A","C","G","T")) {
        seq_vec[idx] <- base_to_use
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # SAFE SCORE ALIGNMENT + MATRIX FILL (this replaces the old mat_ohe[...] <- ...)
  # ---------------------------------------------------------------------------
  
  # Build a score vector for every base in the window.
  # Handles BigWig entries that are >1bp wide and gaps.
  score_vec <- numeric(L)
  if (length(contrib_filt)) {
    # iterate through ranges; assign their (clipped) score across covered indices
    starts <- pmax(1L, start(contrib_filt) - start(region_gr) + 1L)
    ends   <- pmin(L, end(contrib_filt)   - start(region_gr) + 1L)
    scs    <- as.numeric(contrib_filt$score)
    
    for (k in seq_along(scs)) {
      s <- starts[k]; e <- ends[k]
      if (s <= e && s >= 1L && e <= L) {
        # overwrite (typical for single-valued bw); if you prefer sum, use: score_vec[s:e] <- score_vec[s:e] + scs[k]
        score_vec[s:e] <- scs[k]
      }
    }
  }
  
  # Map bases to columns 1..4 (A,C,G,T); non-ACGT (e.g., N) are ignored
  base_levels <- c("A","C","G","T")
  col_idx <- match(toupper(seq_vec), base_levels)  # A=1,C=2,G=3,T=4, others NA
  
  # positions x letters, then transpose for ggseqlogo (letters x positions)
  mat_pos_by_base <- matrix(0, nrow = L, ncol = 4,
                            dimnames = list(NULL, base_levels))
  ok <- !is.na(col_idx)
  if (any(ok)) {
    mat_pos_by_base[cbind(which(ok), col_idx[ok])] <- score_vec[ok]
  }
  contribs_ohe <- t(mat_pos_by_base)
  
  # --- plot ---
  plot <- ggseqlogo::ggseqlogo(contribs_ohe, method = "custom", seq_type = "dna") +
    ggplot2::geom_hline(yintercept = 0, color = "gray90") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = "Genomic Position (bp)", y = track_label) +
    ggplot2::guides(y = "none", fill = "none") +
    BPCells:::trackplot_theme() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  
  trackplot <- BPCells:::wrap_trackplot(plot, ggplot2::unit(rel_height, "null"), takes_sideplot = FALSE)
  if (!is.null(facet_label)) {
    trackplot <- BPCells:::set_trackplot_label(trackplot, labels = facet_label)
  }
  return(trackplot)
}


trackplot_contribs3 <- function(bw, region, genome,
                                                   track_label = "Contributions",
                                                   facet_label = NULL,
                                                   clip_quantile = 0.999,
                                                   rel_height = 0.6) {
  # Requires: GenomicRanges, IRanges, Biostrings, ggseqlogo, ggplot2, scales, BPCells
  stopifnot(requireNamespace("GenomicRanges", quietly = TRUE))
  stopifnot(requireNamespace("IRanges",        quietly = TRUE))
  stopifnot(requireNamespace("Biostrings",     quietly = TRUE))
  stopifnot(requireNamespace("ggseqlogo",      quietly = TRUE))
  stopifnot(requireNamespace("ggplot2",        quietly = TRUE))
  stopifnot(requireNamespace("scales",         quietly = TRUE))
  
  # Region as GRanges (assumes you have str_to_gr() in scope)
  region_gr <- str_to_gr(region)
  if (length(region_gr) != 1L) stop("`region` must resolve to a single interval.")
  L <- GenomicRanges::width(region_gr)
  if (L > 500) {
    warning("@ not recommended to plot basepair-level contribs for >500 bp")
  } else {
    message("@ plotting basepair-level contribs for width ", L)
  }
  
  # Filter contributions to region (FIX: use IRanges::overlapsAny)
  if (!inherits(bw, "GRanges")) stop("`bw` must be a GRanges with a numeric 'score' column.")
  contrib_filt <- bw[IRanges::overlapsAny(bw, region_gr)]
  
  # Robust ymax and axis-label accuracy
  scores_vec <- if (length(contrib_filt) > 0) as.numeric(contrib_filt$score) else numeric(0)
  ymax <- if (length(scores_vec) > 0) {
    unname(stats::quantile(scores_vec, clip_quantile, na.rm = TRUE))
  } else 0
  
  if (!is.finite(ymax) || ymax <= 0) {
    acc <- NULL
  } else {
    acc <- 10^floor(log10(0.01 * ymax))  # ~1% of ymax
    if (!is.finite(acc) || acc <= 0) acc <- NULL
  }
  range_label <- sprintf("[0-%s]", scales::label_number(accuracy = acc, big.mark = " ")(ymax))
  
  # Clip scores at ymax if positive
  if (is.finite(ymax) && ymax > 0 && length(scores_vec) > 0) {
    contrib_filt$score <- pmin(scores_vec, ymax)
  }
  
  # Build per-base score vector aligned to the region
  start0 <- GenomicRanges::start(region_gr)
  scores_per_base <- numeric(L)  # default zeros
  if (length(contrib_filt) > 0) {
    pos_idx <- as.integer(GenomicRanges::start(contrib_filt)) - start0 + 1L
    keep <- !is.na(pos_idx) & pos_idx >= 1L & pos_idx <= L & is.finite(contrib_filt$score)
    if (any(keep)) {
      pos_idx <- pos_idx[keep]
      vals    <- as.numeric(contrib_filt$score)[keep]
      # aggregate if multiple records hit the same base; use max by default
      agg <- tapply(vals, pos_idx, max, na.rm = TRUE)
      scores_per_base[as.integer(names(agg))] <- as.numeric(agg)
    }
  }
  
  # Extract the reference sequence for the region
  region_seq_str <- as.character(Biostrings::getSeq(genome, region_gr))
  bases <- toupper(strsplit(region_seq_str, "", fixed = TRUE)[[1]])
  if (length(bases) != L) stop("@ sequence and region length are not equal.")
  
  # Build one-hot (A/C/G/T) x position matrix with contribution heights
  base_map <- c("A","C","G","T")
  col_idx  <- match(bases, base_map)  # A=1,C=2,G=3,T=4, others (e.g., N)=NA
  
  mat_ohe <- matrix(0, nrow = L, ncol = 4, dimnames = list(NULL, base_map))
  valid <- !is.na(col_idx)
  if (any(valid)) {
    mat_ohe[cbind(seq_len(L)[valid], col_idx[valid])] <- scores_per_base[valid]
  }
  contribs_ohe <- t(mat_ohe)  # ggseqlogo expects rows = letters, cols = positions
  
  # Plot with ggseqlogo
  plot <- ggseqlogo::ggseqlogo(contribs_ohe, method = "custom", seq_type = "dna") +
    ggplot2::geom_hline(yintercept = 0, color = "gray90") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = "Genomic Position (bp)",
                  y = paste0(track_label)) +
    ggplot2::guides(y = "none", fill = "none") +
    BPCells:::trackplot_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
  
  # Wrap as a track
  trackplot <- BPCells:::wrap_trackplot(plot, ggplot2::unit(rel_height, "null"), takes_sideplot = FALSE)
  if (!is.null(facet_label)) {
    trackplot <- BPCells:::set_trackplot_label(trackplot, labels = facet_label)
  }
  trackplot
}


#' Make a paired-allele trackplot from two external bigWig/bedGraph tracks
#'
#' @param bw_ref GRanges or file path to ref-allele signal (bigWig/bedGraph)
#' @param bw_alt GRanges or file path to alt-allele signal (bigWig/bedGraph)
#' @param gene character, gene symbol for region lookup (requires transcripts)
#' @param region character, "chr:start-end"
#' @param track_label character, y-axis label
#' @param facet_label character, facet label (shown at left)
#' @param transcripts transcripts from a BPCells object; required if using `gene`
#' @param flank numeric, bases to extend around `gene` (if `gene` used)
#' @param alpha numeric in [0,1], line alpha
#' @param ymin_zero logical, if TRUE y-axis starts at 0
#' @param rasterize logical, rasterize ggplot layer (requires ggrastr)
#' @param score_cmap named character vector, colors for c("ref","alt")
#' @param return_data logical, return long data.frame instead of plotting
#'
#' @return BPCells trackplot (or data.frame if return_data=TRUE)
trackplot_bw_allelic <- function(bw_ref,
                                 bw_alt,
                                 gene = NULL,
                                 region = NULL,
                                 track_label,
                                 facet_label,
                                 transcripts = NULL,
                                 flank = NULL,
                                 alpha = 0.9,
                                 ymin_zero = TRUE,
                                 rasterize = FALSE,
                                 score_cmap = c("ref" = "blue", "alt" = "red"),
                                 return_data = FALSE,
                                 score_shift = 1) {   
  # --- deps
  stopifnot(!is.null(track_label), !is.null(facet_label))
  suppressPackageStartupMessages({
    library(GenomicRanges); library(IRanges)
    library(dplyr); library(tidyr); library(tibble)
    library(ggplot2); library(glue); library(scales)
  })
  
  # --- small helpers
  .is_path <- function(x) is.character(x) && length(x) == 1 && file.exists(x)
  .import_any <- function(x) {
    if (.is_path(x)) {
      suppressPackageStartupMessages(library(rtracklayer))
      ext <- tools::file_ext(x)
      if (tolower(ext) %in% c("bw","bigwig")) return(rtracklayer::import.bw(x))
      if (tolower(ext) %in% c("bedgraph","bg","bed")) return(rtracklayer::import(x))
      stop("Unsupported file extension: ", ext)
    } else if (inherits(x, "GRanges")) {
      return(x)
    } else stop("Provide GRanges or a valid file path for bigWig/bedGraph.")
  }
  
  message("@ preparing data...")
  gr_ref <- .import_any(bw_ref)
  gr_alt <- .import_any(bw_alt)
  
  # Resolve region
  if (!is.null(gene)) {
    stopifnot(!is.null(transcripts))
    region_gr <- gene_to_gr(transcripts = transcripts, gene = gene, flank = flank)
  } else if (!is.null(region)) {
    region_gr <- str_to_gr(region)
  } else stop("Must provide either 'gene' or 'region'.")
  
  # Keep seqlevels compatible (optional safety)
  seqlevelsStyle(region_gr) <- seqlevelsStyle(gr_ref)[1] %||% seqlevelsStyle(region_gr)[1]
  seqlevelsStyle(gr_alt)    <- seqlevelsStyle(region_gr)[1]
  
  message("@ plotting in region with width ", width(region_gr))
  
  # Subset each track to region
  r_ref <- gr_ref[gr_ref %over% region_gr]
  r_alt <- gr_alt[gr_alt %over% region_gr]
  r_ref$score <- r_ref$score + score_shift
  r_alt$score <- r_alt$score + score_shift
  
  # Score columns: try common names, otherwise "score"
  sco_name <- function(gr) {
    cand <- intersect(c("score","value","signal","coverage"), names(mcols(gr)))
    if (length(cand) == 0) stop("No numeric score column found in GRanges mcols.")
    cand[1]
  }
  s_ref <- sco_name(r_ref)
  s_alt <- sco_name(r_alt)
  
  # Base positions in region (half-open right like bigWig: [start, end))
  positions <- seq.int(start(region_gr), end(region_gr) - 1L)
  
  # Convert to data frame keyed by position
  df_ref <- tibble(pos = start(r_ref), ref = as.numeric(mcols(r_ref)[[s_ref]]))
  df_alt <- tibble(pos = start(r_alt), alt = as.numeric(mcols(r_alt)[[s_alt]]))
  
  # Merge, complete, and arrange
  bw_data <- full_join(df_ref, df_alt, by = "pos") %>%
    complete(pos = positions, fill = list(ref = 0, alt = 0)) %>%
    arrange(pos) %>%
    mutate(facet_label = facet_label)
  
  # Long format
  bw_long <- bw_data %>%
    pivot_longer(cols = c(ref, alt), names_to = "plot_group", values_to = "signal") %>%
    mutate(plot_group = factor(plot_group, levels = c("ref","alt")))
  
  if (return_data) return(bw_long)
  
  # y-range from the actual plotted data
  ymax <- max(bw_long$signal, na.rm = TRUE)
  if (!is.finite(ymax)) ymax <- 0
  ymax_accuracy <- 10^as.integer(log10(pmax(1e-8, 0.01 * (ymax %||% 1))))
  
  if (!ymin_zero) {
    ymin <- min(bw_long$signal, na.rm = TRUE)
    ...
  } else {
    ymin <- 0
    range_label <- sprintf(
      "[0-%s]",
      scales::label_comma(accuracy = ymax_accuracy, big.mark = " ")(ymax)
    )
  }
  
  # Plot
  message("@ plotting...")
  plt <- ggplot(bw_long, aes(group = plot_group)) +
    geom_line(aes(x = pos, y = signal, color = plot_group), alpha = alpha, linewidth=0.7 ) +
    scale_color_manual(values = score_cmap, name = NULL) +
    scale_x_continuous(
      limits = c(start(region_gr), end(region_gr)),
      expand = c(0, 0),
      labels = scales::label_comma(big.mark = " ")
    ) +
    scale_y_continuous(limits = c(ymin, ymax), expand = c(0, 0)) +
    annotate("text",
             x = start(region_gr), y = ymax,
             label = range_label, vjust = 1.5, hjust = -0.1, size = 11 * .8 / .pt) +
    labs(x = "Genomic Position (bp)", y = track_label) +
    guides(y = "none", fill = "none") +
    facet_wrap("facet_label", strip.position = "left") +
    BPCells:::trackplot_theme()
  
  if (!ymin_zero) plt <- plt + geom_hline(yintercept = 0, linewidth = 0.3)
  
  if (isTRUE(rasterize)) {
    suppressPackageStartupMessages(library(ggrastr))
    plt <- ggrastr::rasterize(plt, dpi = 400)
  }
  
  trackplot <- BPCells:::wrap_trackplot(plt, ggplot2::unit(1, "null"), takes_sideplot = FALSE) |>
    BPCells:::set_trackplot_label(labels = facet_label)
  
  return(trackplot)
}


#' Computed binned maximums for a numeric values associated with GRanges in provided bins,
#' similar to GenomicRanges::binnedAverage.
#' 
#' Adapted from 
#' https://divingintogeneticsandgenomics.com/post/compute-averages-sums-on-granges-or-equal-length-bins/
#' and https://stackoverflow.com/a/45683839
#' 
binned_max <- function(bins, numvar, mcolname) {
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
  stopifnot(identical(seqlevels(bins), names(numvar)))
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  sums_list <- lapply(names(numvar),
                      function(seqname) {
                        views <- Views(numvar[[seqname]],
                                       bins_per_chrom[[seqname]])
                        viewMaxs(views)
                      })
  new_mcol <- unsplit(sums_list, as.factor(seqnames(bins)))
  mcols(bins)[[mcolname]] <- new_mcol
  bins
}

#' Wrapper function for \code{binned_max}, to calculate the binned maximums
#' given the region, the bigwig data, and the GRanges metadata column to use for
#' the aggregation.
#' 
#' @param region_gr GRanges
#' @param bw GRanges object obtained from \code{rtracklayer::import.bw(bigwig_path)}
#' @param mcolname character, metadata column name in \code{bw} to aggergate
#' @param tile_width numeric, size of bins (in bp)
#' 
#' @value
#' Returns a new GRanges, where each range is a bin and the metadata contains the
#' aggregated score.
calculate_binned_max <- function(region_gr, bw, mcolname = "score", tile_width) {
  
  bins <- GenomicRanges::tile(region_gr, width = tile_width)[[1]]
  signal_rle <- GenomicRanges::coverage(bw, weight = mcolname)
  seqlevels(bins, pruning.mode="coarse") <- names(signal_rle)
  binned_max <- binned_max(bins, numvar = signal_rle, mcolname = paste0("max_", mcolname))
  
  return(binned_max)
  
}


##### plot causal links instrument->phenotypeA->phenotypeB

plot_region_links_hic <- function(causality_bed,
                                  chr,
                                  x_start,
                                  x_end,
                                  instruments = NULL,   # e.g. c("5:1062047[b38]G,A")
                                  color_map = c(forward = "#56B4E9", reverse = "#E69F00"),
                                  alpha_links = 0.85,
                                  show_rails = TRUE,
                                  draw_peaks_on_A = FALSE,  # draws thin peak segments on baseline
                                  track_label = "Causality",
                                  facet_label = NULL,
                                  rel_height = 0.6,
                                  label_size = 4,
                                  # --- Hi-C arc controls ---
                                  arc_n = 80,
                                  h_min = 0.08,
                                  h_max = 0.90,
                                  scale_bp = 8000,
                                  dist_mode = c("sqrt", "linear", "log"),
                                  baseline = 0,
                                  # --- cosmetics ---
                                  point_color = "red",
                                  point_size = 2.6,
                                  arc_linewidth = 0.7,
                                  pad_frac = 0.03) {
  
  suppressPackageStartupMessages({
    library(dplyr); library(ggplot2); library(scales); library(stringr); library(tidyr)
  })
  dist_mode <- match.arg(dist_mode)
  
  required_cols <- c(
    "chr","pos","exposure_variable",
    "peak_chr","peak_start","peak_end",
    "gene_chr","gene_tss","type"
  )
  stopifnot(all(required_cols %in% names(causality_bed)))
  
  # ---------- helper: build Hi-C arcs (points for geom_path) ----------
  make_arc_df <- function(x1, x2, type, id_prefix, seg) {
    stopifnot(length(x1) == length(x2), length(x1) == length(type))
    if (!length(x1)) return(tibble())
    
    d <- abs(x2 - x1)
    raw_h <- switch(
      dist_mode,
      sqrt   = sqrt(d) / scale_bp,
      linear = d / (scale_bp * 1000),
      log    = log1p(d) / log1p(scale_bp * 1000)
    )
    h <- pmin(h_max, pmax(h_min, h_min + raw_h))
    
    out <- Map(function(a, b, hh, tt, i) {
      t <- seq(0, 1, length.out = arc_n)
      tibble(
        id   = paste0(id_prefix, i),
        seg  = seg,
        type = tt,
        x    = a + (b - a) * t,
        # IMPORTANT: draw arcs DOWN from baseline (Hi-C "hanging" look)
        y    = baseline - hh * sin(pi * t)
      )
    }, x1, x2, h, type, seq_along(x1))
    
    bind_rows(out)
  }
  
  # ---------- prep + region subset ----------
  df0 <- causality_bed %>%
    as_tibble() %>%
    mutate(
      type       = tolower(type),
      type       = ifelse(type %in% c("forward","reverse"), type, NA_character_),
      pos        = as.numeric(pos),
      peak_start = as.numeric(peak_start),
      peak_end   = as.numeric(peak_end),
      gene_tss   = as.numeric(gene_tss),
      peak_mid   = (peak_start + peak_end) / 2
    ) %>%
    filter(chr == !!chr, peak_chr == !!chr, gene_chr == !!chr) %>%
    filter(
      between(pos, x_start, x_end) |
        (peak_end >= x_start & peak_start <= x_end) |
        between(gene_tss, x_start, x_end)
    )
  
  # subset instruments (optional)
  if (!is.null(instruments)) {
    df0 <- df0 %>% filter(exposure_variable %in% instruments)
    if (nrow(df0) == 0) stop("No rows for requested instrument(s) in region on ", chr, ".")
  }
  
  # phenotype A/B targets by direction
  links <- df0 %>%
    mutate(
      x_A_raw = if_else(type == "forward", peak_mid, gene_tss),  # phenotype A
      x_B_raw = if_else(type == "forward", gene_tss,  peak_mid)  # phenotype B
    ) %>%
    mutate(
      pos_clamp = pmin(pmax(pos,     x_start), x_end),
      x_A       = pmin(pmax(x_A_raw, x_start), x_end),
      x_B       = pmin(pmax(x_B_raw, x_start), x_end)
    ) %>%
    filter(is.finite(pos_clamp), is.finite(x_A), is.finite(x_B), !is.na(type))
  
  # unique variants (for points)
  variants_unique <- df0 %>%
    distinct(exposure_variable, chr, pos) %>%
    mutate(pos_clamp = pmin(pmax(pos, x_start), x_end)) %>%
    filter(is.finite(pos_clamp))
  
  # optional: peak segments rendered on baseline (colored by type)
  peaks_on_A <- df0 %>%
    distinct(peak_chr, peak_start, peak_end, type) %>%
    transmute(
      type,
      x1 = pmin(pmax(peak_start, x_start), x_end),
      x2 = pmin(pmax(peak_end,   x_start), x_end)
    ) %>%
    filter(is.finite(x1), is.finite(x2), x2 >= x1)
  
  # ---------- build arc point data ----------
  arc_vA <- make_arc_df(links$pos_clamp, links$x_A, links$type, "vA_", "v_to_A")
  arc_AB <- make_arc_df(links$x_A,       links$x_B, links$type, "AB_", "A_to_B")
  
  # ---------- tight y-lims (baseline at top, arcs below) ----------
  y_vals <- c(baseline, arc_vA$y, arc_AB$y)
  y_vals <- y_vals[is.finite(y_vals)]
  if (length(y_vals) == 0) y_vals <- baseline
  
  pad <- pad_frac * max(1e-6, diff(range(y_vals)))
  y_low  <- min(y_vals) - pad
  y_high <- baseline + 0  # no extra whitespace above baseline
  
  # ---------- plot (Hi-C style) ----------
  p <- ggplot() +
    { if (show_rails) geom_hline(yintercept = baseline, linewidth = 0.2, color = "grey65") else NULL } +
    
    { if (draw_peaks_on_A && nrow(peaks_on_A) > 0)
      geom_segment(data = peaks_on_A,
                   aes(x = x1, xend = x2, y = baseline, yend = baseline, color = type),
                   linewidth = 1, alpha = 0.35)
      else NULL } +
    
    geom_point(data = variants_unique,
               aes(x = pos_clamp, y = baseline),
               color = point_color, size = point_size) +
    
    { if (nrow(arc_vA) > 0)
      geom_path(data = arc_vA,
                aes(x = x, y = y, group = id, color = type),
                linewidth = arc_linewidth, alpha = alpha_links)
      else NULL } +
    
    { if (nrow(arc_AB) > 0)
      geom_path(data = arc_AB,
                aes(x = x, y = y, group = id, color = type),
                linewidth = arc_linewidth, alpha = alpha_links, linetype = "22")
      else NULL } +
    
    scale_x_continuous(
      limits   = c(x_start, x_end),
      labels   = scales::label_comma(accuracy = 1),
      expand   = c(0.01, 0.01),
      position = "top"
    ) +
    # IMPORTANT: reverse y so baseline sits at top and arcs "hang" downward
    scale_y_reverse(
      limits = c(y_high, y_low),
      breaks = NULL,
      expand = c(0, 0)
    ) +
    scale_color_manual(values = color_map, breaks = names(color_map), limits = names(color_map)) +
    labs(x = paste0(chr, " position (bp)"),
         y = NULL,
         color = "Direction") +
    theme_minimal(base_size = 12) +
    BPCells:::trackplot_theme() +
    theme(
      axis.title.x = element_text(vjust = 0.5),
      axis.text.x  = element_text(vjust = 0.5)
    )
  
  # ---------- optional BPCells wrapping ----------
  if (requireNamespace("BPCells", quietly = TRUE) &&
      exists("wrap_trackplot", where = asNamespace("BPCells"), inherits = FALSE)) {
    p <- BPCells:::wrap_trackplot(p, grid::unit(rel_height, "null"), takes_sideplot = FALSE)
    if (!is.null(facet_label) &&
        exists("set_trackplot_label", where = asNamespace("BPCells"), inherits = FALSE)) {
      p <- BPCells:::set_trackplot_label(p, labels = facet_label)
    }
  }
  
  return(p)
}




gtf_tbl_to_exon_gr <- function(gtf_tbl) {
  ex <- gtf_tbl %>%
    filter(feature == "exon") %>%
    mutate(
      start  = as.integer(start),
      end    = as.integer(end),
      strand = ifelse(strand %in% c("+","-"), strand, "*")
    )
  gr <- makeGRangesFromDataFrame(
    ex,
    seqnames.field     = "chr",
    start.field        = "start",
    end.field          = "end",
    strand.field       = "strand",
    keep.extra.columns = TRUE
  )
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  seqlevelsStyle(gr) <- "UCSC"   # ensures 'chr7' style
  gr
}


################# ############################################
################ PARSE INFO FROM VARIANT NAME ###############
#############################################################
# ---- strict parser (only 2 formats) ------------------------------------------
# Accepts ONLY:
#   - "chr6_113306741_G_A"
#   - "6:113306741[b38]G,A"
parse_variant_id <- function(variant_id) {
  x <- trimws(variant_id)
  
  # A) chr6_113306741_G_A
  mA <- stringr::str_match(x, "^chr([0-9XYMT]+)_([0-9]+)_[ACGT]_[ACGT]$")
  if (!is.na(mA[1, 1])) {
    return(list(chr = paste0("chr", mA[1, 2]), pos = as.integer(mA[1, 3])))
  }
  
  # B) 6:113306741[b38]G,A
  mB <- stringr::str_match(x, "^([0-9XYMT]+):([0-9]+)\\[b\\d+\\][ACGT],[ACGT]$")
  if (!is.na(mB[1, 1])) {
    return(list(chr = paste0("chr", mB[1, 1 + 0]), pos = as.integer(mB[1, 2 + 0])))
  }
  
  stop("Unrecognized variant_id format. Allowed only: ",
       '"chr6_113306741_G_A" or "6:113306741[b38]G,A". Got: ', variant_id)
}

#######################################################################################
############################ Genes, for gene track #############################
######################################################################################
# ---- Gene GRanges for a region ---------------------------------------------
get_subset_genes_for_region <- function(
    region,                                    # "chr6:113306661-113306821" OR list(chr=,start=,end=)
    txdb   = TxDb.Hsapiens.UCSC.hg38.knownGene,
    orgdb  = org.Hs.eg.db,
    keep_std = TRUE,
    style  = "UCSC"                            # UCSC to match "chr*" style regions
) {
  # deps
  require(GenomicFeatures)
  require(GenomeInfoDb)
  require(GenomicRanges)
  require(AnnotationDbi)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  require(org.Hs.eg.db)
  require(stringr)
  
  # parse region
  if (is.character(region)) {
    m <- stringr::str_match(region, "^([A-Za-z0-9._]+):([0-9,]+)-([0-9,]+)$")
    if (is.na(m[1,1])) stop("Region string must be like 'chr6:113306661-113306821'")
    r_chr   <- m[1,2]
    r_start <- as.integer(gsub(",", "", m[1,3]))
    r_end   <- as.integer(gsub(",", "", m[1,4]))
  } else if (is.list(region) && all(c("chr","start","end") %in% names(region))) {
    r_chr <- region$chr; r_start <- as.integer(region$start); r_end <- as.integer(region$end)
  } else {
    stop("`region` must be a 'chr:start-end' string or a list(chr=, start=, end=).")
  }
  
  # 1) one GRanges per gene (not transcripts)
  g <- genes(txdb, columns = "gene_id")  # may drop multi-seq/strand genes
  seqlevelsStyle(g) <- style
  if (keep_std) g <- keepStandardChromosomes(g, pruning.mode = "coarse")
  
  # 2) Entrez -> SYMBOL
  sym <- AnnotationDbi::mapIds(
    orgdb,
    keys      = as.character(mcols(g)$gene_id),
    keytype   = "ENTREZID",
    column    = "SYMBOL",
    multiVals = "first"
  )
  
  # 3) columns expected by trackplot_gene()
  mcols(g)$feature       <- "gene"  # initial tag (we'll switch to 'transcript' after subsetting)
  mcols(g)$gene_name     <- unname(sym[as.character(mcols(g)$gene_id)])
  mcols(g)$transcript_id <- mcols(g)$gene_name
  
  # 4) subset to region
  region_gr <- GRanges(seqnames = r_chr, ranges = IRanges(start = r_start, end = r_end))
  subset_g  <- subsetByOverlaps(g, region_gr, type = "any")
  
  # 5) tweak labels/feature as you showed
  if (length(subset_g)) {
    mcols(subset_g)$feature       <- "transcript"
    mcols(subset_g)$transcript_id <- as.character(mcols(subset_g)$gene_name)
    if (keep_std) subset_g <- keepStandardChromosomes(subset_g, pruning.mode = "coarse")
  }
  
  subset_g
}

trackplot_bw_paired_2bw <- function(bw_ref,
                                    bw_alt,
                                    gene = NULL,
                                    region = NULL,
                                    track_label,
                                    facet_label,
                                    transcripts = NULL,
                                    flank = NULL,
                                    alpha = 0.6,
                                    ymin_zero = TRUE,
                                    rasterize = FALSE,
                                    score_cmap = c("ref" = "blue", "alt" = "red"),
                                    return_data = FALSE) {
  
  message("@ preparing data...")
  
  # get region GRanges
  if (!is.null(gene)) {
    region_gr <- gene_to_gr(transcripts = transcripts, gene = gene, flank = flank)
  } else if (!is.null(region)) {
    region_gr <- str_to_gr(region)
  } else if (is.null(gene) & is.null(region)) {
    stop("Must provide either gene or region")
  }
  
  message("@ plotting in region with width ", width(region_gr))
  
  # subset each bigwig to the region
  bw_ref_filt <- bw_ref[bw_ref %over% region_gr]
  bw_alt_filt <- bw_alt[bw_alt %over% region_gr]
  
  # very light sanity check: same ranges & seqnames
  if (!identical(seqnames(bw_ref_filt), seqnames(bw_alt_filt)) ||
      !identical(ranges(bw_ref_filt), ranges(bw_alt_filt))) {
    stop("bw_ref and bw_alt must have identical seqnames and ranges in the plotted region.")
  }
  
  # merge into a single GRanges with score1 / score2, like the original function expects
  bw_filt <- bw_ref_filt
  bw_filt$score1 <- bw_ref_filt$score
  bw_filt$score2 <- bw_alt_filt$score
  
  message("@ plotting without binning.")
  
  # basepair positions in region
  positions <- seq(start(region_gr), end(region_gr) - 1)
  
  # wide format (same as original, but score1/score2 come from ref/alt)
  bw_data <- tibble::tibble(
    pos         = start(bw_filt),
    signal1     = bw_filt$score1,
    signal2     = bw_filt$score2,
    facet_label = facet_label
  ) %>%
    tidyr::complete(
      pos = positions,
      fill = list(
        signal1     = 0,
        signal2     = 0,
        facet_label = facet_label
      )
    ) %>%
    dplyr::arrange(pos)
  
  # transform to long format (no set_colnames; minimal change)
  bw_data <- bw_data %>%
    tidyr::pivot_longer(
      cols      = c(signal1, signal2),
      names_to  = "plot_group",
      values_to = "signal"
    )
  
  # rename groups -> use names(score_cmap)
  bw_data$plot_group[bw_data$plot_group == "signal1"] <- names(score_cmap)[1]
  bw_data$plot_group[bw_data$plot_group == "signal2"] <- names(score_cmap)[2]
  
  # optional check (same as original)
  nrow(bw_data) == length(positions)
  
  if (return_data) return(bw_data)
  
  # y-range from original scores (before completion)
  ymax <- max(bw_filt$score1, bw_filt$score2)
  ymax_accuracy <- 10^as.integer(log10(0.01 * ymax))
  
  if (!ymin_zero) {
    ymin <- min(bw_filt$score1, bw_filt$score2)
    ymin_accuracy <- 10^as.integer(log10(0.01 * abs(ymin)))
    range_label <- glue::glue(
      "[{scales::label_comma(accuracy = ymin_accuracy, big.mark=' ')(ymin)}-",
      "{scales::label_comma(accuracy = ymax_accuracy, big.mark=' ')(ymax)}]"
    )
  } else {
    ymin <- 0
    range_label <- sprintf(
      "[0-%s]",
      scales::label_comma(accuracy = ymax_accuracy, big.mark = " ")(ymax)
    )
  }
  
  # plot track (identical style)
  message("@ plotting...")
  plot <- ggplot2::ggplot(bw_data, ggplot2::aes(group = plot_group)) +
    ggplot2::geom_line(
      ggplot2::aes(x = pos, y = signal, color = plot_group),
      alpha = alpha, linewidth=1
    ) +
    ggplot2::scale_color_manual(values = score_cmap) +
    ggplot2::scale_x_continuous(
      limits = c(start(region_gr), end(region_gr)),
      expand = c(0, 0),
      labels = scales::label_comma(big.mark = " ")
    ) +
    ggplot2::scale_y_continuous(
      limits = c(ymin, ymax),
      expand = c(0, 0)
    ) +
    ggplot2::annotate(
      "text",
      x     = start(region_gr),
      y     = ymax,
      label = range_label,
      vjust = 1.5,
      hjust = -0.1,
      size  = 11 * .8 / ggplot2::.pt
    ) +
    ggplot2::labs(x = "Genomic Position (bp)", y = track_label) +
    ggplot2::guides(y = "none", fill = "none") +
    ggplot2::facet_wrap("facet_label", strip.position = "left") +
    BPCells:::trackplot_theme()
  
  if (!ymin_zero) {
    plot <- plot + ggplot2::geom_hline(yintercept = 0)
  }
  
  if (rasterize) {
    suppressPackageStartupMessages(library(ggrastr))
    # original code used layer = plot_as (which doesn't exist); rasterize whole plot instead
    plot <- ggrastr::rasterize(plot, dpi = 400)
  }
  
  trackplot <- BPCells:::wrap_trackplot(
    plot,
    ggplot2::unit(1, "null"),
    takes_sideplot = FALSE
  ) %>%
    BPCells:::set_trackplot_label(labels = facet_label)
  
  return(trackplot)
}


#' Make a trackplot from an external bw or bedGraph file, that's compatible with BPCells plots
#'
#' This function calculates smoothed data (taking the maximum value within bins)
#' or basepair-resolution data and plots it along the genome.
#'
#' @param bw GRanges object obtained from \code{rtracklayer::import.bw(bigwig_path)} or
#' \code{rtracklayer::import.bedGraph(bg_path)}
#' @param gene character, gene symbol corresponding to the genomic region which should be plot 
#' @param region character, genomic coordinates in the form chr:start-end
#' @param track_label character, the y-axis label for the track
#' @param facet_label character, the facet label for the plot
#' @param transcripts transcripts from a BPCells object. Required if plotting based
#' on a gene.
#' @param flank numeric, specifying how much to extend on either side of the gene for plotting 
#' @param tile_width numeric, specifis bin size for aggregating data. Signal within
#' bins will be aggregated by taking the max in each bin
#' @param plot_as character, one of "area" or "bar", controlling whether signal should
#' be plot as an area or bar plot. Bar plot recommended for smaller regions.
#' @param clip_quantile numeric, quantile for clipping max values. Default: 0.999
#' @param color character, track color
#' @param return_data logical, whether to return input data before plotting.
#' 
#' @value
#' BPCells trackplot as returned by \code{BPCells:::wrap_trackplot}
trackplot_bw <- function(bw,
                         gene = NULL,
                         region = NULL,
                         track_label,
                         facet_label,
                         transcripts = NULL,
                         flank = NULL,
                         tile_width = 100,
                         plot_as = "area",
                         clip_quantile = 0.999,
                         ymin_zero = TRUE,
                         rasterize = TRUE,
                         return_data = FALSE,
                         color = "black") {
  
  message("@ preparing data...")
  
  # get a GRanges containing the region of interest, either from gene or region string
  if (!is.null(gene)) region_gr <- gene_to_gr(transcripts = transcripts, gene = gene, flank = flank)
  else if (!is.null(region)) {
    
    region_gr <- str_to_gr(region)
    
  } else if (is.null(gene) & is.null(region)) stop("Must provide either gene or region")
  
  message("@ plotting in region with width ", width(region_gr))
  
  # subset the bigwig
  bw_filt <- bw[bw %over% region_gr]
  
  # bw_filt <- BRGenomics::makeGRangesBRG(bw_filt)
  
  if (tile_width > 1) {
    
    # get binned values
    message("@ binning data with tile width ", tile_width)
    data_binned_max <- calculate_binned_max(region_gr, bw_filt, mcolname = "score", tile_width = tile_width)
    
    # calculate bin centers follwing BPCells::trackplot_coverage
    bin_centers <- seq(start(region_gr), end(region_gr) - 1, tile_width) + (tile_width/2)
    bin_centers <- pmin(bin_centers, end(region_gr) - 1)
    
    bw_data <- tibble::tibble(
      pos = bin_centers,
      signal = data_binned_max$max_score,
      facet_label = facet_label
    )
    
  } else {
    
    message("@ plotting without binning.")
    
    # get the coordinates of the basepair positions in the region
    positions <- seq(start(region_gr), end(region_gr)-1)
    
    # convert the data over the specified ranges to a dataframe
    bw_data <- tibble::tibble(pos = start(bw_filt), signal = bw_filt$score, facet_label = facet_label) %>%
      # complete missing positions for the regions with 0
      complete(pos = positions,
               fill = list("signal" = 0,
                           "facet_label" = facet_label)) %>%
      arrange(pos)
    
    # double check that lengths are the same
    nrow(bw_data) == length(positions)
    
  }
  
  if (return_data) return(bw_data)
  
  # get clipped ymax: note that we do this on the data in the bigwig *prior* to
  # completing the positions with 0s
  ymax <- quantile(bw_filt$score, clip_quantile)
  ymax_accuracy <- 10^as.integer(log10(0.01 * ymax))
  
  if (!ymin_zero) {
    
    ymin <- quantile(bw_filt$score, 1-clip_quantile)
    ymin_accuracy <- 10^as.integer(log10(0.01 * ymin))
    
  } else {
    
    ymin = 0
    
  }
  
  range_label <- sprintf("[0-%s]", scales::label_comma(accuracy = ymax_accuracy, big.mark=" ")(ymax))
  
  # clip values if needed
  bw_data$signal <- pmin(bw_data$signal, ymax)
  
  # plot track
  message("@ plotting...")
  plot <- ggplot2::ggplot(bw_data)
  
  if (plot_as == "area") {
    
    plot <- plot +
      ggplot2::geom_area(ggplot2::aes(x = pos, y = signal), fill = color)
    
  } else if (plot_as == "bar") {
    
    plot <- plot +
      ggplot2::geom_bar(ggplot2::aes(x = pos, y = signal), fill = color,
                        # don't leave spaces between bars
                        width = 1.1,
                        stat = "identity")
    
  }
  
  plot <- plot +
    ggplot2::scale_x_continuous(limits = c(start(region_gr), end(region_gr)), expand = c(0, 0), labels = scales::label_comma(big.mark=" ")) +
    ggplot2::scale_y_continuous(limits = c(ymin, ymax), expand = c(0, 0)) +
    ggplot2::annotate("text", x = start(region_gr), y = ymax, label = range_label, vjust=1.5, hjust=-0.1, size=11*.8/ggplot2::.pt) +
    ggplot2::labs(x = "Genomic Position (bp)", y = track_label) +
    ggplot2::guides(y="none", fill="none") +
    # facetting is used to get a strip label on the left side, even if only one track plotted
    ggplot2::facet_wrap("facet_label", strip.position="left") +
    BPCells:::trackplot_theme()
  
  # plot track
  message("@ plotting...")
  
  # 1) build the heavy layer (optionally rasterised)
  if (plot_as == "area") {
    lyr <- ggplot2::geom_area(ggplot2::aes(x = pos, y = signal), fill = color)
  } else if (plot_as == "bar") {
    lyr <- ggplot2::geom_bar(
      ggplot2::aes(x = pos, y = signal),
      fill = color,
      width = 1.1,
      stat = "identity"
    )
  } else {
    stop("plot_as must be 'area' or 'bar'")
  }
  
  if (rasterize) lyr <- ggrastr::rasterise(lyr, dpi = 400)
  
  # 2) assemble the full plot with theme/scales/etc
  plot <- ggplot2::ggplot(bw_data) +
    lyr +
    ggplot2::scale_x_continuous(
      limits = c(start(region_gr), end(region_gr)),
      expand = c(0, 0),
      labels = scales::label_comma(big.mark = " ")
    ) +
    ggplot2::scale_y_continuous(limits = c(ymin, ymax), expand = c(0, 0)) +
    ggplot2::annotate(
      "text",
      x = start(region_gr), y = ymax,
      label = range_label,
      vjust = 1.5, hjust = -0.1,
      size = 11 * .8 / ggplot2::.pt
    ) +
    ggplot2::labs(x = "Genomic Position (bp)", y = track_label) +
    ggplot2::guides(y = "none", fill = "none") +
    ggplot2::facet_wrap("facet_label", strip.position = "left") +
    BPCells:::trackplot_theme()
  
  if (plot_as == "area") {
    lyr <- ggplot2::geom_area(ggplot2::aes(x = pos, y = signal), fill = color)
  } else if (plot_as == "bar") {
    lyr <- ggplot2::geom_bar(
      ggplot2::aes(x = pos, y = signal),
      fill = color,
      width = 1.1,
      stat = "identity"
    )
  }

  trackplot <- BPCells:::wrap_trackplot(plot, ggplot2::unit(1, "null"), takes_sideplot = FALSE) %>%
    BPCells:::set_trackplot_label(labels = facet_label)
  
  return(trackplot)
  
}


trackplot_gene_custom <- function(
    transcripts,
    region,
    exon_size = 2.5,
    gene_size = 0.5,
    label_size = 11 * 0.8 / ggplot2::.pt,
    track_label = "Genes",
    return_data = FALSE
) {
  # BPCells internals
  region <- BPCells:::normalize_ranges(region)
  
  transcripts <- BPCells:::normalize_ranges(
    transcripts,
    metadata_cols = c("strand", "feature", "gene_id", "gene_name")
  )
  
  size_range <- c(min(exon_size, gene_size, label_size),
                  max(exon_size, gene_size, label_size))
  linewidth_range <- size_range / 0.75
  exon_size <- exon_size / 0.75
  gene_size <- gene_size / 0.75
  
  data <- transcripts %>%
    tibble::as_tibble() %>%
    dplyr::filter(
      as.character(chr) == as.character(region$chr),
      end > region$start,
      start < region$end,
      feature %in% c("gene", "transcript", "exon")
    ) %>%
    dplyr::mutate(
      strand = dplyr::case_when(
        strand %in% c("+", "plus", "TRUE", TRUE) ~ TRUE,
        strand %in% c("-", "minus", "FALSE", FALSE) ~ FALSE,
        TRUE ~ NA
      ),
      gene_key = dplyr::coalesce(gene_name, gene_id),
      size = dplyr::if_else(feature == "exon", exon_size, gene_size)
    )
  
  if (nrow(data) == 0) {
    if (return_data) {
      data$y <- numeric(0)
      data$facet_label <- character(0)
      return(list(data = data, arrows = data))
    } else {
      return(BPCells:::trackplot_empty(region, label = track_label))
    }
  }
  
  feature_main <- if (any(data$feature == "gene")) "gene" else "transcript"
  
  main_data <- dplyr::filter(data, feature == feature_main) %>%
    dplyr::mutate(y = BPCells:::trackplot_calculate_segment_height(.))
  
  if (anyDuplicated(main_data$gene_key)) {
    dup <- main_data$gene_key[anyDuplicated(main_data$gene_key)][1]
    rlang::abort(sprintf(
      "Found multiple feature == \"%s\" rows with same gene_key (e.g. %s). Ensure 1 gene row per gene.",
      feature_main, dup
    ))
  }
  
  data <- dplyr::left_join(
    data,
    dplyr::select(main_data, gene_key, y),
    by = "gene_key"
  )
  
  if (anyNA(data$y)) {
    rlang::abort(sprintf(
      "Found rows with gene_key missing a feature == \"%s\" row: %s",
      feature_main, data$gene_key[which(is.na(data$y))[1]]
    ))
  }
  
  main_coords <- dplyr::filter(data, feature == feature_main)
  arrows <- BPCells:::trackplot_create_arrow_segs(main_coords, region, 50)
  
  data <- dplyr::mutate(
    data,
    start = pmax(region$start, pmin(region$end, start)),
    end   = pmax(region$start, pmin(region$end, end)),
    facet_label = track_label
  )
  
  if (return_data) return(list(data = data, arrows = arrows))
  
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x    = dplyr::if_else(strand, start, end),
      xend = dplyr::if_else(strand, end, start),
      y = y, yend = y,
      linewidth = size
    )
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(color = factor(strand, levels = c(TRUE, FALSE), labels = c("+", "-"))),
      show.legend = TRUE
    ) +
    ggplot2::geom_segment(
      data = arrows,
      ggplot2::aes(color = factor(strand, levels = c(TRUE, FALSE), labels = c("+", "-"))),
      arrow = grid::arrow(length = grid::unit(0.4 * exon_size, "mm")),
      show.legend = TRUE
    ) +
    ggrepel::geom_text_repel(
      data = dplyr::filter(data, feature == feature_main),
      ggplot2::aes(label = gene_name),
      size = label_size,
      position = ggrepel::position_nudge_repel(y = 0.25)
    ) +
    ggplot2::scale_size(range = size_range, limits = size_range) +
    ggplot2::scale_linewidth(range = linewidth_range, limits = linewidth_range) +
    ggplot2::scale_color_manual(values = c(`+` = "black", `-` = "darkgrey"), drop = FALSE) +
    ggplot2::scale_x_continuous(
      limits = c(region$start, region$end),
      expand = c(0, 0),
      labels = scales::label_number()
    ) +
    ggplot2::scale_y_discrete(labels = NULL, breaks = NULL) +
    ggplot2::labs(x = "Genomic Position (bp)", y = NULL, color = "strand") +
    ggplot2::guides(size = "none", linewidth = "none") +
    ggplot2::facet_wrap("facet_label", strip.position = "left") +
    BPCells:::trackplot_theme()
  
  BPCells:::wrap_trackplot(plot, height = ggplot2::unit(1, "null"), region = region)
}


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


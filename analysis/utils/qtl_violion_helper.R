library(genio)
library(dplyr)
library(ggplot2)
library(stringr)
library(viridis)

plot_eqtl_violin_from_string <- function(
    gene,
    variant_string,   # e.g. "17:39946957[b38]G,A"
    celltype = "allcells",
    expression_bed_dir = "~/cd4_QTL_analysis/02_Gene_expression/analysis/004_ProcessExpression",
    plink_path = "~/cd4_QTL_analysis/01_genotype_snps_covar/02_genotype_covariates/analysis/002.v4_calculating_pcs_MAF5/per_chr_plink_files/",
    out_dir = "~/cd4_qtl_paper_figures/figure_2/plotting/plots",
    point_size = 0.3
) {
  # -----------------------
  # 1) Parse chromosome
  # -----------------------
  info <- unlist(stringr::str_split(variant_string, ":"))
  chr  <- info[[1]]
  
  # -----------------------
  # 2) Build paths
  # -----------------------
  expression_bed <- file.path(
    expression_bed_dir,
    paste0(celltype, "_pseudo_cells_mean_mx.bed")
  )
  
  bim_file <- paste0(plink_path, "chr", chr, "_ashkenazi.367.AF1.QC.BA.king2.hwe.annot.bim")
  fam_file <- paste0(plink_path, "chr", chr, "_ashkenazi.367.AF1.QC.BA.king2.hwe.annot.fam")
  bed_file <- paste0(plink_path, "chr", chr, "_ashkenazi.367.AF1.QC.BA.king2.hwe.annot.bed")
  
  # -----------------------
  # 3) Read PLINK geno
  # -----------------------
  bim <- read_bim(bim_file)
  fam <- read_fam(fam_file)
  gen <- read_bed(
    bed_file,
    bim$id,
    gsub("-", "_", fam$id)
  )
  
  variant <- as.character(variant_string)
  
  if (!variant %in% rownames(gen)) {
    stop("Variant ", variant, " not found in genotype matrix.\n",
         "First few rownames:\n",
         paste(head(rownames(gen)), collapse = ", "))
  }
  
  # Same pattern that works for you
  snp <- as.data.frame(gen[variant, ])
  colnames(snp) <- c(variant)
  snp$ind <- rownames(snp)
  colnames(snp) <- c("variant", "ind")
  
  # -----------------------
  # 4) Read expression
  # -----------------------
  expr <- read.delim(expression_bed)
  colnames(expr) <- gsub("^X", "", colnames(expr))
  
  if (!"gene_name" %in% colnames(expr)) {
    stop("Expression BED must have a 'gene_name' column.")
  }
  
  expr_row <- expr[expr$gene_name == gene, ]
  if (nrow(expr_row) == 0) {
    stop("Gene ", gene, " not found in expression matrix.")
  }
  
  expr_g <- as.data.frame(t(expr_row[, 5:ncol(expr_row)]))
  colnames(expr_g) <- "expr"
  expr_g$ind <- gsub("^X", "", rownames(expr_g))
  
  # -----------------------
  # 5) Merge & summarize
  # -----------------------
  all <- merge(snp, expr_g, by = "ind")
  all <- na.omit(all)
  
  # variant as factor (0,1,2) and numeric dosage
  all$variant   <- as.factor(all$variant)
  all$geno_num  <- as.numeric(as.character(all$variant))  # assumes 0/1/2
  
  Means <- all %>%
    group_by(variant, geno_num) %>%
    summarise(Avg = median(expr), .groups = "drop")
  
  # Fit regression line using the three medians
  fit <- lm(Avg ~ geno_num, data = Means)
  pred_df <- data.frame(
    geno_num = range(Means$geno_num)
  )
  pred_df$fit <- predict(fit, newdata = pred_df)
  
  # -----------------------
  # 6) Plot
  # -----------------------
  p <- ggplot(all, aes(x = geno_num, y = expr)) +
    geom_violin(
      aes(group = variant, fill = variant),
      alpha = 0.7
    ) +
    geom_point(
      size = point_size,
      position = position_jitter(width = 0.1, height = 0)
    ) +
    # regression line through the group medians
    geom_line(
      data = pred_df,
      aes(x = geno_num, y = fit),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "black"
    ) +
    scale_x_continuous(
      breaks = sort(unique(Means$geno_num)),
      labels = sort(as.character(unique(Means$variant)))
    ) +
    ggtitle(paste0(variant, " in ", celltype)) +
    ylab(paste0("Residualized Expression ",gene)) +
    xlab("Genotype") +
    scale_fill_viridis(discrete = TRUE) +
    theme_classic()
  
  invisible(p)
}


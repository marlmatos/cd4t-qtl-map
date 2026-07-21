rm(list = ls())
options(bitmapType="cairo")
library(data.table)
library(dplyr)
library(cowplot)
library(geni.plots)

setwd("/gchm_lab/collab/marlis_pj/coloc")

file_name = "coloc_results/eqtl_gwas_coloc/2021_34594039_IBD_EUR-EAS_preprocessed_coloc_results.csv"
gwas_file = gsub("_coloc_results","",gsub(".csv","",basename(file_name)))

# Load colocalization results
coloc_results = fread(file_name) %>%
  filter(pval_nominal < 2e-04 & pval < 5e-08)

make_locus_plot = function(x){
  
  print("making plot")
  
  eqtl_variant = coloc_results$variant_id_GRC38[x]
  eqtl_name = "All_CD4T_cells"
  disease_name = str_extract(file_name, "[^_/]+_[^_/]+(?=_preprocessed)")
  region = str_split_fixed(coloc_results$region_GRC38[x], ":",2)[,2]
  chr = strsplit(eqtl_variant,":", 4)[[1]][1]
  gene_name = coloc_results$gene[x]
  
  # Load GWAS sumstats
  gwas_stats = fread(paste0('/gchm_lab/collab_gwas/finemapping_autoimmune/data/preprocessed_v1/',gwas_file,'.tsv')) %>%
    mutate(variant_id_GRC37 = paste(chromosome, position,sep=":"))
  
  # Load eqtl
  eqtl_lbf = fread(paste0("SuSiE_finemap_lbf_GRC37/",eqtl_name,"/",eqtl_name,"_chr",chr,"_lbfs_GRC37.txt")) %>%
    mutate(variant_id_GRC37 = paste(gsub("chr","",chr), pos_GRC37, sep = ":")) %>%
    filter(gene == gene_name)
  
  # Load LD matrix
  ld_dir = paste0('Ashkenazi_LD_matrices/chr',chr,'/',region,'/',region,'_imputed')
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
  
  LD = LD_clean %>% select(-variant_id_GRC38)
  
  markers <- markers[match(colnames(LD), markers$marker), ]
  
  identical(colnames(LD), markers$marker)
  
  stack_plot <- fig_region_stack(data = markers, 
                                 trait = c(disease_name, gene_name), 
                                 corr = as.matrix(LD), 
                                 build =38, title_center = TRUE, 
                                 top_marker = eqtl_variant)
  
  ggplot2::ggsave(stack_plot, filename = paste0("Locus_plots/",disease_name,"_",gene_name,".png"), width = 9, height = 3+3*1, dpi = 300, units = "in", limitsize = F)
  
  print("Plot finished")
  
  print(stack_plot)
}

lapply(1:nrow(coloc_results), function(x) make_locus_plot(x))

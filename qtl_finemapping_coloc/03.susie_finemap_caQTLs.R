rm(list = ls())
options(bitmapType="cairo")
library(data.table)
library(susieR)
library(tidyverse)
library(Rfast)
library(bigsnpr)
library(arrow)
library(coloc)

setDTthreads(8)

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("supply in order: \n
       1. eQTL summary statistics file path \n
       2. LD regions file path \n
       3. Path to folder containing LD matrices by chr (i.e data/UKBB_LDmatrices/) \n
       4. Chromosome \n
       5. N")
} else {
  eqtl_path = args[1]
  ld_regions = args[2]
  ld_path = args[3]
  name = args[4]
  N = as.numeric(args[5])
}

cat("This script is designed to work with the output of tensorQTL. \n
    The following column names are expected (Cap sensitive): \n
    phenotype_id \n
    variant_id \n
    start_distance \n
    af \n
    ma_samples \n
    ma_count \n
    pval_nominal \n
    beta \n
    slope_se \n")

wd = "/gchm_lab/collab/marlis_pj/coloc/" # Working directory
setwd(wd)

#Uncomment the following variables to run interactively
# eqtl_path = "/gchm/cd4_caQTL_analysis/variant_to_peak_QTL/run_012625_qc_aware_qsmooth_CPM_MAF5_FDR5_1MB/results/006_caQTLs/filtered_qsmooth_norm_cpm_1mb/cd4_qsmooth_cpm_chromatin_narrowpeaks.cis_qtl_pairs.chr21.csv"
# ld_regions ="regions/CD4T_chromatin.chr21.regions.txt"
# ld_path = "Ashkenazi_LD_matrices/"
# name="CD4T_chromatin"
# N = 367
NCORES <- nb_cores()

cat("Loading eQTL summary statistcs\n")

if(grepl("parquet",eqtl_path)){
  eqtl = as.data.table(arrow::read_parquet(eqtl_path))
  # The position needs to be extracted from the variant ID
  # Use str_extract to extract the numbers between ":" and "["
  eqtl$pos <- as.integer(str_extract(eqtl$variant_id, "(?<=:)[0-9]+(?=\\[)"))
  eqtl$chr = as.integer(str_extract(eqtl$variant_id, "[0-9]+(?=\\:)"))
  eqtl$a0 = str_extract(eqtl$variant_id, "(?<=,)[A-Za-z]+") # a0  = reference allele
  eqtl$a1 = str_extract(eqtl$variant_id, "(?<=\\])[A-Za-z]+(?=,)") # a1 alt allele & effect allele
  colnames(eqtl)[colnames(eqtl) == "slope"] = "beta"
} else{
  eqtl = fread(eqtl_path)
  eqtl$pos <- as.integer(str_extract(eqtl$variant_id, "(?<=:)[0-9]+(?=\\[)"))
  eqtl$chr = as.integer(str_extract(eqtl$variant_id, "[0-9]+(?=\\:)"))
  eqtl$a0 = str_extract(eqtl$variant_id, "(?<=,)[A-Za-z]+") # a0  = reference allele
  eqtl$a1 = str_extract(eqtl$variant_id, "(?<=\\])[A-Za-z]+(?=,)") # a1 alt allele & effect allele
  colnames(eqtl)[colnames(eqtl) == "slope"] = "beta"
}

# Load primary lead variant loci and regions
cat("Loading primary regions\n")

sig_regions <- fread(file = ld_regions, header = T, stringsAsFactors = FALSE, fill=T)

load_LD = function(CHR,LOWER,UPPER,ld_path,sumstats){
  
  cat("Loading LD and BIM files\n")
  
  REGION = paste(LOWER,UPPER,sep="-")
  file_path = paste0(ld_path, "chr",CHR, "/",REGION,"/",REGION)
  
  if(file.exists(paste0(file_path,".bk"))){
    # Attach the "bigSNP" object in R session
    obj.bigSNP <- snp_attach(paste0(file_path,".rds"))
  } else {
    cat("Generating bk & rds files\n")
    # Read from bed/bim/fam, it generates .bk and .rds files.
    snp_readBed(paste0(file_path,".bed"))
    # Attach the "bigSNP" object in R session
    obj.bigSNP <- snp_attach(paste0(file_path,".rds"))
  }
  
  G   <- obj.bigSNP$genotypes
  POS2 <- snp_asGeneticPos(obj.bigSNP$map$chromosome, obj.bigSNP$map$physical.pos)
  
  # Impute missing genotypes
  G2 = snp_fastImputeSimple(G, method = "mode")
  obj.bigSNP$genotypes = G2
  
  if(!file.exists(paste0(file_path,"_imputed.bed"))){
    #write imputed genotypes to bed file
    snp_writeBed(obj.bigSNP, paste0(file_path,"_imputed.bed"), ind.row = rows_along(G2), ind.col = cols_along(G2))
  }
  
  system(paste0("/nfs/sw/plink/plink-1.90b6.24/plink --bfile ",file_path,"_imputed --keep-allele-order --r square --out ",file_path,"_imputed"))
  
  # Load LD and bim file data
  LD <- fread(paste0(file_path,"_imputed.ld"))
  BIM <- fread(paste0(file_path,"_imputed.bim"))
  
  setnames(BIM, c("chr", "variant_id", "dk", "pos", "minor", "major"))
  colnames(LD) = BIM$variant_id
  rownames(LD) = BIM$variant_id
  
  if(file.exists(paste0(file_path,"_imputed.bk"))){
    # Attach the "bigSNP" object in R session
    obj.bigSNP <- snp_attach(paste0(file_path,"_imputed.rds"))
  } else {
    cat("Generating bk & rds files\n")
    # Read from bed/bim/fam, it generates .bk and .rds files.
    snp_readBed(paste0(file_path,"_imputed.bed"))
    # Attach the "bigSNP" object in R session
    obj.bigSNP <- snp_attach(paste0(file_path,"_imputed.rds"))
  }
  
  # Match variant effects to allele LD was calcualted for
  map <- setNames(obj.bigSNP$map[-3], c("chr", "variant_id", "pos", "a0", "a1"))
  sumstats_matched <- snp_match(sumstats, map, join_by_pos = T, remove_dups = F, strand_flip = F)
  colnames(sumstats_matched)[colnames(sumstats_matched) == "a0"] = "non_effect_allele"
  colnames(sumstats_matched)[colnames(sumstats_matched) == "a1"] = "effect_allele"
  
  cat("Finished \n")
  
  return(list(LD, sumstats_matched))
}

# Function to remove rows/columns with the most NAs iteratively
clean_LD_matrix <- function(mat) {
  # Continue removing rows/columns with NAs until no NA remains
  while (any(is.na(mat))) {
    # Count NAs per row and column
    na_count_rows <- rowSums(is.na(mat))
    na_count_cols <- colSums(is.na(mat))
    
    # Find the row or column with the maximum number of NAs
    max_na_row <- which.max(na_count_rows)
    max_na_col <- which.max(na_count_cols)
    
    # Determine whether to remove row or column (prioritize rows)
    if (na_count_rows[max_na_row] >= na_count_cols[max_na_col]) {
      # Remove the row and corresponding column (symmetric matrix)
      mat <- mat[-max_na_row, -max_na_row]
    } else {
      # Remove the column and corresponding row (symmetric matrix)
      mat <- mat[-max_na_col, -max_na_col]
    }
  }
  return(mat)
}

cat("Analysing GWAS", name, "\n")

run_susie_finemap = function(x){
  
  CHR = sig_regions$chr[x]
  LOWER = sig_regions$lower[x]
  UPPER = sig_regions$upper[x]
  REGION = sig_regions$finemap_regions[x]
  PHENO = str_split(sig_regions$phenotype_ids[x],",")[[1]]
  
  cat("Preparing to run SuSiE FINEMAP for region",sig_regions$finemap_regions[x], "index number",x, "\n")
  
  # Filter GWAS sumstats for SNPs in region
  snp_region = eqtl[chr == CHR & inrange(pos, LOWER, UPPER) & phenotype_id %in% PHENO]
  
  cat("Imputing missing genotypes and computing correlations \n")
  
  temp = load_LD(CHR = CHR,LOWER = LOWER,UPPER = UPPER,
                 ld_path = ld_path, sumstats = snp_region)
  
  LD_clean = temp[[1]]
  snp_region = as.data.table(temp[[2]])
  
  genes = unique(snp_region$phenotype_id)
  
  lbf_res_out = data.frame()
  cs_res_out = data.frame()
  no_convergence_out = data.frame()
  
  # Loop through genes
  for (name in genes){
    
    # Subset snp_region by the given name
    snp_region_gene <- snp_region[phenotype_id == name,]
    
    if(nrow(snp_region_gene) < 100){
      print("Not enough variants in region skipping")
      next
    }
    
    # Remove NAs in variant_id
    snp_region_gene <- snp_region_gene[!is.na(variant_id), ]
    
    # Calculate varbeta and z
    snp_region_gene[, varbeta := slope_se^2]
    snp_region_gene[, z := beta / slope_se]
    setkey(snp_region_gene, variant_id)
    
    # Remove duplicated variant IDs
    snp_region_gene <- snp_region_gene[!duplicated(snp_region_gene$variant_id),]
    
    # Get common variant IDs between LD_clean and snp_region_gene
    common_variants <- intersect(rownames(LD_clean), snp_region_gene$variant_id)
    
    # Subset the LD matrix (assuming LD_clean is a data.table)
    LD_clean2 <- LD_clean[rownames(LD_clean) %in% common_variants, , drop = FALSE]
    LD_clean2 <- LD_clean2[, colnames(LD_clean) %in% common_variants, with = FALSE]
    
    # Ensure the row names and column names are consistent
    rownames(LD_clean2) <- colnames(LD_clean2)
    
    # Convert to matrix
    LD_mat <- as.matrix(LD_clean2)
    rownames(LD_mat) = colnames(LD_mat)
    
    # Run the function on your LD matrix
    LD_mat_clean <- clean_LD_matrix(LD_mat)
    
    # Order snp_region_gene to match the LD matrix variant IDs
    snp_region_gene <- snp_region_gene[match(rownames(LD_mat_clean), snp_region_gene$variant_id),]
    
    # Here we loop through all the molecular_trait_ids within the region and perform susie
    cat("Running susie for region Chr",CHR,REGION,name,"\n")
    
    # Format SuSiE input
    D2 = list(
      snp = snp_region_gene$variant_id,
      position = snp_region_gene$pos,
      beta = snp_region_gene$beta,
      varbeta = snp_region_gene$varbeta,
      MAF = snp_region_gene$af,
      N = N,
      type = "quant",
      LD = LD_mat_clean
    )
    print(check_alignment(D2))
    #estimate_s_rss(snp_region_gene$z, LD_mat, n=N)
    
    S2 <- try(
      susie_rss(z = snp_region_gene$z,
                R = LD_mat_clean,
                n = N,
                L = 10,
                coverage = 0.95,
                estimate_residual_variance = TRUE,
                prior_variance = 50,
                check_prior = TRUE,
                max_iter = 1000,
                verbose = FALSE),
      silent = TRUE
    )
    
    lbf_res = data.frame()
    cs_res = data.frame()
    no_convergence = data.frame()
    
    if (inherits(S2, "try-error") || is.null(S2$sets$cs)) {
      
      warning("SuSiE failed or no credible sets found. Lowering coverage to 0.1.")
      S2 <- try(
        susie_rss(z = snp_region_gene$z,
                  R = LD_mat_clean,
                  n = N,
                  L = 10,
                  coverage = 0.10,
                  estimate_residual_variance = TRUE,
                  prior_variance = 50,
                  check_prior = TRUE,
                  max_iter = 1000,
                  verbose = FALSE),
        silent = TRUE
      )
      
      
      if (inherits(S2, "try-error") || is.null(S2$sets$cs)) {
        
        warning("Still no credible sets found for gene ", name)
        no_convergence_out = bind_rows(no_convergence_out,
                                       data.frame(region = paste0(CHR, ":", REGION), gene = name))
        next
        
      } else {
        message("SuSiE converged after lowering coverage.")
      }
    } else {
      message("Fine-mapping successful for gene ", name)
    }
    
    # Save lbf + PIP results
    pr = as.data.frame(t(S2$lbf_variable)) %>%
      mutate(variant_id = D2$snp, PIP = S2$pip,
             region = REGION, 
             peak = name,
             nsnps = length(D2$snp)) %>%
      left_join(snp_region_gene[,c("variant_id","variant_id.ss")], "variant_id") # rejoin unflipped variant ids

    colnames(pr) = c(paste0("lbf_", 1:10), "variant_id", "PIP", "region", "peak", "nsnps", "variant_id.ss")
    lbf_res_out = bind_rows(lbf_res_out, pr)
    
    # Save credible sets
    for (i in seq_along(S2$sets$cs_index)) {
      cs = as.data.frame(cbind(S2$pip[S2$sets$cs[[i]]], t(S2$lbf_variable[, S2$sets$cs[[i]]]))) %>%
        mutate(variant_id = D2$snp[S2$sets$cs[[i]]],
               chr = CHR,
               peak = name,
               region = REGION,
               cs = i,
               coverage = S2$sets$coverage[[i]],
               nsnps = length(D2$snp)) %>%
        left_join(snp_region_gene[,c("variant_id","variant_id.ss")], "variant_id") # rejoin unflipped variant ids
      
      colnames(cs) = c("pip", paste0("lbf_cs", 1:10), "variant_id", "chr", "peak", "region", "cs", "coverage", "nsnps", "variant_id.ss")
      cs_res_out = bind_rows(cs_res_out, cs)
    }
    
  } # end of gene loop
  
  message("Finished processing region ", REGION)
  
  return(list(lbf = lbf_res_out, cs = cs_res_out, no_convergence = no_convergence_out))
}

options(warn = 1)

run_susie_finemap2 = purrr::possibly(run_susie_finemap, otherwise = NA, quiet = FALSE)

results = lapply(1:nrow(sig_regions), function(x) run_susie_finemap2(x))

#cd4_atac_summits_peak_3433b

eqtl$unique_id = paste0(eqtl$variant_id, eqtl$phenotype_id)

safe_extract = function(lst, i) {
  if (is.list(lst) && length(lst) >= i) return(lst[[i]])
  return(NULL)
}

lbf_res_combined <- do.call(rbind, lapply(results, safe_extract, i = 1)) %>%
  mutate(unique_id = paste0(variant_id.ss, peak)) %>%
  left_join(eqtl, by = "unique_id") %>%
  dplyr::select(-phenotype_id, -unique_id)

cs_res_combined <- do.call(rbind, lapply(results, safe_extract, i = 2))

non_converged_regions <- do.call(rbind, lapply(results, safe_extract, i = 3))

if (nrow(cs_res_combined) > 0) {
  chr_tag = cs_res_combined$chr[1]
  
  dir.create(file.path("SuSiE_finemap_credible_sets", name), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path("SuSiE_finemap_lbf", name), showWarnings = FALSE, recursive = TRUE)
  
  fwrite(cs_res_combined, file.path("SuSiE_finemap_credible_sets", name, paste0(name, "_chr", chr_tag, "_credible_sets.txt")),
         sep = "\t", row.names = FALSE, quote = FALSE)
  
  fwrite(non_converged_regions, file.path("SuSiE_finemap_credible_sets", name, paste0(name, "_chr", chr_tag, "_non_converged_regions.txt")),
         sep = "\t", row.names = FALSE, quote = FALSE)
  
  fwrite(lbf_res_combined, file.path("SuSiE_finemap_lbf", name, paste0(name, "_chr", chr_tag, "_lbfs.txt")),
         sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  warning("cs_res_combined is empty — skipping file writing")
}

rm(list=ls())
libraries = c("data.table",
              "tidyverse",
              "coloc")

invisible(suppressMessages(lapply(libraries, require, character.only = TRUE)))

setDTthreads(8)

args = commandArgs(trailingOnly = TRUE)
 
if (length(args) < 2) {
 stop("Supply GWAS name and eQTL name")
} else {
 gwas_id = args[1]
 eqtl_name = args[2]
}

# gwas_id = "2020_32888494_lymph-count_EUR_preprocessed"
# eqtl_name = "All_CD4T_cells"

wd = "/gchm_lab/collab/marlis_pj/coloc/"

setwd(wd)

gwas_dir = '/gchm_lab/collab_gwas/finemapping_autoimmune/'

# Obtain list of gwas regions
regions <- system(
  paste0('ls ', gwas_dir, 'output/nathan_completed/', gwas_id, '/*_1.5Mb/susie/*.rds'),
  intern = TRUE
)
regions = regions[grepl(".rds", regions)]

# Load GWAS summary statistics
gwas_stats = fread(paste0(gwas_dir,"data/preprocessed_v1/",gwas_id,".tsv")) %>%
  mutate(variant_id_GRC37 = paste(chromosome, position, allele2, allele1,sep=":")) %>%
  mutate(variant_id_GRC37_alt = paste(chromosome, position, allele1, allele2,sep=":"))

if(nrow(gwas_stats) < 5e-06){
  stop("GWAS summary statistics contains less than 5 million variants")
}

region = regions[7]

run_coloc = function(region){
  
  cat("Loading region",region,"\n")
  
  # extract chromosome from region
  split = str_split(region,"_")[[1]]
  chr = str_split(split[length(split)],"\\.")[[1]][1]
  region_name = str_split(split[length(split)],"\\.")[[1]][2]
  
  res <- readRDS(region)  # or loop over all if needed
  
  # Extract the number of credible sets
  cs_gwas = res$sets
  idx_gwas = sort(cs_gwas$cs_index)
  
  if (length(idx_gwas) == 0 | res$converged == F) {
    idx_gwas = 1
  }
  
  cs_cols = paste0("lbf_",idx_gwas)
  
  # Extract GWAS lbfs from susie object
  # Left join GWAS summary statistic information
  lbf = as.data.frame(t(res$lbf_variable)) %>%
    rownames_to_column("variant_id_GRC37") %>%
    left_join(gwas_stats, "variant_id_GRC37") %>%
    select(-variant_id_GRC37_alt)
  
  # Extract PIP
  lbf$PIP = res$pip

  # Flip variants allele naming and left join GWAS sumstat info
  lbf2 = lbf %>% 
    filter(is.na(pval) == T) %>% # keep rows which didn't match variants
    select(-c(chromosome,position,allele1,allele2,beta,se,pval)) %>% # remove sumstat information
    dplyr::rename(variant_id_GRC37_alt = variant_id_GRC37) %>% # rename unmatched variants to alt
    left_join(gwas_stats, "variant_id_GRC37_alt") %>%
    select(-variant_id_GRC37) %>%
    dplyr::rename(variant_id_GRC37 = variant_id_GRC37_alt)

  gwas_lbf = rbind(lbf[is.na(lbf$pval) == F,], lbf2)
  colnames(gwas_lbf)[2:11] = paste0("lbf_",seq(1,10))
  # Filter lbfs to number of credible sets
  gwas_lbf = gwas_lbf[, c("variant_id_GRC37", cs_cols, "PIP", "chromosome",  "position", "allele1", "allele2",  "beta", "se", "pval")]
  
  # Read eQTL lbfs and cre variant ids
  eqtl_lbf = fread(paste0("SuSiE_finemap_lbf_GRC37/",eqtl_name,"/",eqtl_name,"_chr",chr,"_lbfs_GRC37.txt")) %>%
    mutate(variant_id_GRC37 = paste(gsub("chr","",chr), pos_GRC37, a0, a1, sep = ":"),
           variant_id_GRC37_alt = paste(gsub("chr","",chr), pos_GRC37, a1, a0, sep = ":")) %>%
    filter(variant_id_GRC37 %in% gwas_lbf$variant_id_GRC37 | 
             variant_id_GRC37_alt %in% gwas_lbf$variant_id_GRC37)
  
  if(nrow(eqtl_lbf) < 1){
    print("no significant eQTLs in this region, skipping")
    return(NA)
  }
  
  cat("Number of eQTL variants after filtering:", nrow(eqtl_lbf), "\n")
  
  # Load eQTL credible sets
  cs_df = fread(paste0("SuSiE_finemap_credible_sets/",eqtl_name,"/",eqtl_name,"_chr",chr,"_credible_sets.txt")) 
  
  # Filter GWAS variants to overlapping eQTL variants
  gwas_lbf2 = gwas_lbf[gwas_lbf$variant_id_GRC37 %in% eqtl_lbf$variant_id_GRC37 | 
                         gwas_lbf$variant_id_GRC37 %in% eqtl_lbf$variant_id_GRC37_alt,]
  
  genes = unique(eqtl_lbf$gene)
  
  coloc_res = data.frame()
  
  for (mol_id in genes){
    
    cat("Running coloc for", mol_id,"\n")
    
    # Format eQTL lbf df
    eqtl_lbf2 = eqtl_lbf[gene == mol_id,] %>%
      mutate(ref_alt = ifelse(variant_id_GRC37 %in% gwas_lbf2$variant_id_GRC37,1,0)) %>%
      mutate(variant_id_GRC37 = ifelse(ref_alt == 0, variant_id_GRC37_alt, variant_id_GRC37)) %>%
      select(-variant_id_GRC37_alt) %>% 
      distinct(variant_id_GRC37, .keep_all = T)
    
    gwas_lbf3 = gwas_lbf2[gwas_lbf2$variant_id_GRC37 %in% eqtl_lbf2$variant_id_GRC37,]
    
    if(nrow(gwas_lbf3) < 300){
      cat("Number of overlapping variants too low",nrow(gwas_lbf3), "\n")
      next
    }
    
    # Make DFs have the same order
    eqtl_lbf2 = eqtl_lbf2[match(eqtl_lbf2$variant_id_GRC37, gwas_lbf3$variant_id_GRC37),]
    gwas_lbf3 = gwas_lbf3[match(gwas_lbf3$variant_id_GRC37, eqtl_lbf2$variant_id_GRC37),]
    
    # Extract credible sets
    idx_eqtl = cs_df %>% 
      filter(region %in% eqtl_lbf2$region_GRC38 & gene %in% eqtl_lbf2$gene) %>%
      filter(gene %in% eqtl_lbf2$gene) %>%
      distinct(cs) %>% unlist()
    
    if (length(idx_eqtl) == 0) {
      idx_eqtl = 1
      eqtl_converged  = FALSE
    } else{
      eqtl_converged  = T
    }
    
    
    # remove lbfs which have no credible sets found
    cs_eqtl_col = paste0("lbf_",seq(1,10)[!seq(1,10) %in% idx_eqtl])
    eqtl_lbf2 = eqtl_lbf2[, setdiff(colnames(eqtl_lbf2), cs_eqtl_col), with = FALSE]
    
    ## Convert LBF data frames into matrices suitable for the coloc.bf_bf() method
    eqtl_mat = as.matrix(eqtl_lbf2[, .SD, .SDcols = grepl("lbf", colnames(eqtl_lbf2))])
    row.names(eqtl_mat) = eqtl_lbf2$variant_id_GRC37
    eqtl_mat = t(eqtl_mat)
    
    gwas_mat = as.matrix(gwas_lbf3[, grepl("lbf", colnames(gwas_lbf3))])
    row.names(gwas_mat) = gwas_lbf3$variant_id_GRC37
    gwas_mat = t(gwas_mat)
    
    if( res$converged == F | eqtl_converged == F ){
      
      print("Running single variant coloc")
      
      lbf1 <- gwas_mat[1, ]
      lbf2 <- eqtl_mat[1, ]
      
      # (a) single-SNP posterior within H4
      sumlbf <- lbf1 + lbf2
      snp_pp  <- exp(sumlbf - coloc:::logsum(sumlbf))   # same as SNP.PP.H4
      
      # (b) region-level hypothesis posteriors
      pp <- coloc:::combine.abf(
        lbf1, lbf2,
        p1  = 1e-4,
        p2  = 1e-4,
        p12 = 1e-5
      )
      
      temp = data.frame(nsnps = ncol(gwas_mat), 
                        hit1 = names(snp_pp)[snp_pp == max(snp_pp)],
                        hit2 = names(snp_pp)[snp_pp == max(snp_pp)])
      
      # Convert pp to a 1-row data.frame
      pp_df <- as.data.frame(as.list(pp))
      
      # Combine
      coloc <- cbind(temp, pp_df)
                          
      coloc$idx1 = 1
      coloc$idx2 = 1
      
    } else{
      
      print("Running multi variant coloc")
      
      coloc = coloc::coloc.bf_bf(gwas_mat,eqtl_mat)
      coloc = dplyr::as_tibble(coloc$summary)
      
    }

    if(max(coloc$PP.H4.abf,na.rm=T) < 0.5 ){
      cat(mol_id,"No significant colocalisation, PP.H4 = ",max(coloc$PP.H4.abf,na.rm=T), "\n")
      next
    }
    
    coloc = coloc %>% dplyr::filter(PP.H4.abf > 0.5) %>% 
      rename(GWAS_variant_GRC37 = hit1, eQTL_variant_GRC37 = hit2)
    
    cat(mol_id,"Wooo Significant colocalisation! \n")
    
    temp_df = gwas_lbf3 %>% rename(GWAS_variant_GRC37 = variant_id_GRC37)
    temp_df2 = eqtl_lbf2 %>% rename(eQTL_variant_GRC37 = variant_id_GRC37)
    
    coloc = coloc %>%
      left_join(temp_df, "GWAS_variant_GRC37") %>%
      left_join(temp_df2,"eQTL_variant_GRC37") %>%
      mutate(region = paste(chr, region_name,sep = "."))
    
    coloc_res = bind_rows(coloc_res,coloc)
    
  }
  
  return(coloc_res)
  
}

run_coloc = purrr::possibly(run_coloc, otherwise = NA, quiet = F)
results = lapply(1:length(regions), function(x) run_coloc(regions[x]))

# Filter results to remove NA and 0-column data frames
results_filtered <- results[!sapply(results, function(x) {
  is.null(x) || (is.atomic(x) && any(is.na(x))) || (is.data.frame(x) && ncol(x) == 0)
})]

if (length(results_filtered) > 0) {
  results2 <- dplyr::bind_rows(results_filtered)
} else {
  results2 <- data.frame()  # or NA if that's more appropriate
}

if(nrow(results2) > 0){
  
  print("Writing file")
  
  fwrite(results2, paste0("coloc_results/",gwas_id, "_coloc_results.csv"), quote = F, row.names = F, sep = "\t")
}

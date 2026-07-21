rm(list=ls())
libraries = c("data.table",
              "tidyverse",
              "coloc")

invisible(suppressMessages(lapply(libraries, require, character.only = TRUE)))

setDTthreads(8)

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Supply caQTL name and eQTL name")
} else {
  caqtl_name = args[1]
  eqtl_name = args[2]
  chr = args[3]
}

# caqtl_name = "CD4T_chromatin"
# eqtl_name = "All_CD4T_cells"
# chr = 19

wd = "/gchm_lab/collab/marlis_pj/coloc/"

setwd(wd)

# Load eQTL lbfs
eqtl_lbf = fread(paste0("SuSiE_finemap_lbf/",eqtl_name,"/",eqtl_name,"_chr",chr,"_lbfs.txt"))
# extract regions
regions = unique(eqtl_lbf$region)
# load eQTL credible sets
eqtl_cs = fread(paste0("SuSiE_finemap_credible_sets/",eqtl_name,"/",eqtl_name,"_chr",chr,"_credible_sets.txt"))

# load caQTL lbfs
caqtl_lbfs = fread(paste0("SuSiE_finemap_lbf/",caqtl_name,"/",caqtl_name,"_chr",chr,"_lbfs.txt"))
# load caQTLs credible sets
caqtl_cs = fread(paste0("SuSiE_finemap_credible_sets/",caqtl_name,"/",caqtl_name,"_chr",chr,"_credible_sets.txt"))

# Start function here
coloc_region = regions[1]

run_coloc = function(coloc_region){
  
  # filter eQTL lbfs to those within region
  eqtl_region = eqtl_lbf[region == coloc_region,]
  
  # extract positions from region
  split = str_split(coloc_region,":")[[1]][2]
  lower = as.integer(str_split(split, "-")[[1]][1])
  upper = as.integer(str_split(split, "-")[[1]][2])
  
  # filter caQTL lbfs to region
  caqtl_region = caqtl_lbfs[inrange(pos, lower, upper)]
  
  # loop through all genes in eQTL region
  genes = unique(eqtl_region$gene)
  
  coloc_res_genes = data.frame()
  
  for (mol_id in genes){
    
    cat("Running coloc for", mol_id,"\n")
    
    # filter eqtl lbf to gene
    eqtl_gene = eqtl_lbf[gene == mol_id,]
    # filter to eqtl credible sets
    eqtl_gene_cs = sort(eqtl_cs[region == coloc_region & gene == mol_id,cs])
    cs_eqtl_col = paste0("lbf_",seq(1,10)[!seq(1,10) %in% eqtl_gene_cs])
    eqtl_gene = eqtl_gene[, setdiff(colnames(eqtl_gene), cs_eqtl_col), with = FALSE]
    
    caqtl_gene = caqtl_region[variant_id.x %in% eqtl_gene$variant_id.x,]
    
    if(nrow(caqtl_gene) < 200){
      cat("Number of overlapping variants too low",nrow(caqtl_gene), "\n")
      next
    }
    
    peaks = unique(caqtl_gene$peak)
    
    coloc_res_peaks = data.frame()
    
    # loop through peaks
    for (peak_name in peaks){
      
      cat("Running coloc for", mol_id, "and",peak_name,"\n")
      
      caqtl_peak = caqtl_gene[peak == peak_name,]
      # obtain credible sets
      caqtl_cs_peak = sort(unique(caqtl_cs[peak == peak_name,cs]))
      cs_caqtl_col = paste0("lbf_",seq(1,10)[!seq(1,10) %in% caqtl_cs_peak])
      caqtl_peak = caqtl_peak[, setdiff(colnames(caqtl_peak), cs_caqtl_col), with = FALSE]
      
      eqtl_peak = eqtl_gene[variant_id.x %in% caqtl_peak$variant_id.x,]
      
      # Make DFs have the same order
      caqtl_peak = caqtl_peak[match(caqtl_peak$variant_id.x, eqtl_peak$variant_id.x),]
      eqtl_peak = eqtl_peak[match(eqtl_peak$variant_id.x, caqtl_peak$variant_id.x),]
      
      ## Convert LBF data frames into matrixes suitable for the coloc.bf_bf() method
      eqtl_mat = as.matrix(eqtl_peak[, .SD, .SDcols = grepl("lbf", colnames(eqtl_peak))])
      row.names(eqtl_mat) = eqtl_peak$variant_id.x
      eqtl_mat = t(eqtl_mat)
      
      caqtl_mat = as.matrix(caqtl_peak[, .SD, .SDcols = grepl("lbf", colnames(caqtl_peak))])
      row.names(caqtl_mat) = caqtl_peak$variant_id.x
      caqtl_mat = t(caqtl_mat)
      
      coloc = coloc::coloc.bf_bf(caqtl_mat,eqtl_mat, overlap.min = 1)
      coloc = dplyr::as_tibble(coloc$summary)
      
      if(max(coloc$PP.H4.abf,na.rm=T) < 0.5 ){
        cat(mol_id,"No significant colocalisation, PP.H4 = ",max(coloc$PP.H4.abf,na.rm=T), "\n")
        next
      }
      
      coloc = coloc %>% dplyr::filter(PP.H4.abf > 0.5) %>% 
        rename(caQTL_variant = hit1, eQTL_variant = hit2)
      
      cat(mol_id,"Wooo Significant colocalisation! \n")
      
      temp_df = caqtl_peak %>% rename(caQTL_variant = variant_id.x)
      temp_df2 = eqtl_peak %>% rename(eQTL_variant = variant_id.x)
      
      coloc = coloc %>%
        left_join(temp_df, "caQTL_variant") %>%
        left_join(temp_df2,"eQTL_variant")
      
      coloc_res_peaks = bind_rows(coloc_res_peaks,coloc)
      
    }
    
    coloc_res_genes = bind_rows(coloc_res_genes, coloc_res_peaks)
    
  }
  
  return(coloc_res_genes)

}

run_coloc = purrr::possibly(run_coloc, otherwise = NA, quiet = F)
results = lapply(1:length(regions), function(x) run_coloc(regions[x]))
results2 <- dplyr::bind_rows(results)
# remove duplicate columns
fwrite(results2, paste0("coloc_results/ca_eqtl_coloc/chr",chr, "_coloc_results.csv"), quote = F, row.names = F, sep = "\t")

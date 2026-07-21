libraries <- c("coloc",
               "susieR",
               "data.table",
               "dplyr",
               "bigsnpr",
               "geni.plots")

invisible(suppressMessages(lapply(libraries, require, character.only = TRUE)))

setDTthreads(8)

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Supply necessary files")
} else {
  gwas_id = args[1]
  quant = args[2]
}

wd = "/gchm_lab/collab/marlis_pj/coloc/"

setwd(wd)

if(quant == T){
  type = "quant"
} else {
  type = "cc"
}

dataset_id = "GCST90435144"
type= "quant"

load_LD = function(CHR,LOWER,UPPER,ld_path,sumstats){
  
  cat("Loading LD and BIM files\n")
  
  REGION = paste(LOWER,UPPER,sep=".")
  
  if(file.exists(paste0(ld_path, "chr",CHR, "/",REGION,"/",REGION,".bk"))){
    # Attach the "bigSNP" object in R session
    obj.bigSNP <- snp_attach(paste0(ld_path, "chr",CHR, "/",REGION,"/",REGION,".rds"))
  } else {
    cat("Generating bk & rds files\n")
    # Read from bed/bim/fam, it generates .bk and .rds files.
    snp_readBed(paste0(ld_path, "chr",CHR, "/",REGION,"/",REGION,".bed"))
    # Attach the "bigSNP" object in R session
    obj.bigSNP <- snp_attach(paste0(ld_path, "chr",CHR, "/",REGION,"/",REGION,".rds"))
  }
  
  # Match variant effects to allele LD was calcualted for
  map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a0", "a1"))
  sumstats_matched <- snp_match(sumstats, map, join_by_pos = T, remove_dups = F)
  
  # Define filenames
  LDfilename <- paste0(ld_path, "chr",CHR, "/",REGION,"/",REGION)
  BIMfilename <- paste0(ld_path, "chr",CHR, "/",REGION,"/",REGION)
  # Load LD and bim file data
  LD <- fread(paste(LDfilename, "ld", sep = "." ))
  BIM <- fread( paste(BIMfilename, "bim", sep = ".") )
  # assign headers
  setnames(BIM, c("chr", "variant_id", "dk", "pos", "a0", "a1"))
  # assign SNP labels to LD matrix
  setnames(LD, BIM$variant_id)
  LD[, variant_id := BIM$variant_id]
  
  # Remove columns from the LD matrix that are all NAs
  na_columns = colSums(is.na(LD)) == nrow(LD)
  nacol_names = names(na_columns[na_columns == TRUE])
  LD = LD[, !names(LD) %in% nacol_names, with = FALSE]
  # Filter rows with complete cases
  LD = LD[complete.cases(LD)]
  # Filter snps in LD matrix for those in GWAS sum stats
  LD <- LD[LD$variant_id %in% sumstats_matched$variant_id,]
  # Make sure LD only contains columns in LD$SNP and "SNP"
  required_cols <- unique(c(LD$variant_id,"variant_id"))
  LD <- LD[, ..required_cols]
  # Remove duplicated SNPs
  LD_clean <- unique(LD, by = "variant_id")
  # Remove duplicates from BIM file
  BIM_clean <- BIM[!duplicated(BIM, by = "variant_id"), ]
  
  cat("Finished \n")
  return(list(LD_clean, sumstats_matched))
}

run_coloc = function(dataset_id){
  
  #translated_gwas_name = dict[gwas_id]
  
  ## Import lbf variable values from eQTL catalogue into R
  gwas_file = fread(paste0("GWAS_sumstats/",dataset_id,".h.tsv.gz"))
  colnames(gwas_file)[1:4] = c("chr","pos","a1","a0")
  gwas_file$variant_id = paste0(gwas_file$chr,":", gwas_file$pos, "[b38]",gwas_file$a0,",",gwas_file$a1)
  gwas_id = stringr::str_split_fixed(dataset_id,"\\.",4)[,1]
  
  coloc_res = data.frame()
  
  for (i in 1:22){
    
    cat("Running coloc between GWAS",gwas_id,"&","All_CD4T_cells chr",i,"\n")
    
    eqtl_lbf = fread(paste0("SuSiE_finemap_lbf/All_CD4T_cells/All_CD4T_cells_chr",i,"_lbfs.txt"),fill=T)
  
    gwas = gwas_file[variant_id %in% eqtl_lbf$variant_id.x,]
    
    if (min(gwas$p_value, na.rm=T) > 1e-08){
      next("No GWS significant variants in this region for this GWAS")
    }
    
    # Loop through regions
    regions = unique(eqtl_lbf$region)
    
    for (ld_region in regions){
      
      cat("Region",ld_region,"\n")
      
      lower = stringr::str_split_fixed(stringr::str_extract(ld_region, "(?<=:).*"),"-",2)[,1]
      upper = stringr::str_split_fixed(stringr::str_extract(ld_region, "(?<=:).*"),"-",2)[,2]
      
      eqtl_lbf2 = eqtl_lbf[region == ld_region,]
      gwas2 = gwas[variant_id %in% eqtl_lbf2$variant_id.x,]
      
      if (min(gwas2$p_value, na.rm=T) > 1e-08){
        next("No GWS significant variants in this region for this GWAS")
      }
      
      if(nrow(gwas2) < 500){
        cat("Number of overlapping variants too low",nrow(gwas2), "\n")
        next
      }
      
      system(paste0("mkdir -p ","1kg_LD/chr",i,"/",lower,".",upper))
      
      # Write variants to txt file
      fwrite(list(eqtl_lbf2$variant_id.x), paste0("1kg_LD/chr",i,"/",lower,".",upper,"/variants_to_extract.txt"), 
             quote = F, row.names = F, sep = " ", col.names = F)
      
      if(!file.exists(paste0("1kg_LD/", "chr",i, "/",lower,".",upper,"/",lower,".",upper,".ld"))){
        # Submit the job and capture the job ID
        job_submission <- system(paste("sbatch 05.Generate_LD_1000genomes.sh", i, lower, upper, sep = " "), intern = TRUE)
        job_id <- sub("Submitted batch job ", "", job_submission)
        
        # Wait for Job Completion
        # Function to check if the job is still running
        is_job_running <- function(job_id) {
          job_status <- system(paste("squeue -j", job_id), intern = TRUE)
          return(length(job_status) > 1)  # If there's more than one line, the job is still running
        }
        
        # Wait for the job to finish
        while (is_job_running(job_id)) {
          Sys.sleep(10)  # Wait for 10 seconds before checking again
        }
      }
      
      # Once the job is done, proceed with the next command
      temp <- load_LD(CHR = i, LOWER = lower, UPPER = upper, ld_path = "1kg_LD/", sumstats = gwas2)
      
      gwas_LDmat = temp[[1]] %>% tibble::column_to_rownames("variant_id")
      gwas_sumstats = as.data.table(temp[[2]])
      
      # Finemap GWAS region
      gwas_sumstats[, varbeta := standard_error^2]
      gwas_sumstats[, z := beta / standard_error]
      
      S2 = try(susie_rss(z = gwas_sumstats$z, R = as.matrix(gwas_LDmat), gwas_sumstats$n[1],
                         estimate_residual_variance = FALSE,
                         prior_variance = 50,
                         check_prior = TRUE))
      
      if (inherits(S2, "try-error") || is.null(S2$sets) || S2$converged == F) {
        
        print("No credible sets found")
        next
        
      } else {
        
        gwas_lbf = as.data.frame(t(S2$lbf_variable)) %>% tibble::rownames_to_column("variant_id") %>%
          left_join(gwas_sumstats,"variant_id")
        colnames(gwas_lbf)[1:11] = c( "variant_id" ,paste("lbf_",seq(1,10), sep=""))
        
        genes = unique(eqtl_lbf2$molecular_id)
        
        for (gene in genes){
          
          cat("Gene",gene,"\n")
          
          eqtl_gene = eqtl_lbf2[molecular_id == gene]
          
          gwas_gene_region = gwas_lbf[gwas_lbf$variant_id %in% eqtl_gene$variant_id.x,]
          
          eqtl_gene2 = eqtl_gene[variant_id.x %in% gwas_gene_region$variant_id,]
          
          if(nrow(eqtl_gene2) < 500){
            cat("Number of overlapping variants too low",nrow(eqtl_gene2), "\n")
            next
          }
          
          ## Convert LBF data frames into matrixes suitable for the coloc.bf_bf() method
          eqtl_mat = as.matrix(dplyr::select(eqtl_gene2, lbf_1:lbf_10))
          row.names(eqtl_mat) = eqtl_gene2$variant_id.x
          eqtl_mat = t(eqtl_mat)
          
          gwas_mat = as.matrix(dplyr::select(gwas_gene_region, lbf_1:lbf_10))
          row.names(gwas_mat) = gwas_gene_region$variant_id
          gwas_mat = t(gwas_mat)
          
          coloc = coloc::coloc.bf_bf(gwas_mat,eqtl_mat)
          coloc = dplyr::as_tibble(coloc$summary) %>% 
            rename(variant_id = hit1, variant_id.x = hit2) %>% dplyr::filter(PP.H4.abf > 0.5) 
          
          coloc = coloc %>% mutate(dataset_id = dataset_id) %>%
            left_join(eqtl_gene2[,c("variant_id.x","region","pval_nominal", "PIP")], "variant_id.x") %>% 
            left_join(gwas_sumstats[,c("variant_id","p_value")],"variant_id") %>%
            rename(gwas_variant = variant_id, eqtl_variant = variant_id.x, eqtl_PIP = PIP)%>%
            mutate(gene_id = gene) %>% filter(pval_nominal < 1e-03 & p_value < 1e-05 & eqtl_PIP > 0.01)
          
          if(nrow(coloc) == 0){
            cat(gene,"No significant colocalisation \n")
            next
          }
          
          cat("Sucessful Colocalization Ploting region",ld_region,gene,"\n")
          
          identical(gwas_gene_region$variant_id, eqtl_gene2$variant_id.x)
          
          markers <- data.frame(marker = gwas_gene_region$variant_id,
                                chr = i,
                                pos = eqtl_gene2$pos,
                                pvalue_1 = eqtl_gene2$pval_nominal,
                                pvalue_2 = gwas_gene_region$p_value)
          
          required_cols <- colnames(gwas_LDmat)[colnames(gwas_LDmat) %in% gwas_gene_region$variant_id]
          LD_plot <- gwas_LDmat[, required_cols]
          LD_plot = LD_plot[required_cols,]
          
          markers <- markers[match(colnames(LD_plot), markers$marker), ]
          
          snp = coloc$eqtl_variant[which(coloc$PP.H4.abf == max(coloc$PP.H4.abf, na.rm = T))]
          
          stack_plot <- fig_region_stack(data = markers, trait = c(gene, dataset_id), corr = LD_plot,
                                         top_marker = snp, build =38, title_center = TRUE)
          
          dir_path <- paste0("Locus_plots/All_CD4T_cells/", dataset_id)
          
          # Check if the directory exists before attempting to create it
          if (!file.exists(dir_path)) {
            dir.create(dir_path, recursive = TRUE)
          }
          
          ggplot2::ggsave(stack_plot, filename = paste0(dir_path,"/chr",ld_region,"_",gene,".png"), width = 9, height = 4+4*1, dpi = 300, units = "in", limitsize = F)
          
          coloc_res = bind_rows(coloc_res,coloc)
          
        }
        
      }
      
    }
    
  }
  return(coloc_res)
}

options(warn=1)
### In case of an error return NA ### 
#run_coloc = purrr::possibly(run_coloc, otherwise = NA, quiet = F)
result = run_coloc(gwas_id)
# Write results
res_dir = paste0("coloc_results/")
system(paste0("bash -c 'mkdir -p ", res_dir, "'"))
fwrite(result, paste0(res_dir,gwas_id,"_coloc_results.txt"), sep = "\t", row.names = F, quote=F)
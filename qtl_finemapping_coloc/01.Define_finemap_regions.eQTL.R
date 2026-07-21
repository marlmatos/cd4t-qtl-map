libraries <- c("data.table",
               "tidyverse",
               "arrow")

invisible(suppressMessages(lapply(libraries, require, character.only = TRUE)))

setDTthreads(8)

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("supply cond_indep file path & name \n")
} else {
  cond_indep_path = args[1]
  name = args[2]
  dir = args[3]
}

cond_indep_path = "/gchm/cd4_QTL_analysis/04_Run_interaction_eQTLs/interaction_TF_disrupted_chromBPNet/results/interaction_chrombpnet_tfs/TF_activity_interaction_chrombpnet_tfs/ATF1/ATF1_chr20.cis_qtl_top_assoc.txt.gz"
#cond_indep_path = "/gchm/cd4_QTL_analysis/03_Run_cisQTL_perchr/analysis/cis_eQTLs/all_CD4T_cells_MAF5/allcells.cis_qtl_pairs.chr1.csv"
name="ATF1"
dir = "/gchm_lab/collab/marlis_pj/coloc"

setwd(dir)

cat("Loading", name, "\n")

# Read conditionally independent SNPs
if(grepl("parquet",cond_indep_path)){
  cond_indep = as.data.table(arrow::read_parquet(cond_indep_path))
} else{
  cond_indep = fread(cond_indep_path)
}

cond_indep$pos <- as.integer(str_extract(cond_indep$variant_id, "(?<=:)[0-9]+(?=\\[)"))
cond_indep$chr = as.integer(str_extract(cond_indep$variant_id, "[0-9]+(?=\\:)"))
cond_indep = cond_indep[order(pos),]
cond_indep[, lower := pos - 250000]
cond_indep[, upper := pos + 250000]
# Set negative values to zero
cond_indep$lower[cond_indep$lower < 0] = 0
# Calculate differences between regions
cond_indep[, difference := lower - shift(upper, 1, type = "lag"),]
cond_indep$difference[is.na(cond_indep$difference) == T] = 1
# Group variants with overlaping regions
cond_indep[, group := cumsum(difference > 0) + 1]
# Calculate lower and upper bounds of grouped region
cond_indep[, `:=` (regions = paste0(min(lower),"-", max(upper))), by = group]
cond_indep$lower[cond_indep$lower < 0] = 0

# Aggregate molecular ids per region
cond_indep_unique <- cond_indep[, .(
  phenotype_ids = paste(unique(phenotype_id), collapse = ",")
), by = .(regions)]

cond_indep = cond_indep %>% 
  distinct(regions, .keep_all=T) %>%
  left_join(cond_indep_unique, by = "regions") %>%
  mutate(finemap_regions = paste0(chr,":", regions))

cond_indep$lower = as.integer(str_split_fixed(cond_indep$regions,"-",2)[,1])
cond_indep$upper = as.integer(str_split_fixed(cond_indep$regions,"-",2)[,2])

setcolorder(cond_indep, c("chr", "lower", "upper", "finemap_regions","phenotype_ids",
                          setdiff(names(cond_indep), c("chr", "lower", "upper", "finemap_regions", "phenotype_ids"))))

if("V1" %in% colnames(cond_indep)){
  cond_indep[, c("V1","phenotype_id") := NULL]
} else{
  cond_indep[, c("phenotype_id") := NULL]
}

finemap_regions = unique(cond_indep$regions)

if (nrow(cond_indep[chr == 6 & upper >= 28500000 & 33500000 >= lower,]) > 0) {
  
  print("Truncating MHC region")
  
  # Filter to 
  tmp = cond_indep[chr == 6 & upper >= 28500000 & 33500000 >= lower,]
  
  tmp[, upper := 28500000 - 1]
  tmp[, lower := 33500000 + 1]
  
  tmp2 = cond_indep[!finemap_regions %in% tmp$finemap_regions,]
  
  cond_indep = bind_rows(tmp,tmp2)
  
}

print("Finished, writing regions file")

system("mkdir -p regions")

print(paste0("regions/", name,".chr",cond_indep$chr[1],".regions.txt"))

fwrite(cond_indep, paste0("regions/", name,".chr",cond_indep$chr[1],".regions.txt"), sep = "\t", row.names = F, quote=F)
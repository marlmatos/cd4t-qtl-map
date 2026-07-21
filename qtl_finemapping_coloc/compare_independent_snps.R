library(data.table)
library(tidyverse)

all = data.table()

for (i in c(1:22)){
  tmp = fread(paste0("/gchm/cd4_QTL_analysis/03_Run_cisQTL_perchr/analysis/cis_eQTLs/all_CD4T_cells_MAF5/allcells.independent_cis_qtl_pairs.chr",i,".csv"))
  all = rbind(all,tmp)
}

all2 = all %>%
  filter(pval_perm < 0.05) %>%
  distinct()

all_finemap = data.table()

for (i in c(1:22)){
  tmp = fread(paste0("/gchm_lab/collab/marlis_pj/coloc/SuSiE_finemap_credible_sets/All_CD4T_cells/All_CD4T_cells_chr",i,"_credible_sets.txt"))
  all_finemap = rbind(all_finemap,tmp)
}

all_finemap2 = all_finemap %>%
  mutate(unique_id = paste(region, gene, cs)) %>%
  distinct(unique_id, .keep_all = T)
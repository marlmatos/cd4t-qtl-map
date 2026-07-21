library(data.table)
library(tidyverse)

setwd("/gchm_lab/collab/marlis_pj/coloc")

# first get average size of regions

region_size = data.frame()

for (i in 1:22){
  # Load credible set
  region = fread(paste0("regions/All_CD4T_cells.chr",i,".regions.txt")) %>%
    mutate(diff = (upper-lower)/1000000)
  
  region_size = rbind(region_size, region)
}

max(region_size$diff)
median(region_size$diff)

size = data.frame()

for (i in 1:22){
  # Load credible set
  cred_set = fread(paste0("SuSiE_finemap_credible_sets/All_CD4T_cells/All_CD4T_cells_chr",i,"_credible_sets.txt")) %>%
    group_by(gene, region, cs) %>%
    summarise(n = n())
  
  size = rbind(size, cred_set)
}

median(size$n)
IQR(size$n)

hist(size$n)

size2 = size[size$n < 1000,]

median(size2$n)
IQR(size2$n)

hist(size2$n)

## Same for chromain QTLs

region_size = data.frame()

for (i in 1:22){
  # Load credible set
  region = fread(paste0("regions/CD4T_chromatin.chr",i,".regions.txt")) %>%
    mutate(diff = (upper-lower)/1000000)
  
  region_size = rbind(region_size, region)
}

min(region_size$diff)
max(region_size$diff)
median(region_size$diff)

size = data.frame()

for (i in 1:22){
  # Load credible set
  cred_set = fread(paste0("SuSiE_finemap_credible_sets/CD4T_chromatin/CD4T_chromatin_chr",i,"_credible_sets.txt")) %>%
    group_by(peak, region, cs) %>%
    summarise(n = n())
  
  size = rbind(size, cred_set)
}

median(size$n)
IQR(size$n)

hist(size$n)

size2 = size[size$n < 1000,]

median(size2$n)
IQR(size2$n)

hist(size2$n)
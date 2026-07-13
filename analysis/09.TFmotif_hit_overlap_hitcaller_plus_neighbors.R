################################
## Script to Process overlap MosDisco Motif Instances
## (plus neighbors)
## Author: Marliette Matos
## Date: 10/13/2025
################################s
library(readr)

#variant_hits<-read_tsv("~/cd4_chrombpnet/chrombpnet_model_b7/motif_hit_calls/motifs_hits_FILTVARS/variant_hit_calls.tsv")
variant_hits_filt<-read_tsv("~/cd4_qtl_paper_figures/figure_4/data/motif_variant_overlap_hit_caller/filtered_variant_hits.tsv") 

variant_hits_filt<-variant_hits_filt %>% mutate(hit_id=paste0(peak_id,"_", seqnames,"_", start, "_",end, "_",pattern_name))
vars<- variant_hits_filt %>% distinct(hit_id, variant_loc) 

#read unique hits and strait away filter out peaks that did not make the filtering in variant_hit_filt
#also filter other motifs instaces that do not met the reqs
unique_hits<-read_tsv("~/cd4_chrombpnet/chrombpnet_model_b7/motif_hit_calls/motifs_hits_FILTVARS/hits_unique.tsv") %>% filter(peak_id %in% variant_hits_filt$peak_id) %>% filter(hit_coefficient>5) #instance version

modisco_out40<-read_tsv("~/cd4_qtl_paper_figures/figure_4/data/motif_variant_overlap_hit_caller/modisco_motids_minseq40.tsv") %>% dplyr::select(pattern, num_seqlets, match_confidence, TF_best_match_code_friendly_unique_unique, family)

unique_hits<-unique_hits %>% right_join(modisco_out40, join_by("motif_name"=="pattern")) # I only want to keep high quality patterns

#annotate variants hits
unique_hits<-unique_hits %>% 
  mutate(hit_id=paste0(peak_id,"_",chr,"_", start, "_",end, "_",motif_name)) %>% 
  left_join(vars, by="hit_id")  %>% 
  mutate(variant_hit=ifelse(is.na(variant_loc), FALSE, TRUE))
  

table(unique_hits$variant_hit)
write_tsv(unique_hits, "~/cd4_qtl_paper_figures/figure_4/data/motif_variant_overlap_hit_caller/unique_hits_filtered_motif_instances.tsv")

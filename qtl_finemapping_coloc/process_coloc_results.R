library(data.table)
library(ggplot2)
library(readxl)
library(dplyr)
library(cowplot)

setwd("/gchm_lab/collab/marlis_pj/coloc/coloc_results/")

meta.data <- read_excel("../Autoimmune_GWAS_STING-seq.xlsx", sheet = "Completed_Meta")

files = system("ls *results.csv", intern = T)

# remove commas
merged = data.frame()

for(i in 1:length(files)){
  
  print(i)
  
  system(paste("sed 's/\\([0-9]\\+\\[b38\\][A-Z]\\),\\([A-Z]\\)/\\1|\\2/g' -i", files[i]))
  coloc_res <- fread(files[i], sep = ",")
  
  if(is.na(coloc_res[1,1]) == T){
    next("DF is empty")
  }
  
  file_name = gsub("_coloc_results.csv", "", files[i])
  
  coloc_res$file_name = file_name
  coloc_res$disease_phenotype = meta.data$disease_phenotype[meta.data$gcloud_name == file_name]
  
  merged = rbind(merged, coloc_res)
}

merged2 = merged %>% filter(is.na(nsnps) == F)

fwrite(merged2, "processed_coloc_results.txt", sep = "\t", row.names = F, quote = F)

# Summarize and order by n
plot_df <- merged2 %>% 
  group_by(disease_phenotype) %>%
  summarise(n = n()) %>%
  arrange(n)  # Ensure data is sorted

# Convert disease_phenotype into a factor with levels in ascending order
plot_df$disease_phenotype <- factor(plot_df$disease_phenotype, levels = plot_df$disease_phenotype)

# Save the plot
png("../plots/lollipop_plot.png", width = 12, height = 8, units = "in", res = 300)

ggplot(plot_df, aes(x=disease_phenotype, y=n)) +
  geom_segment(aes(x=disease_phenotype, xend=disease_phenotype, y=0, yend=n), color="skyblue") +
  geom_point(color="blue", size=4, alpha=0.6) +
  theme_cowplot() + ylab("number of colocalizations") +
  coord_flip() +  # Flip for horizontal bars
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

dev.off()
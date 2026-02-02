#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
#usage: Rscript library_batch.r your_file.csv

## Load the required libraries
library(dplyr)
library(stringr)

#Read metadata dataframe
meta  <-read.csv(args[1], sep=",")

subset_meta <- meta %>%
  filter(!is.na(scRNA_batch_Dave)) %>%
  filter(scRNAseq==TRUE)


# there is some spaces behind batch names in scRNA_batch_Dave
subset_meta$scRNA_batch_Dave <- str_trim(subset_meta$scRNA_batch_Dave, side = "right")

subset_meta$scRNA_batch_Dave <- str_trim(subset_meta$scRNA_batch_Dave, side = "left")

# Assuming 'subset_meta' is your data frame and 'scRNA_batch_Dave' is the column you want to use for subsetting
# Extract the specific column as a data frame
column_to_extract <- subset_meta[, "WGS_sampleID"]

# Get unique values from 'scRNA_batch_Dave'
unique_values <- unique(subset_meta$scRNA_batch_Dave)

# Create a list of data frames, each containing the extracted column for a unique value
subsets_list <- lapply(unique_values, function(x) subset(column_to_extract, subset_meta$scRNA_batch_Dave == x))


# Save each data frame in the list as a separate variable and write to a TXT file

for (i in 1:length(unique_values)) {
  subset_df <- subset(column_to_extract, subset_meta$scRNA_batch_Dave == unique_values[i])
  variable_name <- paste0("scGEX-", unique_values[i])
  assign(variable_name, subset_df)
  write.table(subset_df, file = paste0(variable_name, ".txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(tibble)
library(ggrepel)

in_dir <- "/gcgl/mmatos/cd4_aging_project/data/ATAC-seq/atac_preprocessing/STAR/Picard_Metrics"
files  <- list.files(in_dir, pattern = ".CollectMultipleMetrics.insert_size_metrics", full.names = TRUE)

# --- Parse ONLY the HISTOGRAM numbers into a df ---
extract_histogram_block <- function(lines) {
  i <- which(grepl("^##\\s*HISTOGRAM", lines))
  if (length(i) == 0) return(NULL)
  
  header_idx <- i[1] + 1
  data_start <- header_idx + 1
  if (data_start > length(lines)) return(NULL)
  
  rest <- lines[data_start:length(lines)]
  end_rel <- which(rest == "" | grepl("^##", rest) | grepl("^#", rest))
  data_end <- if (length(end_rel) == 0) length(rest) else (min(end_rel) - 1)
  if (data_end < 1) return(NULL)
  
  rest[1:data_end]
}

read_picard_hist_numbers <- function(path) {
  lines <- readLines(path, warn = FALSE)
  data_lines <- extract_histogram_block(lines)
  if (is.null(data_lines)) return(NULL)
  
  df <- read.table(
    text = paste(data_lines, collapse = "\n"),
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(df) < 2) return(NULL)
  
  colnames(df)[1:2] <- c("insert_size", "count")
  df$insert_size <- as.numeric(df$insert_size)
  df$count       <- as.numeric(df$count)
  
  df <- df %>% filter(!is.na(insert_size), !is.na(count), count > 0)
  if (nrow(df) == 0) return(NULL)
  
  sample_id <- str_match(basename(path), "^([^\\.]+)\\.")[,2]
  if (is.na(sample_id) || sample_id == "") sample_id <- basename(path)
  
  as_tibble(df) %>% mutate(sample = sample_id)
}

hist_df <- map(files, read_picard_hist_numbers) %>% compact() %>% bind_rows()

# --- “Density-like” y without converting histogram -> kernel density:
# normalize counts per sample (probability mass per bp; since insert_size is in 1-bp steps)
plot_df <- hist_df %>%
  group_by(sample) %>%
  mutate(y = count / sum(count)) %>%
  ungroup()

# label each curve near the right side (no legend)
label_df <- plot_df %>%
  group_by(sample) %>%
  filter(insert_size >= quantile(insert_size, 0.90, na.rm = TRUE)) %>%
  slice_max(order_by = y, n = 1, with_ties = FALSE) %>%
  ungroup()

size_dist<-ggplot(plot_df, aes(x = insert_size, y = y, group = sample)) +
  geom_line(linewidth = 0.1, alpha = 0.6) + theme_classic(base_size = 13) +
  xlim(c(0,1000))+
  labs(
    x = "Insert size (bp)",
    y = "Proportion of read pairs per 1-bp bin"  ) +
  theme(legend.position = "none",
        axis.text = element_text(size=16),
        axis.title = element_text(size=18))
  

ggsave("~/cd4_qtl_paper_figures/supplements/plots/atac_insert_size_density.pdf",size_dist, width = 6, height = 4)

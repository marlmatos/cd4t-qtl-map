library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(viridisLite)

tmp_root  <- "/gcgl/mmatos/cd4_aging_project/data/ATAC-seq/atac_preprocessing/nextflow_tmp"   # parent dir containing d2/0d4dad* etc
list_file <- "/gchm/ATAC-seq_analysis/ATACseq_nf-preprocessing/scripts/sleepy_murdock_locations.log"

# 1) Read list.txt and keep only Deeptools_QC_Scaled_Region jobs
jobs <- read.delim(list_file, header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE) %>%
  setNames(c("dir_stub","job_name","exit_code","status")) %>%
  filter(grepl("Deeptools_QC_Scaled_Region", job_name), status == "COMPLETED") %>%
  mutate(sample_from_job = str_match(job_name, "\\((T\\d+)\\)")[,2])

# 2) Expand wildcard dir stub and locate *.scaleregions.plotProfile.tab inside
find_profile_tab <- function(dir_stub, tmp_root) {
  matches <- Sys.glob(file.path(tmp_root, paste0(dir_stub, "*")))
  if (length(matches) == 0) return(character(0))
  
  tabs <- unlist(lapply(matches, function(d) {
    list.files(d, pattern = "\\.scaleregions\\.plotProfile\\.tab$", full.names = TRUE)
  }), use.names = FALSE)
  
  tabs
}

jobs <- jobs %>%
  mutate(tab_paths = map(dir_stub, ~find_profile_tab(.x, tmp_root))) %>%
  tidyr::unnest(tab_paths)

# 3) Read one deepTools plotProfile.tab into long df
read_plotProfile_tab <- function(path) {
  hdr <- readLines(path, n = 2, warn = FALSE)
  if (length(hdr) < 2) return(NULL)
  
  # x-axis bins from line 2 (keep numeric tokens)
  bin_tokens <- strsplit(hdr[2], "\t", fixed = TRUE)[[1]]
  x_bins <- suppressWarnings(as.numeric(bin_tokens))
  x_bins <- x_bins[!is.na(x_bins)]
  
  dat <- read.delim(path, header = FALSE, sep = "\t", quote = "", skip = 2,
                    stringsAsFactors = FALSE, fill = TRUE, check.names = FALSE)
  if (nrow(dat) == 0 || ncol(dat) < 3) return(NULL)
  
  keep <- !vapply(dat, function(col) all(is.na(col) | col == ""), logical(1))
  dat <- dat[, keep, drop = FALSE]
  if (ncol(dat) < 3) return(NULL)
  
  colnames(dat)[1:2] <- c("sample","group")
  
  ymat <- dat[, 3:ncol(dat), drop = FALSE]
  ymat[] <- lapply(ymat, function(z) suppressWarnings(as.numeric(z)))
  
  n_bins <- ncol(ymat)
  x <- if (length(x_bins) == n_bins) x_bins else seq_len(n_bins)
  
  bind_cols(dat[,1:2, drop = FALSE], as_tibble(ymat)) %>%
    mutate(.row = row_number()) %>%
    pivot_longer(cols = -(c(sample, group, .row)),
                 names_to = "bin_col", values_to = "signal") %>%
    group_by(.row) %>%
    mutate(bin = x) %>%
    ungroup() %>%
    select(sample, group, bin, signal) %>%
    filter(!is.na(signal))
}

plot_df <- jobs %>%
  mutate(df = map(tab_paths, read_plotProfile_tab)) %>%
  select(sample_from_job, tab_paths, df) %>%
  unnest(df) %>%
  mutate(sample = ifelse(is.na(sample) | sample == "", sample_from_job, sample))

# 4) Plot: all samples on one graph, no legend, label lines
label_df <- plot_df %>%
  group_by(sample) %>%
  filter(bin >= quantile(bin, 0.98, na.rm = TRUE)) %>%
  slice_max(order_by = signal, n = 1, with_ties = FALSE) %>%
  ungroup()

# If your deepTools run used: 3kb upstream + 5kb body + 3kb downstream (common),
# and you have N bins total, then:
# - upstream = 3/11 of bins
# - body     = 5/11 of bins
# - downstream = 3/11 of bins
# TSS is at the upstream/body boundary; TES at the body/downstream boundary.

n_bins <- length(sort(unique(plot_df$bin)))  # works whether bin is 1..N or 1.0..N

tss_pos <- n_bins * (3/11)
tes_pos <- n_bins * (8/11)

# (optional) snap to the nearest actual bin value
bin_vals <- sort(unique(plot_df$bin))
tss_pos <- bin_vals[which.min(abs(bin_vals - tss_pos))]
tes_pos <- bin_vals[which.min(abs(bin_vals - tes_pos))]

# Now plot with vertical lines + labeled axis breaks
profile.pt<-ggplot(plot_df, aes(x = bin, y = signal, group = sample)) +
  geom_line(linewidth = 0.1, alpha = 0.85) +
  geom_vline(xintercept = tss_pos, linetype = 2, linewidth = 0.4) +
  geom_vline(xintercept = tes_pos, linetype = 2, linewidth = 0.4) +
  scale_x_continuous(
    breaks = c(min(bin_vals), tss_pos, tes_pos, max(bin_vals)),
    labels = c("-3 kb", "TSS", "TES", "+3 kb")
  ) +
  labs(x = "Bins (scaled regions)", y = "Average signal") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.text = element_text(size=16),
        axis.title = element_text(size=18))

ggsave("~/cd4_qtl_paper_figures/supplements/plots/atac_insert_size_deeptools_scaled_regions.pdf",profile.pt, width = 6, height = 4)


library(dplyr)
library(readr)


# ---- MultiQC -> total reads per sample (FastQC) ----
FASTQC_DIR <- "/gcgl/mmatos/cd4_aging_project/data/scRNA-seq/results/scrnaseq_lib_preprocessing/fastQC/multiqc_data_1"

stopifnot(dir.exists(FASTQC_DIR))

# Helpers
to_numeric <- function(x) {
  # handles "1,234,567" or "1234567" or numeric already
  suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
}

pick_existing <- function(paths) {
  paths[file.exists(paths)]
}

# Candidate files
f_fastqc   <- file.path(FASTQC_DIR, "multiqc_fastqc.txt")
f_general  <- file.path(FASTQC_DIR, "multiqc_general_stats.txt")
f_json     <- file.path(FASTQC_DIR, "multiqc_data.json")

existing <- pick_existing(c(f_fastqc, f_general, f_json))
if (length(existing) == 0) {
  stop("No expected MultiQC files found in: ", FASTQC_DIR,
       "\nLooked for: multiqc_fastqc.txt, multiqc_general_stats.txt, multiqc_data.json")
}
message("Found: ", paste(basename(existing), collapse = ", "))

# 1) Prefer multiqc_fastqc.txt because it's directly FastQC metrics
if (file.exists(f_fastqc)) {
  fastqc <- read.delim(f_fastqc, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  
  # Find a sample column (MultiQC is usually 'Sample')
  sample_col <- intersect(c("Sample", "sample", "Sample Name", "sample_name"), colnames(fastqc))
  if (length(sample_col) == 0) {
    stop("Couldn't find a sample column in multiqc_fastqc.txt. Columns are:\n",
         paste(colnames(fastqc), collapse = ", "))
  }
  sample_col <- sample_col[1]
  
  # Find Total Sequences column
  total_col <- intersect(c("Total Sequences", "Total sequences", "total_sequences"), colnames(fastqc))
  if (length(total_col) == 0) {
    stop("Couldn't find 'Total Sequences' in multiqc_fastqc.txt. Columns are:\n",
         paste(colnames(fastqc), collapse = ", "))
  }
  total_col <- total_col[1]
  
  df <- data.frame(
    Sample = fastqc[[sample_col]],
    Total_Sequences = to_numeric(fastqc[[total_col]]),
    stringsAsFactors = FALSE
  )
  
} else if (file.exists(f_general)) {
  # 2) Fallback: multiqc_general_stats.txt sometimes includes the FastQC total reads columns
  gen <- read.delim(f_general, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  
  # Sample column is usually "Sample"
  if (!("Sample" %in% colnames(gen))) {
    stop("Couldn't find 'Sample' column in multiqc_general_stats.txt. Columns are:\n",
         paste(colnames(gen), collapse = ", "))
  }
  
  total_candidates <- c(
    "Total Sequences",
    "FastQC: Total Sequences",
    "FastQC: Total sequences",
    "FastQC_mqc-generalstats_fastqc-total_sequences"
  )
  total_col <- intersect(total_candidates, colnames(gen))
  
  if (length(total_col) == 0) {
    idx <- grep("FastQC.*Total.*Sequen", colnames(gen), ignore.case = TRUE)
    if (length(idx) == 0) {
      stop("Couldn't find a Total Sequences column in multiqc_general_stats.txt.\nColumns are:\n",
           paste(colnames(gen), collapse = ", "))
    }
    total_col <- colnames(gen)[idx[1]]
  } else {
    total_col <- total_col[1]
  }
  
  df <- data.frame(
    Sample = gen[["Sample"]],
    Total_Sequences = to_numeric(gen[[total_col]]),
    stringsAsFactors = FALSE
  )
  
} else {
  # 3) Last resort: parse multiqc_data.json (more annoying but doable without extra packages)
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Need jsonlite to parse multiqc_data.json. Install with install.packages('jsonlite').")
  }
  library(jsonlite)
  
  j <- jsonlite::fromJSON(f_json)
  
  # MultiQC JSON structure varies by version/modules.
  # j$report_saved_raw_data$fastqc OR j$report_data_sources / etc.
  # We'll search recursively for a data.frame/list that has "Total Sequences".
  find_total_sequences <- function(x) {
    out <- list()
    recurse <- function(obj, path = "") {
      if (is.data.frame(obj) && ("Total Sequences" %in% colnames(obj) || "Total sequences" %in% colnames(obj))) {
        out[[path]] <<- obj
      } else if (is.list(obj)) {
        for (nm in names(obj)) recurse(obj[[nm]], paste0(path, "/", nm))
      }
    }
    recurse(x, "")
    out
  }
  
  hits <- find_total_sequences(j)
  if (length(hits) == 0) stop("Parsed multiqc_data.json but couldn't find a table with 'Total Sequences'.")
  
  # Pick first hit
  tab <- hits[[1]]
  total_col <- if ("Total Sequences" %in% colnames(tab)) "Total Sequences" else "Total sequences"
  
  # Need sample names: often rownames are sample IDs
  sample <- rownames(tab)
  if (is.null(sample) || all(sample == "")) {
    # try a column called Sample
    if ("Sample" %in% colnames(tab)) sample <- tab[["Sample"]] else stop("Couldn't infer sample IDs from JSON hit.")
  }
  
  df <- data.frame(
    Sample = sample,
    Total_Sequences = to_numeric(tab[[total_col]]),
    stringsAsFactors = FALSE
  )
}

# Clean & sanity checks
df <- df[!is.na(df$Total_Sequences), ]
df <- df[order(df$Sample), ]

# Write per-sample table
out_tsv <- file.path(FASTQC_DIR, "multiqc_total_sequences.tsv")
write.table(df, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote: ", out_tsv)

# Summaries for paper
n_samples <- nrow(df)
total_reads <- sum(df$Total_Sequences, na.rm = TRUE)

cat("\n--- Summary ---\n")
cat("N samples: ", n_samples, "\n", sep = "")
cat("Total sequences (sum across samples): ", format(total_reads, big.mark = ","), "\n", sep = "")

# Only report pairs if you are confident about what your FastQC is counting.
cat("If FastQC Total Sequences are reads (common), total read pairs ≈ ",
    format(round(total_reads / 2), big.mark = ","), "\n", sep = "")


library(stringr)

per_sample_R2 <- fastqc %>%
  mutate(
    sample_id = str_replace(Sample, "_S\\d+_L\\d{3}_R2_001$", "")
  ) %>%
  group_by(sample_id) %>%
  summarise(
    n_lanes = n(),
    R2_total_sequences = sum(`Total Sequences`, na.rm = TRUE),
    mean_gc = weighted.mean(`%GC`, w = `Total Sequences`, na.rm = TRUE),
    mean_dedup_pct = weighted.mean(total_deduplicated_percentage, w = `Total Sequences`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(R2_total_sequences))

per_sample_R2


median(per_sample_R2$R2_total_sequences)
quantile(per_sample_R2$R2_total_sequences, c(0.25, 0.5, 0.75))

library(ggplot2)
library(scales)

# basic checks
stopifnot("R2_total_sequences" %in% colnames(per_sample_R2))

x <- per_sample_R2$R2_total_sequences
q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)

df <- data.frame(R2_total_sequences = x)

p_hist <- ggplot(df, aes(x = R2_total_sequences)) +
  geom_histogram(bins = 40, color = "grey20", linewidth = 0.25) +
  geom_vline(xintercept = q[2], linewidth = 0.8, linetype = 1) +   # median
  geom_vline(xintercept = q[c(1,3)], linewidth = 0.6, linetype = 2) + # IQR
  annotate(
    "text",
    x = q[2],
    y = Inf,
    vjust = 1.3,
    label = sprintf(
      "Median = %s\nIQR = %s–%s",
      scientific(q[2], digits = 3),
      scientific(q[1], digits = 3),
      scientific(q[3], digits = 3)
    ),
    size = 3.3
  ) +
  scale_x_continuous(labels = label_number(scale = 1e-6, suffix = "M")) +
  labs(
    x = "Read 2 total sequences per library (millions; summed across lanes)",
    y = "Number of libraries",
    title = "Sequencing depth distribution across scRNA-seq libraries (R2)"
  ) +
  theme_classic(base_size = 12)

ggsave("supp_R2_depth_histogram.pdf", p_hist, width = 6.5, height = 4.2, useDingbats = FALSE)
ggsave("supp_R2_depth_histogram.png", p_hist, width = 6.5, height = 4.2, dpi = 300)
p_hist


library(data.table)
library(stringr)

# ----------------------------
# Inputs
# ----------------------------
caqtl_dir <- "/gcgl/sghatan/marlis_pj/coloc/coloc_results/ca_gwas_coloc"
eqtl_dir  <- "/gcgl/sghatan/marlis_pj/coloc/coloc_results/eqtl_gwas_coloc"

out_csv   <- "~/cd4_qtl_paper_figures/figure_6/data/gwas_colocalization_study_ids.csv"

# trait_labels: named character vector where names are trait codes, values are pretty labels
# (you already have this in your environment)

# ----------------------------
# Helpers
# ----------------------------
has_data_rows <- function(f) {
  # TRUE if there is at least 1 row of data (beyond header); FALSE otherwise
  tryCatch({
    dt <- data.table::fread(f, nrows = 1, showProgress = FALSE)
    nrow(dt) > 0
  }, error = function(e) {
    FALSE
  })
}

list_csv_nonempty <- function(dir_path) {
  files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0) return(data.table(file = character(), full_path = character(), nonempty = logical()))
  data.table(
    file      = basename(files),
    full_path = files,
    nonempty  = vapply(files, has_data_rows, logical(1))
  )
}

# ----------------------------
# 1) List + filter non-empty files
# ----------------------------
caqtl_files <- list_csv_nonempty(caqtl_dir)[nonempty == TRUE]
eqtl_files  <- list_csv_nonempty(eqtl_dir)[nonempty == TRUE]

# (optional) see what got dropped
caqtl_dropped <- list_csv_nonempty(caqtl_dir)[nonempty == FALSE]
eqtl_dropped  <- list_csv_nonempty(eqtl_dir)[nonempty == FALSE]

# ----------------------------
# 2) Parse filenames into fields
# ----------------------------

# caQTL: CD4T_chromatin_<YEAR>_<ACCESSION>_<TRAIT>_<POP>_preprocessed...
m_caqtl <- str_match(
  caqtl_files$file,
  "^CD4T_chromatin_(\\d{4})_([^_]+)_([^_]+)_([^_]+)_preprocessed"
)
caqtl_std <- data.table(
  source    = "caqtl",
  file      = caqtl_files$file,
  full_path = caqtl_files$full_path,
  year      = m_caqtl[,2],
  accession = m_caqtl[,3],
  trait     = m_caqtl[,4],
  ethnicity = m_caqtl[,5]
)

# eQTL: <YEAR>_<ACCESSION>_<TRAIT>_<POP>_preprocessed...
m_eqtl <- str_match(
  eqtl_files$file,
  "^(\\d{4})_([^_]+)_([^_]+)_([^_]+)_preprocessed"
)
eqtl_std <- data.table(
  source    = "eqtl",
  file      = eqtl_files$file,
  full_path = eqtl_files$full_path,
  year      = m_eqtl[,2],
  accession = m_eqtl[,3],
  trait     = m_eqtl[,4],
  ethnicity = m_eqtl[,5]
)

# Keep only successfully parsed rows (prevents NA_NA_NA_NA)
caqtl_std <- caqtl_std[!is.na(year) & !is.na(accession) & !is.na(trait) & !is.na(ethnicity)]
eqtl_std  <- eqtl_std[!is.na(year) & !is.na(accession) & !is.na(trait) & !is.na(ethnicity)]

caqtl_std[, gwas_id := paste(year, accession, trait, ethnicity, sep = "_")]
eqtl_std[,  gwas_id := paste(year, accession, trait, ethnicity, sep = "_")]

# ----------------------------
# 3) UNION reference table (one row per gwas_id)
# ----------------------------
combined <- rbindlist(list(eqtl_std, caqtl_std), use.names = TRUE, fill = TRUE)

gwas_ref <- combined[, .(
  year       = year[1],
  accession  = accession[1],
  trait      = trait[1],
  ethnicity  = ethnicity[1],
  has_eqtl   = any(source == "eqtl"),
  has_caqtl  = any(source == "caqtl"),
  n_eqtl     = sum(source == "eqtl"),
  n_caqtl    = sum(source == "caqtl"),
  eqtl_files  = paste(file[source == "eqtl"],  collapse = ";"),
  caqtl_files = paste(file[source == "caqtl"], collapse = ";"),
  eqtl_paths  = paste(full_path[source == "eqtl"],  collapse = ";"),
  caqtl_paths = paste(full_path[source == "caqtl"], collapse = ";")
), by = gwas_id]

# ----------------------------
# 4) Add pretty trait labels
# ----------------------------
gwas_ref[, trait_label := unname(trait_labels[trait])]
gwas_ref[is.na(trait_label), trait_label := trait]

# ----------------------------
# 5) Save
# ----------------------------
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
fwrite(gwas_ref, file = out_csv)

# Optional: quick report
cat("Saved:", out_csv, "\n")
cat("Non-empty caQTL files:", nrow(caqtl_files), "\n")
cat("Non-empty eQTL files :", nrow(eqtl_files),  "\n")
cat("Union gwas_id rows   :", nrow(gwas_ref),    "\n")

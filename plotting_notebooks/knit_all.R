#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rmarkdown)
})

args <- commandArgs(trailingOnly = TRUE)

# Optional first arg: directory containing the figure notebooks.
# Defaults to the directory where this script lives.
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
script_dir <- dirname(script_path)
work_dir <- if (length(args) >= 1 && nzchar(args[1])) normalizePath(args[1], mustWork = TRUE) else script_dir

setwd(work_dir)

pattern <- "^plotting_notebook_figure_[0-9]+\\.Rmd$"
rmd_files <- list.files(path = work_dir, pattern = pattern, full.names = TRUE)

if (length(rmd_files) == 0) {
  stop("No figure Rmd files found in: ", work_dir)
}

# Render in numeric figure order.
figure_id <- as.integer(sub("^.*_figure_([0-9]+)\\.Rmd$", "\\1", basename(rmd_files)))
rmd_files <- rmd_files[order(figure_id)]

message("Found ", length(rmd_files), " notebooks to knit in ", work_dir)

results <- vector("list", length(rmd_files))

for (i in seq_along(rmd_files)) {
  input <- rmd_files[[i]]
  label <- basename(input)
  message("\n[", i, "/", length(rmd_files), "] Rendering ", label)

  ok <- TRUE
  out <- NA_character_
  err <- NULL

  tryCatch({
    out <- rmarkdown::render(
      input = input,
      quiet = FALSE,
      envir = new.env(parent = globalenv())
    )
  }, error = function(e) {
    ok <<- FALSE
    err <<- conditionMessage(e)
  })

  results[[i]] <- list(
    file = label,
    ok = ok,
    output = out,
    error = err
  )

  if (ok) {
    message("  OK: ", out)
  } else {
    message("  FAILED: ", err)
  }
}

ok_vec <- vapply(results, function(x) x$ok, logical(1))
n_ok <- sum(ok_vec)
n_fail <- length(results) - n_ok

message("\n========== Knit Summary ==========")
for (res in results) {
  status <- if (res$ok) "OK" else "FAILED"
  detail <- if (res$ok) res$output else res$error
  message("- ", res$file, " [", status, "] ", detail)
}
message("==================================")
message("Succeeded: ", n_ok, " | Failed: ", n_fail)

if (n_fail > 0) {
  quit(status = 1)
}

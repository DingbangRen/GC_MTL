#!/usr/bin/env Rscript

script_path <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]),
  winslash = "/",
  mustWork = FALSE
)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
source(file.path(repo_root, "R", "gcvs_repro.R"), local = globalenv())

parsed <- gcvs_parse_args(commandArgs(trailingOnly = TRUE))
report_id <- gcvs_default_if_null(parsed[["report-id"]], gcvs_default_if_null(parsed[["exp-id"]], ""))

if (!nzchar(report_id)) {
  stop("Missing required --report-id for report registration.")
}

as_index_string <- function(value){
  value <- gcvs_default_if_null(value, "")
  if (length(value) == 0 || is.na(value)) {
    return("")
  }
  as.character(value[[1]])
}

new_row <- data.frame(
  date = as_index_string(gcvs_default_if_null(parsed$date, format(Sys.Date(), "%Y-%m-%d"))),
  report_id = as_index_string(report_id),
  report_type = tolower(as_index_string(gcvs_default_if_null(parsed[["report-type"]], "daily"))),
  experiment_id = as_index_string(gcvs_default_if_null(parsed[["experiment-id"]], gcvs_default_if_null(parsed[["exp-id"]], ""))),
  dataset = as_index_string(gcvs_default_if_null(parsed$dataset, "")),
  profile = as_index_string(gcvs_default_if_null(parsed$profile, "")),
  M1_rule = as_index_string(gcvs_default_if_null(parsed[["m1-rule"]], "")),
  M2_rule = as_index_string(gcvs_default_if_null(parsed[["m2-rule"]], "")),
  theta_mean = as_index_string(gcvs_default_if_null(parsed[["theta-mean"]], "")),
  theta_peaks = as_index_string(gcvs_default_if_null(parsed[["theta-peaks"]], "")),
  theta_ess = as_index_string(gcvs_default_if_null(parsed[["theta-ess"]], "")),
  M_mean = as_index_string(gcvs_default_if_null(parsed[["m-mean"]], "")),
  decision = as_index_string(gcvs_default_if_null(parsed$decision, "INCONCLUSIVE")),
  report_path = gcvs_repo_relative_path(repo_root, as_index_string(gcvs_default_if_null(parsed[["report-path"]], ""))),
  html_path = gcvs_repo_relative_path(repo_root, as_index_string(gcvs_default_if_null(parsed[["html-path"]], ""))),
  pdf_path = gcvs_repo_relative_path(repo_root, as_index_string(gcvs_default_if_null(parsed[["pdf-path"]], ""))),
  meta_path = gcvs_repo_relative_path(repo_root, as_index_string(gcvs_default_if_null(parsed[["meta-path"]], ""))),
  stringsAsFactors = FALSE
)

report_index_csv <- file.path(repo_root, "results", "reports", "index.csv")
gcvs_register_report_entry(report_index_csv, new_row)

cat("Registered report:", new_row$report_id[[1]], "\n")
cat("Report type:", new_row$report_type[[1]], "\n")
cat("Report index:", gcvs_repo_relative_path(repo_root, report_index_csv), "\n")

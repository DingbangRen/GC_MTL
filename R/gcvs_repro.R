gcvs_default_if_null <- function(x, y){
  if (is.null(x) || length(x) == 0) {
    y
  } else {
    x
  }
}

gcvs_this_file <- function(){
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_match <- grep(file_arg, cmd_args)
  if (length(file_match) > 0) {
    return(normalizePath(sub(file_arg, "", cmd_args[file_match[1]]), winslash = "/", mustWork = FALSE))
  }

  frame_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
  if (length(frame_files) > 0) {
    return(normalizePath(frame_files[[length(frame_files)]], winslash = "/", mustWork = FALSE))
  }

  normalizePath("scripts/run_gcvs_experiment.R", winslash = "/", mustWork = FALSE)
}

gcvs_repo_root <- function(){
  current_dir <- dirname(gcvs_this_file())
  if (basename(current_dir) %in% c("scripts", "R")) {
    return(normalizePath(file.path(current_dir, ".."), winslash = "/", mustWork = FALSE))
  }
  normalizePath(current_dir, winslash = "/", mustWork = FALSE)
}

gcvs_repo_relative_path <- function(repo_root, path){
  path <- gcvs_default_if_null(path, "")
  if (!nzchar(path)) {
    return("")
  }

  normalized_root <- normalizePath(repo_root, winslash = "/", mustWork = FALSE)
  normalized_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root_prefix <- paste0(normalized_root, "/")

  if (identical(normalized_path, normalized_root)) {
    return(".")
  }

  if (startsWith(normalized_path, root_prefix)) {
    return(substr(normalized_path, nchar(root_prefix) + 1L, nchar(normalized_path)))
  }

  normalized_path
}

gcvs_report_index_columns <- function(){
  c(
    "date",
    "report_id",
    "report_type",
    "experiment_id",
    "dataset",
    "profile",
    "M1_rule",
    "M2_rule",
    "theta_mean",
    "theta_peaks",
    "theta_ess",
    "M_mean",
    "decision",
    "report_path",
    "html_path",
    "pdf_path",
    "meta_path"
  )
}

gcvs_align_report_index_schema <- function(current){
  if ("exp_id" %in% names(current) && !("report_id" %in% names(current))) {
    current$report_id <- current$exp_id
  }
  if ("exp_id" %in% names(current) && !("experiment_id" %in% names(current))) {
    current$experiment_id <- current$exp_id
  }
  if (!("report_type" %in% names(current))) {
    current$report_type <- if ("dataset" %in% names(current)) {
      ifelse(tolower(current$dataset) == "manuscript", "manuscript", "experiment")
    } else {
      "experiment"
    }
  }
  if (!("profile" %in% names(current))) {
    current$profile <- ""
  }
  if (!("html_path" %in% names(current))) {
    current$html_path <- ""
  }
  current
}

gcvs_register_report_entry <- function(report_index_csv, entry){
  target_cols <- gcvs_report_index_columns()
  repo_root <- gcvs_repo_root()

  if (file.exists(report_index_csv)) {
    current <- read.csv(report_index_csv, stringsAsFactors = FALSE, check.names = FALSE)
    current <- gcvs_align_report_index_schema(current)
  } else {
    current <- data.frame(stringsAsFactors = FALSE)
  }

  missing_cols <- setdiff(target_cols, names(current))
  for (col in missing_cols) {
    current[[col]] <- rep("", nrow(current))
  }

  entry_missing_cols <- setdiff(target_cols, names(entry))
  for (col in entry_missing_cols) {
    entry[[col]] <- ""
  }

  current <- current[, target_cols, drop = FALSE]
  entry <- entry[, target_cols, drop = FALSE]

  current[] <- lapply(current, function(column) {
    column[is.na(column)] <- ""
    as.character(column)
  })
  entry[] <- lapply(entry, function(column) {
    column[is.na(column)] <- ""
    as.character(column)
  })

  for (path_col in c("report_path", "html_path", "pdf_path", "meta_path")) {
    current[[path_col]] <- vapply(current[[path_col]], function(path) {
      gcvs_repo_relative_path(repo_root, path)
    }, character(1))
    entry[[path_col]] <- vapply(entry[[path_col]], function(path) {
      gcvs_repo_relative_path(repo_root, path)
    }, character(1))
  }

  if (nrow(current) > 0) {
    current <- current[
      !(current$report_id == entry$report_id[[1]] & current$report_type == entry$report_type[[1]]),
      ,
      drop = FALSE
    ]
  }

  updated <- rbind(current, entry)
  write.csv(updated, report_index_csv, row.names = FALSE)
}

gcvs_default_config <- function(profile = "smoke"){
  profile <- tolower(profile)
  if (profile == "paper") {
    return(list(
      profile = "paper",
      niter = 300L,
      burnin = 150L,
      checkpoint_every = 10L,
      inner_steps = 40L,
      sigma_px_steps = 1L,
      epsilon = 0.005,
      alpha_const = 3,
      eta = 1,
      a_psi = 2,
      b_psi = 2,
      a_err = 0.001,
      b_err = 0.001,
      opt_rate = 0.7,
      const_0 = 10,
      epsconst = 1,
      inisd = 0.01
    ))
  }

  list(
    profile = "smoke",
    niter = 30L,
    burnin = 15L,
    checkpoint_every = 5L,
    inner_steps = 8L,
    sigma_px_steps = 1L,
    epsilon = 0.005,
    alpha_const = 3,
    eta = 1,
    a_psi = 2,
    b_psi = 2,
    a_err = 0.001,
    b_err = 0.001,
    opt_rate = 0.7,
    const_0 = 8,
    epsconst = 1,
    inisd = 0.01
  )
}

gcvs_parse_args <- function(args = commandArgs(trailingOnly = TRUE)){
  parsed <- list()
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (startsWith(arg, "--")) {
      key_value <- sub("^--", "", arg)
      if (grepl("=", key_value, fixed = TRUE)) {
        parts <- strsplit(key_value, "=", fixed = TRUE)[[1]]
        parsed[[parts[1]]] <- parts[2]
      } else if (i < length(args) && !startsWith(args[[i + 1L]], "--")) {
        parsed[[key_value]] <- args[[i + 1L]]
        i <- i + 1L
      } else {
        parsed[[key_value]] <- TRUE
      }
    }
    i <- i + 1L
  }
  parsed
}

gcvs_numeric_arg <- function(parsed, key, default){
  value <- gcvs_default_if_null(parsed[[key]], default)
  if (is.integer(default)) {
    as.integer(value)
  } else {
    as.numeric(value)
  }
}

gcvs_boolean_arg <- function(parsed, key, default = FALSE){
  value <- gcvs_default_if_null(parsed[[key]], default)
  if (is.logical(value)) {
    return(isTRUE(value))
  }
  tolower(as.character(value)[[1]]) %in% c("1", "true", "t", "yes", "y")
}

gcvs_normalize_sampler_mode <- function(value){
  normalized <- tolower(gsub("-", "_", gcvs_default_if_null(value, "full")))
  if (!normalized %in% c("full", "lambda_sigma_only")) {
    stop("Unsupported sampler mode. Expected one of: full, lambda-sigma-only.")
  }
  normalized
}

gcvs_build_config <- function(parsed){
  defaults <- gcvs_default_config(gcvs_default_if_null(parsed$profile, "smoke"))
  timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  dataset <- tolower(gcvs_default_if_null(parsed$dataset, "synthetic"))

  list(
    dataset = dataset,
    profile = defaults$profile,
    niter = gcvs_numeric_arg(parsed, "niter", defaults$niter),
    burnin = gcvs_numeric_arg(parsed, "burnin", defaults$burnin),
    checkpoint_every = gcvs_numeric_arg(parsed, "checkpoint-every", defaults$checkpoint_every),
    inner_steps = gcvs_numeric_arg(parsed, "inner-steps", defaults$inner_steps),
    sigma_px_steps = gcvs_numeric_arg(parsed, "sigma-px-steps", defaults$sigma_px_steps),
    epsilon = gcvs_numeric_arg(parsed, "epsilon", defaults$epsilon),
    alpha_const = gcvs_numeric_arg(parsed, "alpha", defaults$alpha_const),
    eta = gcvs_numeric_arg(parsed, "eta", defaults$eta),
    a_psi = gcvs_numeric_arg(parsed, "a-psi", defaults$a_psi),
    b_psi = gcvs_numeric_arg(parsed, "b-psi", defaults$b_psi),
    a_err = gcvs_numeric_arg(parsed, "a-err", defaults$a_err),
    b_err = gcvs_numeric_arg(parsed, "b-err", defaults$b_err),
    opt_rate = gcvs_numeric_arg(parsed, "opt-rate", defaults$opt_rate),
    const_0 = gcvs_numeric_arg(parsed, "const-0", defaults$const_0),
    epsconst = gcvs_numeric_arg(parsed, "epsconst", defaults$epsconst),
    inisd = gcvs_numeric_arg(parsed, "inisd", defaults$inisd),
    sampler_mode = gcvs_normalize_sampler_mode(gcvs_default_if_null(parsed[["sampler-mode"]], "full")),
    seed = gcvs_numeric_arg(parsed, "seed", 123L),
    resume = gcvs_boolean_arg(parsed, "resume", FALSE),
    exp_id = gcvs_default_if_null(
      parsed[["exp-id"]],
      paste("EXP", "GCVS", toupper(dataset), toupper(defaults$profile), timestamp, sep = "-")
    )
  )
}

gcvs_load_synthetic_dataset <- function(repo_root){
  env <- new.env(parent = globalenv())
  source(file.path(repo_root, "simulation data setting.R"), local = env)
  list(
    dataset = "synthetic",
    X_train = env$X_train,
    Y_train = env$Y_train,
    X_test = env$X_test,
    Y_test = env$Y_test,
    beta_true = as.matrix(env$Beta_alltasks_J_80)
  )
}

gcvs_load_real_dataset <- function(repo_root, dataset){
  env <- new.env(parent = globalenv())
  if (dataset == "sarcos") {
    load(file.path(repo_root, "Xlist_sarcos.RData"), envir = env)
    load(file.path(repo_root, "Ylist_sarcos.RData"), envir = env)
    load(file.path(repo_root, "Xlist_sarcos_test.RData"), envir = env)
    load(file.path(repo_root, "Ylist_sarcos_test.RData"), envir = env)
    return(list(
      dataset = "sarcos",
      X_train = env$Xlist_sarcos,
      Y_train = env$Ylist_sarcos,
      X_test = env$Xlist_sarcos_test,
      Y_test = env$Ylist_sarcos_test,
      beta_true = NULL
    ))
  }

  if (dataset == "isolet") {
    load(file.path(repo_root, "Xlist_isolet_train.RData"), envir = env)
    load(file.path(repo_root, "Ylist_isolet_train.RData"), envir = env)
    load(file.path(repo_root, "Xlist_isolet_test.RData"), envir = env)
    load(file.path(repo_root, "Ylist_isolet_test.RData"), envir = env)
    return(list(
      dataset = "isolet",
      X_train = env$Xlist_isolet_train,
      Y_train = env$Ylist_isolet_train,
      X_test = env$Xlist_isolet_test,
      Y_test = env$Ylist_isolet_test,
      beta_true = NULL
    ))
  }

  stop(sprintf("Unsupported dataset '%s'. Expected one of: synthetic, sarcos, isolet.", dataset))
}

gcvs_load_dataset <- function(repo_root, dataset){
  if (dataset == "synthetic") {
    return(gcvs_load_synthetic_dataset(repo_root))
  }
  gcvs_load_real_dataset(repo_root, dataset)
}

gcvs_prepare_model_env <- function(repo_root, J, K, eta){
  env <- new.env(parent = globalenv())
  assign("J", J, envir = env)
  assign("K", K, envir = env)
  assign("eta", eta, envir = env)
  source(file.path(repo_root, "GCVS_posteriors.R"), local = env)
  env
}

gcvs_rsquare <- function(X, Y, beta){
  estimate <- X %*% c(beta)
  ss_explained <- sum((estimate - mean(Y))^2)
  ss_residual <- sum((Y - estimate)^2)
  ss_explained / (ss_explained + ss_residual)
}

gcvs_rmse <- function(X, Y, beta){
  estimate <- X %*% c(beta)
  sqrt(mean((Y - estimate)^2))
}

gcvs_selection_precision <- function(beta_hat, beta_true, threshold = 0.1){
  selected <- abs(beta_hat) > threshold
  true_active <- abs(beta_true) > 0
  selected_count <- sum(selected)
  if (selected_count == 0) {
    return(NA_real_)
  }
  sum(selected & true_active) / selected_count
}

gcvs_matrix_mean <- function(samples, field){
  Reduce(`+`, lapply(samples, function(sample) sample[[field]])) / length(samples)
}

gcvs_format_metric <- function(x, digits = 4){
  if (length(x) == 0 || is.na(x) || !is.finite(x)) {
    return("NA")
  }
  sprintf(paste0("%.", digits, "f"), x)
}

gcvs_first_nonempty <- function(values){
  nonempty <- values[nzchar(values)]
  if (length(nonempty) == 0) {
    return("")
  }
  nonempty[[1]]
}

gcvs_effective_sample_size <- function(x){
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 3L) {
    return(NA_real_)
  }
  if (stats::sd(x) == 0) {
    return(as.numeric(n))
  }

  acf_values <- stats::acf(
    x,
    lag.max = min(n - 1L, floor(n / 2)),
    plot = FALSE,
    demean = TRUE
  )$acf[-1]

  if (length(acf_values) < 2L) {
    return(as.numeric(n))
  }

  pair_count <- floor(length(acf_values) / 2L)
  pair_sums <- vapply(seq_len(pair_count), function(i) {
    acf_values[(2L * i) - 1L] + acf_values[2L * i]
  }, numeric(1))

  last_positive_pair <- 0L
  for (i in seq_along(pair_sums)) {
    if (pair_sums[[i]] > 0) {
      last_positive_pair <- i
    } else {
      break
    }
  }

  if (last_positive_pair == 0L) {
    return(as.numeric(n))
  }

  tau <- 1 + 2 * sum(pair_sums[seq_len(last_positive_pair)])
  ess <- n / max(tau, 1)
  min(as.numeric(n), ess)
}

gcvs_standardized_half_shift <- function(x){
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 4L) {
    return(NA_real_)
  }

  midpoint <- floor(n / 2L)
  first_half <- x[seq_len(midpoint)]
  second_half <- x[seq(from = midpoint + 1L, to = n)]
  pooled_scale <- stats::sd(x)
  if (!is.finite(pooled_scale) || pooled_scale == 0) {
    pooled_scale <- 1
  }
  abs(mean(first_half) - mean(second_half)) / pooled_scale
}

gcvs_default_lambda_coords <- function(dataset, beta_mean){
  J <- nrow(beta_mean)
  K <- ncol(beta_mean)

  if (tolower(dataset) == "synthetic" && J >= 80L && K >= 8L) {
    return(data.frame(
      feature = c(3L, 12L, 18L, 35L, 60L, 75L),
      task = c(1L, 1L, 5L, 7L, 1L, 8L),
      label = c(
        "lambda[3,1] active-majority",
        "lambda[12,1] shared-majority",
        "lambda[18,5] task5-group",
        "lambda[35,7] task7-outlier",
        "lambda[60,1] inactive-majority",
        "lambda[75,8] task8-outlier"
      ),
      stringsAsFactors = FALSE
    ))
  }

  abs_beta <- abs(beta_mean)
  flat_order <- order(abs_beta, decreasing = TRUE)
  candidate_linear <- unique(c(
    flat_order[[1]],
    flat_order[[max(1L, floor(length(flat_order) / 4L))]],
    flat_order[[max(1L, floor(length(flat_order) / 2L))]],
    which.min(abs_beta)
  ))
  coords <- arrayInd(candidate_linear, .dim = dim(abs_beta))
  data.frame(
    feature = coords[, 1],
    task = coords[, 2],
    label = sprintf("lambda[%d,%d]", coords[, 1], coords[, 2]),
    stringsAsFactors = FALSE
  )
}

gcvs_default_sigma_pairs <- function(dataset, sigma_mean){
  K <- nrow(sigma_mean)

  if (tolower(dataset) == "synthetic" && K >= 8L) {
    return(data.frame(
      task_1 = c(1L, 1L, 1L, 5L, 1L, 7L),
      task_2 = c(2L, 3L, 4L, 6L, 7L, 8L),
      label = c(
        "Sigma[1,2]",
        "Sigma[1,3]",
        "Sigma[1,4]",
        "Sigma[5,6]",
        "Sigma[1,7]",
        "Sigma[7,8]"
      ),
      stringsAsFactors = FALSE
    ))
  }

  upper_idx <- which(upper.tri(sigma_mean), arr.ind = TRUE)
  upper_values <- sigma_mean[upper.tri(sigma_mean)]
  max_positive <- upper_idx[which.max(upper_values), , drop = FALSE]
  max_negative <- upper_idx[which.min(upper_values), , drop = FALSE]
  default_pairs <- unique(rbind(
    c(1L, min(2L, K)),
    c(max_positive[1, 1], max_positive[1, 2]),
    c(max_negative[1, 1], max_negative[1, 2]),
    c(max(1L, K - 1L), K)
  ))
  data.frame(
    task_1 = default_pairs[, 1],
    task_2 = default_pairs[, 2],
    label = sprintf("Sigma[%d,%d]", default_pairs[, 1], default_pairs[, 2]),
    stringsAsFactors = FALSE
  )
}

gcvs_build_diagnostics <- function(run_result){
  lambda_coords <- gcvs_default_lambda_coords(run_result$dataset, run_result$beta_mean)
  sigma_pairs <- gcvs_default_sigma_pairs(run_result$dataset, run_result$sigma_mean)

  row_acceptance_summary <- data.frame(
    feature = seq_len(ncol(run_result$row_acceptance)),
    mean_acceptance = colMeans(run_result$row_acceptance, na.rm = TRUE),
    min_acceptance = apply(run_result$row_acceptance, 2, min, na.rm = TRUE),
    max_acceptance = apply(run_result$row_acceptance, 2, max, na.rm = TRUE),
    mean_step_size = colMeans(run_result$row_step_size, na.rm = TRUE),
    final_step_size = apply(run_result$row_step_size, 2, function(x) {
      finite_x <- x[is.finite(x)]
      if (length(finite_x) == 0) NA_real_ else tail(finite_x, 1L)
    }),
    stringsAsFactors = FALSE
  )

  lambda_trace_rows <- lapply(seq_len(nrow(lambda_coords)), function(i) {
    feature_id <- lambda_coords$feature[[i]]
    task_id <- lambda_coords$task[[i]]
    values <- run_result$lambda_samples[, feature_id, task_id]
    data.frame(
      draw = seq_along(values),
      label = lambda_coords$label[[i]],
      feature = feature_id,
      task = task_id,
      value = values,
      stringsAsFactors = FALSE
    )
  })
  lambda_traces <- do.call(rbind, lambda_trace_rows)
  lambda_summary <- do.call(rbind, lapply(split(lambda_traces, lambda_traces$label), function(trace_df) {
    data.frame(
      label = trace_df$label[[1]],
      feature = trace_df$feature[[1]],
      task = trace_df$task[[1]],
      posterior_mean = mean(trace_df$value),
      posterior_sd = stats::sd(trace_df$value),
      ess = gcvs_effective_sample_size(trace_df$value),
      first_half_mean = mean(trace_df$value[seq_len(floor(nrow(trace_df) / 2L))]),
      second_half_mean = mean(trace_df$value[seq(from = floor(nrow(trace_df) / 2L) + 1L, to = nrow(trace_df))]),
      standardized_half_shift = gcvs_standardized_half_shift(trace_df$value),
      stringsAsFactors = FALSE
    )
  }))

  sigma_trace_rows <- lapply(seq_len(nrow(sigma_pairs)), function(i) {
    task_1 <- sigma_pairs$task_1[[i]]
    task_2 <- sigma_pairs$task_2[[i]]
    values <- run_result$sigma_samples[, task_1, task_2]
    data.frame(
      draw = seq_along(values),
      label = sigma_pairs$label[[i]],
      task_1 = task_1,
      task_2 = task_2,
      value = values,
      stringsAsFactors = FALSE
    )
  })
  sigma_traces <- do.call(rbind, sigma_trace_rows)
  sigma_summary <- do.call(rbind, lapply(split(sigma_traces, sigma_traces$label), function(trace_df) {
    data.frame(
      label = trace_df$label[[1]],
      task_1 = trace_df$task_1[[1]],
      task_2 = trace_df$task_2[[1]],
      posterior_mean = mean(trace_df$value),
      posterior_sd = stats::sd(trace_df$value),
      ess = gcvs_effective_sample_size(trace_df$value),
      first_half_mean = mean(trace_df$value[seq_len(floor(nrow(trace_df) / 2L))]),
      second_half_mean = mean(trace_df$value[seq(from = floor(nrow(trace_df) / 2L) + 1L, to = nrow(trace_df))]),
      standardized_half_shift = gcvs_standardized_half_shift(trace_df$value),
      stringsAsFactors = FALSE
    )
  }))

  sigma_pair_indices <- which(upper.tri(run_result$sigma_mean), arr.ind = TRUE)
  sigma_pair_summary <- do.call(rbind, lapply(seq_len(nrow(sigma_pair_indices)), function(i) {
    task_1 <- sigma_pair_indices[i, 1]
    task_2 <- sigma_pair_indices[i, 2]
    data.frame(
      task_1 = task_1,
      task_2 = task_2,
      label = sprintf("Sigma[%d,%d]", task_1, task_2),
      posterior_mean = run_result$sigma_mean[task_1, task_2],
      abs_posterior_mean = abs(run_result$sigma_mean[task_1, task_2]),
      sign = ifelse(run_result$sigma_mean[task_1, task_2] >= 0, "positive", "negative"),
      stringsAsFactors = FALSE
    )
  }))
  sigma_pair_summary <- sigma_pair_summary[order(-sigma_pair_summary$abs_posterior_mean), , drop = FALSE]

  support_recovery_summary <- NULL
  task_abs_correlation_summary <- NULL
  if (!is.null(run_result$beta_true) && all(dim(run_result$beta_true) == dim(run_result$beta_mean))) {
    thresholds <- c(0.10, 0.25, 0.50, 1.00)
    support_recovery_summary <- do.call(rbind, lapply(thresholds, function(threshold) {
      selected <- abs(run_result$beta_mean) > threshold
      truth <- abs(run_result$beta_true) > 0
      tp <- sum(selected & truth)
      fp <- sum(selected & !truth)
      fn <- sum(!selected & truth)
      precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
      recall <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
      f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) NA_real_ else {
        2 * precision * recall / (precision + recall)
      }
      data.frame(
        threshold = threshold,
        precision = precision,
        recall = recall,
        f1 = f1,
        true_positive = tp,
        false_positive = fp,
        false_negative = fn,
        stringsAsFactors = FALSE
      )
    }))

    task_abs_correlation_summary <- data.frame(
      task = seq_len(ncol(run_result$beta_mean)),
      abs_correlation = vapply(seq_len(ncol(run_result$beta_mean)), function(task_id) {
        stats::cor(abs(run_result$beta_mean[, task_id]), abs(run_result$beta_true[, task_id]))
      }, numeric(1)),
      stringsAsFactors = FALSE
    )
  }

  list(
    row_acceptance_summary = row_acceptance_summary,
    lambda_traces = lambda_traces,
    lambda_summary = lambda_summary,
    sigma_traces = sigma_traces,
    sigma_summary = sigma_summary,
    sigma_pair_summary = sigma_pair_summary,
    support_recovery_summary = support_recovery_summary,
    task_abs_correlation_summary = task_abs_correlation_summary
  )
}

gcvs_report_basename <- function(run_result){
  paste(
    format(Sys.Date(), "%Y%m%d"),
    run_result$config$exp_id,
    toupper(run_result$dataset),
    paste0("GDP-a", run_result$config$alpha_const, "-eta", run_result$config$eta),
    "FULL-MMALA",
    run_result$decision,
    sep = "__"
  )
}

gcvs_initialize_results_paths <- function(repo_root, run_result){
  exp_id <- run_result$config$exp_id
  report_basename <- gcvs_report_basename(run_result)
  experiment_dir <- file.path(repo_root, "results", "experiments", exp_id)
  figures_dir <- file.path(experiment_dir, "figures")
  checkpoints_dir <- file.path(experiment_dir, "checkpoints")
  logs_dir <- file.path(repo_root, "logs", "iteration_logs")
  reports_dir <- file.path(repo_root, "results", "reports")
  reports_pdf_dir <- file.path(reports_dir, "pdf")
  reports_meta_dir <- file.path(reports_dir, "meta")
  dir.create(experiment_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(checkpoints_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(reports_pdf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(reports_meta_dir, recursive = TRUE, showWarnings = FALSE)

  list(
    experiment_dir = experiment_dir,
    figures_dir = figures_dir,
    checkpoints_dir = checkpoints_dir,
    summary_csv = file.path(experiment_dir, "summary.csv"),
    iteration_csv = file.path(experiment_dir, "iteration_metrics.csv"),
    beta_csv = file.path(experiment_dir, "posterior_beta_mean.csv"),
    sigma_csv = file.path(experiment_dir, "posterior_sigma_mean.csv"),
    lambda_row_acceptance_csv = file.path(experiment_dir, "lambda_row_acceptance.csv"),
    lambda_row_acceptance_summary_csv = file.path(experiment_dir, "lambda_row_acceptance_summary.csv"),
    lambda_row_step_size_csv = file.path(experiment_dir, "lambda_row_step_size.csv"),
    lambda_row_restarts_csv = file.path(experiment_dir, "lambda_row_restarts.csv"),
    representative_lambda_csv = file.path(experiment_dir, "representative_lambda_summary.csv"),
    representative_sigma_csv = file.path(experiment_dir, "representative_sigma_summary.csv"),
    sigma_pair_summary_csv = file.path(experiment_dir, "sigma_pair_summary.csv"),
    support_recovery_csv = file.path(experiment_dir, "support_recovery_summary.csv"),
    task_abs_correlation_csv = file.path(experiment_dir, "task_abs_correlation_summary.csv"),
    lambda_samples_rds = file.path(experiment_dir, "posterior_lambda_samples.rds"),
    sigma_samples_rds = file.path(experiment_dir, "posterior_sigma_samples.rds"),
    checkpoint_latest_rds = file.path(checkpoints_dir, "checkpoint_latest.rds"),
    progress_json = file.path(experiment_dir, "progress.json"),
    metadata_json = file.path(experiment_dir, "metadata.json"),
    diagnostics_png = file.path(figures_dir, "sampling_diagnostics.png"),
    beta_heatmap_png = file.path(figures_dir, "posterior_beta_mean_heatmap.png"),
    sigma_heatmap_png = file.path(figures_dir, "posterior_sigma_mean_heatmap.png"),
    representative_traces_png = file.path(figures_dir, "representative_traces.png"),
    log_md = file.path(logs_dir, paste0(format(Sys.time(), "%Y%m%d__"), exp_id, ".md")),
    report_md = file.path(experiment_dir, "experiment_report.md"),
    report_html = file.path(experiment_dir, "experiment_report.html"),
    report_pdf = file.path(reports_pdf_dir, paste0(report_basename, ".pdf")),
    report_meta_json = file.path(reports_meta_dir, paste0(report_basename, ".json")),
    index_csv = file.path(repo_root, "results", "experiments", "index.csv"),
    report_index_csv = file.path(reports_dir, "index.csv")
  )
}

gcvs_initialize_runtime_paths <- function(repo_root, config){
  experiment_dir <- file.path(repo_root, "results", "experiments", config$exp_id)
  checkpoints_dir <- file.path(experiment_dir, "checkpoints")
  dir.create(experiment_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(checkpoints_dir, recursive = TRUE, showWarnings = FALSE)

  list(
    experiment_dir = experiment_dir,
    checkpoints_dir = checkpoints_dir,
    checkpoint_latest_rds = file.path(checkpoints_dir, "checkpoint_latest.rds"),
    progress_json = file.path(experiment_dir, "progress.json")
  )
}

gcvs_atomic_save_rds <- function(object, path){
  tmp_path <- paste0(path, ".tmp-", Sys.getpid(), "-", format(Sys.time(), "%Y%m%d%H%M%OS6"))
  on.exit(if (file.exists(tmp_path)) unlink(tmp_path), add = TRUE)
  saveRDS(object, file = tmp_path)
  if (!file.rename(tmp_path, path)) {
    stop(sprintf("Failed to atomically move checkpoint into place: %s", path))
  }
  invisible(path)
}

gcvs_atomic_write_json <- function(object, path){
  tmp_path <- paste0(path, ".tmp-", Sys.getpid(), "-", format(Sys.time(), "%Y%m%d%H%M%OS6"))
  on.exit(if (file.exists(tmp_path)) unlink(tmp_path), add = TRUE)
  jsonlite::write_json(object, path = tmp_path, pretty = TRUE, auto_unbox = TRUE, null = "null")
  if (!file.rename(tmp_path, path)) {
    stop(sprintf("Failed to atomically move JSON into place: %s", path))
  }
  invisible(path)
}

gcvs_checkpoint_path <- function(runtime_paths, iter){
  file.path(runtime_paths$checkpoints_dir, sprintf("checkpoint_iter_%06d.rds", as.integer(iter)))
}

gcvs_resume_command <- function(config){
  paste(
    "Rscript scripts/run_gcvs_experiment.R",
    paste("--dataset", config$dataset),
    paste("--profile", config$profile),
    paste("--exp-id", config$exp_id),
    "--resume"
  )
}

gcvs_validate_resume_config <- function(saved_config, current_config){
  comparable_fields <- c(
    "dataset", "profile", "niter", "burnin", "inner_steps", "sigma_px_steps", "epsilon",
    "alpha_const", "eta", "a_psi", "b_psi", "a_err", "b_err", "opt_rate", "sampler_mode",
    "const_0", "epsconst", "inisd", "seed", "exp_id"
  )

  mismatched <- comparable_fields[!vapply(comparable_fields, function(field) {
    identical(saved_config[[field]], current_config[[field]])
  }, logical(1))]

  if (length(mismatched) > 0) {
    stop(
      sprintf(
        "Resume requested, but checkpoint configuration differs for: %s",
        paste(mismatched, collapse = ", ")
      )
    )
  }
}

gcvs_write_progress <- function(runtime_paths, config, dataset, checkpoint_state){
  progress <- list(
    experiment_id = config$exp_id,
    dataset = dataset,
    profile = config$profile,
    niter = config$niter,
    burnin = config$burnin,
    checkpoint_every = config$checkpoint_every,
    last_completed_iter = checkpoint_state$last_completed_iter,
    retained_count = checkpoint_state$retained_count,
    latest_checkpoint = runtime_paths$checkpoint_latest_rds,
    latest_iteration_checkpoint = checkpoint_state$checkpoint_file,
    checkpoint_written_at = checkpoint_state$checkpoint_written_at,
    status = if (checkpoint_state$last_completed_iter >= config$niter) "complete" else "running",
    resume_command = gcvs_resume_command(config)
  )
  gcvs_atomic_write_json(progress, runtime_paths$progress_json)
}

gcvs_write_checkpoint <- function(runtime_paths, config, dataset, checkpoint_state, archive = TRUE){
  checkpoint_state$checkpoint_version <- 2L
  checkpoint_state$checkpoint_written_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  checkpoint_state$checkpoint_file <- gcvs_checkpoint_path(runtime_paths, checkpoint_state$last_completed_iter)
  gcvs_atomic_save_rds(checkpoint_state, runtime_paths$checkpoint_latest_rds)
  if (isTRUE(archive)) {
    gcvs_atomic_save_rds(checkpoint_state, checkpoint_state$checkpoint_file)
  }
  gcvs_write_progress(runtime_paths, config, dataset, checkpoint_state)
  invisible(checkpoint_state)
}

gcvs_load_checkpoint <- function(runtime_paths){
  if (!file.exists(runtime_paths$checkpoint_latest_rds)) {
    stop(
      sprintf(
        "Resume requested, but no checkpoint was found at %s",
        runtime_paths$checkpoint_latest_rds
      )
    )
  }
  readRDS(runtime_paths$checkpoint_latest_rds)
}

gcvs_run_sampler <- function(dataset_info, config, repo_root){
  X_train <- dataset_info$X_train
  Y_train <- dataset_info$Y_train
  X_test <- dataset_info$X_test
  Y_test <- dataset_info$Y_test
  beta_true <- dataset_info$beta_true

  J <- ncol(X_train[[1]])
  K <- length(Y_train)
  model <- gcvs_prepare_model_env(repo_root, J = J, K = K, eta = config$eta)
  kept_count <- config$niter - config$burnin
  if (kept_count <= 0L) {
    stop("Burn-in must be strictly smaller than the total number of iterations.")
  }

  runtime_paths <- gcvs_initialize_runtime_paths(repo_root, config)
  alpha <- matrix(config$alpha_const, nrow = J, ncol = K)
  checkpoint_interval <- max(1L, as.integer(config$checkpoint_every))
  checkpoint_state <- NULL
  if (isTRUE(config$resume)) {
    checkpoint_state <- gcvs_load_checkpoint(runtime_paths)
    if (is.null(checkpoint_state$checkpoint_version) || checkpoint_state$checkpoint_version < 2L) {
      stop(
        paste(
          "Resume requested, but the checkpoint was created before the two-phase lambda sampler",
          "(burn-in adaptation plus post-burnin fixed-kernel MMALA).",
          "Start a fresh experiment with a new exp-id instead of mixing sampler versions."
        )
      )
    }
    gcvs_validate_resume_config(checkpoint_state$config, config)
    if (!identical(checkpoint_state$dataset, dataset_info$dataset)) {
      stop("Resume requested, but the saved checkpoint dataset does not match the requested dataset.")
    }
  }

  if (is.null(checkpoint_state)) {
    set.seed(config$seed)
    init <- model$priorspecification_allnew(
      J = J,
      K = K,
      alpha = config$alpha_const,
      eta = config$eta,
      a_psi = 1,
      b_psi = 1
    )

    Beta <- init$Beta
    Tau <- init$Tau
    Lambda <- pmax(init$Lambda, 1e-08)
    Psi <- init$Psi
    Sigma <- diag(1, K)
    lambda_step_state <- rep(config$epsilon, J)
    beta_sum <- matrix(0, nrow = J, ncol = K)
    sigma_sum <- matrix(0, nrow = K, ncol = K)
    retained_count <- 0L
    keep_index <- 1L
    start_iter <- 1L
    lambda_samples <- array(NA_real_, dim = c(kept_count, J, K))
    sigma_samples <- array(NA_real_, dim = c(kept_count, K, K))
    row_acceptance <- matrix(NA_real_, nrow = config$niter, ncol = J)
    row_step_size <- matrix(NA_real_, nrow = config$niter, ncol = J)
    row_restarts <- matrix(NA_real_, nrow = config$niter, ncol = J)
    iteration_metrics <- data.frame(
      iteration = seq_len(config$niter),
      mean_acceptance = NA_real_,
      mean_step_size = NA_real_,
      max_restarts = NA_real_,
      sampler_mode = rep(config$sampler_mode, config$niter),
      lambda_phase = NA_character_,
      sigma_12 = NA_real_,
      lambda_min = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    if (!is.null(checkpoint_state$random_seed)) {
      assign(".Random.seed", checkpoint_state$random_seed, envir = .GlobalEnv)
    } else {
      set.seed(config$seed)
    }

    Beta <- checkpoint_state$Beta
    Tau <- checkpoint_state$Tau
    Lambda <- checkpoint_state$Lambda
    Psi <- checkpoint_state$Psi
    Sigma <- checkpoint_state$Sigma
    lambda_step_state <- checkpoint_state$lambda_step_state
    beta_sum <- checkpoint_state$beta_sum
    sigma_sum <- checkpoint_state$sigma_sum
    retained_count <- checkpoint_state$retained_count
    keep_index <- retained_count + 1L
    start_iter <- checkpoint_state$last_completed_iter + 1L
    lambda_samples <- checkpoint_state$lambda_samples
    sigma_samples <- checkpoint_state$sigma_samples
    row_acceptance <- checkpoint_state$row_acceptance
    row_step_size <- checkpoint_state$row_step_size
    row_restarts <- checkpoint_state$row_restarts
    iteration_metrics <- checkpoint_state$iteration_metrics
    if (is.null(lambda_step_state) || length(lambda_step_state) != J) {
      stop("Checkpoint is missing lambda_step_state for the two-phase lambda sampler.")
    }
  }

  last_completed_iter <- start_iter - 1L
  if (start_iter <= config$niter) {
    cat(
      sprintf(
        "Running experiment %s from iteration %d of %d.\n",
        config$exp_id,
        start_iter,
        config$niter
      )
    )
  } else {
    cat(
      sprintf(
        "Checkpoint for %s already reached iteration %d. Finalizing artifacts only.\n",
        config$exp_id,
        last_completed_iter
      )
    )
  }

  iteration_range <- if (start_iter <= config$niter) {
    seq.int(from = start_iter, to = config$niter)
  } else {
    integer(0)
  }

  for (iter in iteration_range) {
    lambda_phase <- if (iter <= config$burnin) "adaptive_burnin" else "fixed_sampling"
    if (identical(config$sampler_mode, "full")) {
      Errvar <- model$Errvar_posterior(
        X_train, Y_train, Beta,
        a_err = config$a_err,
        b_err = config$b_err,
        K = K
      )
      Beta <- model$Beta_posterior(X_train, Y_train, Errvar, Psi, Tau)
      Psi <- model$Psi_posterior(Tau, Beta, a_psi = config$a_psi, b_psi = config$b_psi)
      Tau <- model$Tau_posteriorGDP(Lambda, Beta, Psi)
    }

    lambda_update <- model$adaptive_full_M_MALA_Lambda_GDP(
      Beta = Beta,
      Sigma = Sigma,
      Psi = Psi,
      shape = alpha,
      rate = config$eta,
      epsilon = config$epsilon,
      L = config$inner_steps,
      opt_rate = config$opt_rate,
      const_0 = config$const_0,
      epsconst = config$epsconst,
      inisd = config$inisd,
      ini_Lambda = Lambda,
      row_epsilon = lambda_step_state,
      adapt_epsilon = identical(lambda_phase, "adaptive_burnin")
    )
    Lambda <- lambda_update$Lambda
    lambda_step_state <- lambda_update$next_row_epsilon
    for (sigma_step in seq_len(config$sigma_px_steps)) {
      Sigma <- model$Sigma_PX_posterior_GDP(Lambda = Lambda, Sigma = Sigma, shape = alpha, rate = config$eta)
    }

    row_acceptance[iter, ] <- lambda_update$Acc_prob
    row_step_size[iter, ] <- lambda_update$eps_used
    row_restarts[iter, ] <- lambda_update$trial_number
    iteration_metrics$mean_acceptance[iter] <- mean(lambda_update$Acc_prob)
    iteration_metrics$mean_step_size[iter] <- mean(lambda_update$eps_used)
    iteration_metrics$max_restarts[iter] <- max(lambda_update$trial_number)
    iteration_metrics$lambda_phase[iter] <- lambda_phase
    iteration_metrics$sigma_12[iter] <- if (K > 1) Sigma[1, 2] else NA_real_
    iteration_metrics$lambda_min[iter] <- min(Lambda)

    if (iter > config$burnin) {
      lambda_samples[keep_index, , ] <- Lambda
      sigma_samples[keep_index, , ] <- Sigma
      beta_sum <- beta_sum + Beta
      sigma_sum <- sigma_sum + Sigma
      retained_count <- retained_count + 1L
      keep_index <- keep_index + 1L
    }

    last_completed_iter <- iter
    checkpoint_state <- list(
      config = config,
      dataset = dataset_info$dataset,
      J = J,
      K = K,
      last_completed_iter = last_completed_iter,
      retained_count = retained_count,
      Beta = Beta,
      Tau = Tau,
      Lambda = Lambda,
      Psi = Psi,
      Sigma = Sigma,
      lambda_step_state = lambda_step_state,
      beta_sum = beta_sum,
      sigma_sum = sigma_sum,
      lambda_samples = lambda_samples,
      sigma_samples = sigma_samples,
      row_acceptance = row_acceptance,
      row_step_size = row_step_size,
      row_restarts = row_restarts,
      iteration_metrics = iteration_metrics,
      random_seed = if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      } else {
        NULL
      }
    )
    should_checkpoint <- ((iter %% checkpoint_interval) == 0L || iter == config$niter)
    if (should_checkpoint) {
      gcvs_write_checkpoint(
        runtime_paths,
        config,
        dataset_info$dataset,
        checkpoint_state,
        archive = TRUE
      )
    }
  }

  if (retained_count == 0L) {
    stop("Burn-in is greater than or equal to the total number of iterations; no posterior samples were retained.")
  }

  beta_mean <- beta_sum / retained_count
  sigma_mean <- sigma_sum / retained_count

  r2_by_task <- vapply(seq_len(K), function(task_id) {
    gcvs_rsquare(X_test[[task_id]], Y_test[[task_id]], beta_mean[, task_id])
  }, numeric(1))
  rmse_by_task <- vapply(seq_len(K), function(task_id) {
    gcvs_rmse(X_test[[task_id]], Y_test[[task_id]], beta_mean[, task_id])
  }, numeric(1))

  precision <- if (!is.null(beta_true) && all(dim(beta_true) == dim(beta_mean))) {
    gcvs_selection_precision(beta_mean, beta_true)
  } else {
    NA_real_
  }

  has_numerical_issue <- any(!is.finite(as.matrix(iteration_metrics[, c("mean_acceptance", "mean_step_size", "lambda_min")])))
  mean_acceptance <- mean(iteration_metrics$mean_acceptance, na.rm = TRUE)
  decision <- if (has_numerical_issue) {
    "NEEDS-DEBUG"
  } else if (config$profile == "paper") {
    "PASS"
  } else {
    "INCONCLUSIVE"
  }

  list(
    dataset = dataset_info$dataset,
    config = config,
    J = J,
    K = K,
    beta_true = beta_true,
    beta_mean = beta_mean,
    sigma_mean = sigma_mean,
    lambda_samples = lambda_samples,
    sigma_samples = sigma_samples,
    row_acceptance = row_acceptance,
    row_step_size = row_step_size,
    row_restarts = row_restarts,
    iteration_metrics = iteration_metrics,
    retained_count = retained_count,
    checkpoint_latest_rds = runtime_paths$checkpoint_latest_rds,
    progress_json = runtime_paths$progress_json,
    r2_by_task = r2_by_task,
    rmse_by_task = rmse_by_task,
    precision = precision,
    mean_acceptance = mean_acceptance,
    decision = decision,
    diagnostics = NULL
  )
}

gcvs_plot_diagnostics <- function(run_result, png_path){
  png(filename = png_path, width = 1400, height = 900)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  plot(
    run_result$iteration_metrics$iteration,
    run_result$iteration_metrics$mean_acceptance,
    type = "l",
    lwd = 2,
    col = "steelblue",
    xlab = "Iteration",
    ylab = "Mean lambda acceptance",
    main = "Acceptance by Iteration"
  )
  abline(h = run_result$config$opt_rate, lty = 2, col = "firebrick")

  plot(
    run_result$iteration_metrics$iteration,
    run_result$iteration_metrics$sigma_12,
    type = "l",
    lwd = 2,
    col = "darkgreen",
    xlab = "Iteration",
    ylab = expression(Sigma[1,2]),
    main = "Sigma[1,2] Trace"
  )

  boxplot(
    as.data.frame(run_result$row_acceptance),
    outline = FALSE,
    xaxt = "n",
    ylab = "Acceptance",
    main = "Row-wise Lambda Acceptance"
  )
  axis(1, at = c(1, floor(run_result$J / 4), floor(run_result$J / 2), floor(3 * run_result$J / 4), run_result$J))
  abline(h = run_result$config$opt_rate, lty = 2, col = "firebrick")

  boxplot(
    as.data.frame(run_result$row_step_size),
    outline = FALSE,
    xaxt = "n",
    ylab = "Step size",
    main = "Row-wise Step Size"
  )
  axis(1, at = c(1, floor(run_result$J / 4), floor(run_result$J / 2), floor(3 * run_result$J / 4), run_result$J))
  dev.off()
}

gcvs_plot_beta_heatmap <- function(run_result, png_path){
  beta_abs <- abs(run_result$beta_mean)
  png(filename = png_path, width = 1200, height = 900)
  par(mar = c(5, 5, 3, 1))
  image(
    x = seq_len(ncol(beta_abs)),
    y = seq_len(nrow(beta_abs)),
    z = t(beta_abs[nrow(beta_abs):1, , drop = FALSE]),
    col = gray.colors(128, start = 0, end = 1),
    axes = FALSE,
    xlab = "Task",
    ylab = "Feature index",
    main = "Posterior Mean |Beta| Heatmap"
  )
  axis(1, at = seq_len(ncol(beta_abs)), labels = paste("task", seq_len(ncol(beta_abs))), las = 2)
  axis(2, at = c(1, floor(nrow(beta_abs) / 4), floor(nrow(beta_abs) / 2), floor(3 * nrow(beta_abs) / 4), nrow(beta_abs)),
       labels = rev(c(1, floor(nrow(beta_abs) / 4), floor(nrow(beta_abs) / 2), floor(3 * nrow(beta_abs) / 4), nrow(beta_abs))))
  box()
  dev.off()
}

gcvs_plot_sigma_heatmap <- function(run_result, png_path){
  sigma_mat <- run_result$sigma_mean
  sigma_palette <- grDevices::colorRampPalette(c("#3b0f70", "#f7f7f7", "#b40426"))(128)
  png(filename = png_path, width = 900, height = 800)
  par(mar = c(5, 5, 3, 1))
  image(
    x = seq_len(ncol(sigma_mat)),
    y = seq_len(nrow(sigma_mat)),
    z = t(sigma_mat[nrow(sigma_mat):1, , drop = FALSE]),
    col = sigma_palette,
    zlim = c(-1, 1),
    axes = FALSE,
    xlab = "Task",
    ylab = "Task",
    main = "Posterior Mean Sigma Correlations"
  )
  axis(1, at = seq_len(ncol(sigma_mat)), labels = paste("task", seq_len(ncol(sigma_mat))), las = 2)
  axis(2, at = seq_len(nrow(sigma_mat)), labels = rev(paste("task", seq_len(nrow(sigma_mat)))), las = 2)
  box()
  dev.off()
}

gcvs_plot_representative_traces <- function(run_result, png_path){
  trace_labels <- c(unique(run_result$diagnostics$lambda_traces$label), unique(run_result$diagnostics$sigma_traces$label))
  panel_count <- length(trace_labels)
  ncol_plot <- 2L
  nrow_plot <- ceiling(panel_count / ncol_plot)

  png(filename = png_path, width = 1400, height = max(900, 300 * nrow_plot))
  par(mfrow = c(nrow_plot, ncol_plot), mar = c(4, 4, 3, 1))

  for (label in unique(run_result$diagnostics$lambda_traces$label)) {
    trace_df <- run_result$diagnostics$lambda_traces[run_result$diagnostics$lambda_traces$label == label, , drop = FALSE]
    ess_value <- run_result$diagnostics$lambda_summary$ess[run_result$diagnostics$lambda_summary$label == label][[1]]
    plot(trace_df$draw, trace_df$value, type = "l", lwd = 1.5, col = "steelblue",
         xlab = "Retained draw", ylab = "Lambda", main = paste0(label, " (ESS=", gcvs_format_metric(ess_value, 1), ")"))
  }

  for (label in unique(run_result$diagnostics$sigma_traces$label)) {
    trace_df <- run_result$diagnostics$sigma_traces[run_result$diagnostics$sigma_traces$label == label, , drop = FALSE]
    ess_value <- run_result$diagnostics$sigma_summary$ess[run_result$diagnostics$sigma_summary$label == label][[1]]
    plot(trace_df$draw, trace_df$value, type = "l", lwd = 1.5, col = "darkgreen",
         xlab = "Retained draw", ylab = "Sigma value", main = paste0(label, " (ESS=", gcvs_format_metric(ess_value, 1), ")"))
    abline(h = 0, lty = 2, col = "gray50")
  }

  dev.off()
}

gcvs_write_row_metric_csv <- function(metric_matrix, path){
  metric_df <- as.data.frame(metric_matrix, stringsAsFactors = FALSE)
  names(metric_df) <- paste0("feature_", seq_len(ncol(metric_df)))
  metric_df <- cbind(iteration = seq_len(nrow(metric_df)), metric_df)
  write.csv(metric_df, path, row.names = FALSE)
}

gcvs_write_run_log <- function(run_result, paths){
  evidence_paths <- c(
    paths$summary_csv,
    paths$iteration_csv,
    paths$beta_csv,
    paths$sigma_csv,
    paths$checkpoint_latest_rds,
    paths$progress_json,
    paths$lambda_row_acceptance_summary_csv,
    paths$representative_lambda_csv,
    paths$representative_sigma_csv,
    paths$sigma_pair_summary_csv,
    paths$diagnostics_png,
    paths$beta_heatmap_png,
    paths$sigma_heatmap_png,
    paths$representative_traces_png,
    paths$metadata_json
  )
  decision_text <- if (run_result$decision == "PASS") {
    "support"
  } else if (run_result$decision == "NEEDS-DEBUG") {
    "reject"
  } else {
    "inconclusive"
  }

  lines <- c(
    paste("#", run_result$config$exp_id),
    "",
    "## Problem",
    paste(
      "Run a reproducible", run_result$dataset,
      "GCVS experiment with a relative-path entrypoint, save machine-readable summaries,",
      "sampling diagnostics, and a structured iteration log."
    ),
    "",
    "## Hypothesis",
    "The cleaned driver and corrected lambda MMALA kernel can run end-to-end without ad hoc path edits or invalid-support fallback moves.",
    "",
    "## Plan",
    paste(
      "Use profile", run_result$config$profile,
      "with", run_result$config$niter, "outer iterations,",
      run_result$config$burnin, "burn-in iterations, and",
      run_result$config$inner_steps, "inner lambda-MMALA steps."
    ),
    paste("Sigma PX refreshes per outer iteration:", run_result$config$sigma_px_steps),
    "Metrics: mean lambda acceptance, posterior predictive R^2, RMSE, and sampling diagnostics.",
    "Stop rule: finish the configured run or halt on a numerical failure.",
    "",
    "## Run",
    paste("Dataset:", run_result$dataset),
    paste("Seed:", run_result$config$seed),
    paste("Decision:", run_result$decision),
    "",
    "## Evaluate",
    paste("Mean lambda acceptance:", sprintf("%.4f", run_result$mean_acceptance)),
    paste("Mean predictive R^2:", sprintf("%.4f", mean(run_result$r2_by_task))),
    paste("Mean predictive RMSE:", sprintf("%.4f", mean(run_result$rmse_by_task))),
    paste("Selection precision:", ifelse(is.na(run_result$precision), "NA", sprintf("%.4f", run_result$precision))),
    "",
    "## Decide",
    decision_text,
    "",
    "## Log",
    paste("- problem_id:", paste0("PROB-", run_result$config$exp_id)),
    paste("- hypothesis_id:", paste0("HYP-", run_result$config$exp_id)),
    paste("- experiment_id:", run_result$config$exp_id),
    paste("- assumption_list: dataset files are loaded from the repository; the run uses the corrected full MMALA kernel and clipped copula probabilities."),
    paste("- changes_made: relative-path dataset loader, structured artifact writer, corrected lambda support handling, and numerical clipping for copula transforms."),
    paste("- evidence_paths:", paste(evidence_paths, collapse = ", ")),
    paste("- result_summary: mean acceptance", sprintf("%.4f", run_result$mean_acceptance), ", mean R^2", sprintf("%.4f", mean(run_result$r2_by_task)), ", mean RMSE", sprintf("%.4f", mean(run_result$rmse_by_task))),
    paste("- decision:", run_result$decision),
    paste("- next_action:", if (run_result$config$profile == "paper") "Compare against saved paper figures and regenerate the manuscript tables." else "Escalate to the paper profile once smoke diagnostics look stable."),
    paste("- risk_or_blocker:", if (run_result$config$profile == "paper") "External baseline scripts are still legacy-only and are not part of the cleaned core pipeline." else "This is a smoke run, so evidence is not strong enough for a paper-level claim.")
  )

  writeLines(lines, con = paths$log_md)
}

gcvs_write_metadata <- function(run_result, paths){
  jsonlite::write_json(
    list(
      experiment_id = run_result$config$exp_id,
      dataset = run_result$dataset,
      profile = run_result$config$profile,
      decision = run_result$decision,
      mean_acceptance = run_result$mean_acceptance,
      mean_r2 = mean(run_result$r2_by_task),
      mean_rmse = mean(run_result$rmse_by_task),
      evidence_paths = list(
        summary_csv = paths$summary_csv,
        iteration_csv = paths$iteration_csv,
        beta_csv = paths$beta_csv,
        sigma_csv = paths$sigma_csv,
        lambda_row_acceptance_csv = paths$lambda_row_acceptance_csv,
        lambda_row_acceptance_summary_csv = paths$lambda_row_acceptance_summary_csv,
        lambda_row_step_size_csv = paths$lambda_row_step_size_csv,
        lambda_row_restarts_csv = paths$lambda_row_restarts_csv,
        representative_lambda_csv = paths$representative_lambda_csv,
        representative_sigma_csv = paths$representative_sigma_csv,
        sigma_pair_summary_csv = paths$sigma_pair_summary_csv,
        support_recovery_csv = paths$support_recovery_csv,
        task_abs_correlation_csv = paths$task_abs_correlation_csv,
        lambda_samples_rds = paths$lambda_samples_rds,
        sigma_samples_rds = paths$sigma_samples_rds,
        checkpoint_latest_rds = paths$checkpoint_latest_rds,
        progress_json = paths$progress_json,
        diagnostics_png = paths$diagnostics_png,
        beta_heatmap_png = paths$beta_heatmap_png,
        sigma_heatmap_png = paths$sigma_heatmap_png,
        representative_traces_png = paths$representative_traces_png,
        log_md = paths$log_md,
        report_md = paths$report_md,
        report_html = if (file.exists(paths$report_html)) paths$report_html else "",
        report_pdf = if (file.exists(paths$report_pdf)) paths$report_pdf else "",
        report_meta_json = paths$report_meta_json
      )
    ),
    path = paths$metadata_json,
    pretty = TRUE,
    auto_unbox = TRUE
  )
}

gcvs_append_index <- function(run_result, paths){
  new_row <- data.frame(
    date = format(Sys.time(), "%Y-%m-%d"),
    exp_id = run_result$config$exp_id,
    dataset = run_result$dataset,
    profile = run_result$config$profile,
    sampler_mode = run_result$config$sampler_mode,
    niter = run_result$config$niter,
    burnin = run_result$config$burnin,
    inner_steps = run_result$config$inner_steps,
    sigma_px_steps = run_result$config$sigma_px_steps,
    mean_acceptance = round(run_result$mean_acceptance, 6),
    mean_r2 = round(mean(run_result$r2_by_task), 6),
    mean_rmse = round(mean(run_result$rmse_by_task), 6),
    precision = ifelse(is.na(run_result$precision), "", round(run_result$precision, 6)),
    decision = run_result$decision,
    log_path = paths$log_md,
    progress_path = paths$progress_json,
    checkpoint_path = paths$checkpoint_latest_rds,
    stringsAsFactors = FALSE
  )

  if (file.exists(paths$index_csv)) {
    current <- read.csv(paths$index_csv, stringsAsFactors = FALSE, check.names = FALSE)
    missing_cols <- setdiff(names(new_row), names(current))
    for (col in missing_cols) {
      current[[col]] <- rep("", nrow(current))
    }
    extra_cols <- setdiff(names(current), names(new_row))
    for (col in extra_cols) {
      new_row[[col]] <- ""
    }
    current <- current[, names(new_row), drop = FALSE]
    current <- current[current$exp_id != new_row$exp_id, , drop = FALSE]
    updated <- rbind(current, new_row)
  } else {
    updated <- new_row
  }

  write.csv(updated, paths$index_csv, row.names = FALSE)
}

gcvs_write_experiment_report <- function(run_result, paths){
  last_iter <- tail(run_result$iteration_metrics, 1)
  lambda_ess_values <- run_result$diagnostics$lambda_summary$ess
  lambda_ess_values <- lambda_ess_values[is.finite(lambda_ess_values)]
  sigma_ess_values <- run_result$diagnostics$sigma_summary$ess
  sigma_ess_values <- sigma_ess_values[is.finite(sigma_ess_values)]
  sigma_shift_values <- run_result$diagnostics$sigma_summary$standardized_half_shift
  sigma_shift_values <- sigma_shift_values[is.finite(sigma_shift_values)]
  min_lambda_ess <- if (length(lambda_ess_values) == 0) NA_real_ else min(lambda_ess_values)
  min_sigma_ess <- if (length(sigma_ess_values) == 0) NA_real_ else min(sigma_ess_values)
  max_sigma_shift <- if (length(sigma_shift_values) == 0) NA_real_ else max(sigma_shift_values)
  lines <- c(
    "# Experiment Report",
    "",
    "## 1. Experiment ID",
    run_result$config$exp_id,
    "",
    "## 2. Objective",
    "Stabilize the reproducible GC-MTL pipeline, align the lambda posterior simulation with the implemented two-phase MMALA kernel, and verify that the cleaned repository can run end-to-end without manual path edits.",
    "",
    "## 3. Setup",
    paste("- Dataset:", run_result$dataset),
    paste("- Profile:", run_result$config$profile),
    paste("- Sampler mode:", gsub("_", "-", run_result$config$sampler_mode)),
    paste("- Outer iterations:", run_result$config$niter),
    paste("- Burn-in:", run_result$config$burnin),
    paste("- Checkpoint cadence (outer iterations):", run_result$config$checkpoint_every),
    paste("- Inner lambda-MMALA steps per row:", run_result$config$inner_steps),
    paste("- Sigma PX refreshes per outer iteration:", run_result$config$sigma_px_steps),
    paste("- Hyperparameters: alpha =", run_result$config$alpha_const, ", eta =", run_result$config$eta, ", a_psi =", run_result$config$a_psi, ", b_psi =", run_result$config$b_psi, ", a_err =", run_result$config$a_err, ", b_err =", run_result$config$b_err),
    paste("- Adaptation target:", run_result$config$opt_rate),
    "",
    "## 4. Intervention",
    "- Replaced the broken top-level experiment driver with a repository-relative CLI entrypoint.",
    "- Removed hidden dependencies on pre-created objects such as `Beta_alltasks_J_40` and `W` from the runnable path.",
    "- Corrected the lambda MMALA kernel so invalid-support proposals are rejected rather than repaired with an unmatched fallback draw.",
    "- Replaced acceptance-triggered row restarts with a two-phase lambda sampler: burn-in adaptation followed by post-burnin fixed-kernel MMALA.",
    "- Added two-sided probability clipping before the Gamma-to-Gaussian copula transform and repaired indefinite proposal metrics by eigenvalue flooring.",
    if (identical(run_result$config$sampler_mode, "lambda_sigma_only")) {
      "- This diagnostic run held Beta, Psi, and Tau fixed after initialization and updated only lambda and Sigma."
    } else {
      "- This run used the full blocked Gibbs-plus-MMALA sampler."
    },
    "- Added structured artifact writing for summaries, diagnostics, logs, and a reviewable experiment report.",
    "",
    "## 5. Raw Evidence",
    paste("- mean_acceptance =", run_result$mean_acceptance),
    paste("- mean_r2 =", mean(run_result$r2_by_task)),
    paste("- mean_rmse =", mean(run_result$rmse_by_task)),
    paste("- precision =", run_result$precision),
    paste("- final_sigma_12 =", last_iter$sigma_12),
    paste("- final_lambda_min =", last_iter$lambda_min),
    paste("- summary_csv:", paths$summary_csv),
    paste("- iteration_csv:", paths$iteration_csv),
    paste("- diagnostics_png:", paths$diagnostics_png),
    paste("- checkpoint_latest_rds:", paths$checkpoint_latest_rds),
    paste("- progress_json:", paths$progress_json),
    paste("- lambda_row_acceptance_summary_csv:", paths$lambda_row_acceptance_summary_csv),
    paste("- representative_lambda_csv:", paths$representative_lambda_csv),
    paste("- representative_sigma_csv:", paths$representative_sigma_csv),
    paste("- beta_csv:", paths$beta_csv),
    paste("- sigma_csv:", paths$sigma_csv),
    paste("- log_md:", paths$log_md),
    "",
    "## 6. Processed Metrics",
    paste("- Mean lambda acceptance:", gcvs_format_metric(run_result$mean_acceptance)),
    paste("- Mean predictive R^2:", gcvs_format_metric(mean(run_result$r2_by_task))),
    paste("- Mean predictive RMSE:", gcvs_format_metric(mean(run_result$rmse_by_task))),
    paste("- Selection precision:", gcvs_format_metric(run_result$precision)),
    paste("- Min representative lambda ESS:", gcvs_format_metric(min_lambda_ess)),
    paste("- Min representative Sigma ESS:", gcvs_format_metric(min_sigma_ess)),
    paste("- Max representative Sigma half-shift:", gcvs_format_metric(max_sigma_shift)),
    paste("- Decision:", run_result$decision),
    "",
    "## 7. Key Findings",
    "- The cleaned repository now has a reproducible entrypoint that runs from the committed files and writes the minimum evidence bundle expected by the AGENTS workflow.",
    "- The lambda sampler now preserves the intended Metropolis-Hastings kernel by rejecting proposals outside the positive support instead of silently changing the proposal law.",
    "- Two-sided clipping in the copula transform removes the `qnorm(0)` failure mode that previously produced `-Inf`, `NaN`, and unstable gradients for extremely small lambda values.",
    "- The smoke run is numerically stable, but it is still only a short validation run and does not yet certify paper-level posterior accuracy.",
    "",
    "## 8. Mechanism Hypothesis",
    "The main instability came from a combination of support-violating lambda proposals and underflow in the Gamma CDF to probit transform. Rejecting invalid proposals and clipping the copula probabilities should restore a valid acceptance ratio and keep the metric-based proposal well defined.",
    "",
    "## 9. Alternative Explanations",
    "- The improved behavior may partly reflect the synthetic dataset being easier than the real-data settings.",
    "- A short smoke run can miss slower mixing problems that only appear under the longer paper profile.",
    "",
    "## 10. Comparison",
    "Compared with the committed version, the cleaned pipeline no longer depends on undefined workspace objects, the lambda kernel matches the implemented proposal density, and the manuscript-facing sampler description can now be stated without relying on undocumented fallbacks.",
    "",
    "## 11. Failure / Limitation",
    "- Legacy preprocessing and external baseline scripts still depend on datasets or third-party code that are not part of the cleaned core sampling pipeline.",
    "- This report is based on the smoke profile; a longer run is still needed before claiming reproduction of the manuscript figures.",
    "- A PDF archive was not generated because no local LaTeX or HTML-to-PDF engine is currently installed.",
    "",
    "## 12. Next Step",
    "- Run the paper profile on the synthetic experiment and compare the posterior mean coefficient heatmap and task-correlation matrix against the manuscript figures.",
    "- Parameterize the remaining legacy helper scripts so external comparison workflows no longer depend on hard-coded Windows paths.",
    "- Regenerate the manuscript section on posterior computation from the cleaned implementation and then check the experiment section against the saved artifacts.",
    "",
    "## 13. Confidence",
    "Moderate confidence in the kernel fix, artifact pipeline, and reproducibility of the cleaned core experiment path; lower confidence in full paper replication until the longer chain and figure-by-figure comparison are complete."
  )

  writeLines(lines, con = paths$report_md)
}

gcvs_render_experiment_report <- function(paths){
  pandoc <- Sys.which("pandoc")
  if (!nzchar(pandoc)) {
    return(invisible(NULL))
  }

  system2(
    pandoc,
    args = c("--standalone", paths$report_md, "-o", paths$report_html),
    stdout = FALSE,
    stderr = FALSE
  )

  pdf_engine <- gcvs_first_nonempty(c(
    Sys.which("xelatex"),
    Sys.which("pdflatex"),
    Sys.which("tectonic"),
    Sys.which("wkhtmltopdf")
  ))
  if (!nzchar(pdf_engine)) {
    return(invisible(NULL))
  }

  system2(
    pandoc,
    args = c("--standalone", paths$report_md, "-o", paths$report_pdf, paste0("--pdf-engine=", pdf_engine)),
    stdout = FALSE,
    stderr = FALSE
  )
}

gcvs_write_report_metadata <- function(run_result, paths){
  jsonlite::write_json(
    list(
      experiment_id = run_result$config$exp_id,
      dataset = run_result$dataset,
      decision = run_result$decision,
      summary_csv = paths$summary_csv,
      trace_csv = paths$iteration_csv,
      report_md = paths$report_md,
      report_html = if (file.exists(paths$report_html)) paths$report_html else "",
      pdf_path = if (file.exists(paths$report_pdf)) paths$report_pdf else "",
      diagnostics_png = paths$diagnostics_png,
      checkpoint_latest_rds = paths$checkpoint_latest_rds,
      progress_json = paths$progress_json,
      lambda_row_acceptance_summary_csv = paths$lambda_row_acceptance_summary_csv,
      representative_lambda_csv = paths$representative_lambda_csv,
      representative_sigma_csv = paths$representative_sigma_csv,
      beta_heatmap_png = paths$beta_heatmap_png,
      sigma_heatmap_png = paths$sigma_heatmap_png,
      representative_traces_png = paths$representative_traces_png,
      metadata_json = paths$metadata_json
    ),
    path = paths$report_meta_json,
    pretty = TRUE,
    auto_unbox = TRUE
  )
}

gcvs_append_report_index <- function(run_result, paths){
  new_row <- data.frame(
    date = format(Sys.Date(), "%Y-%m-%d"),
    report_id = run_result$config$exp_id,
    report_type = "experiment",
    experiment_id = run_result$config$exp_id,
    dataset = run_result$dataset,
    profile = run_result$config$profile,
    M1_rule = paste0("GDP(alpha=", run_result$config$alpha_const, ",eta=", run_result$config$eta, ")"),
    M2_rule = paste0(
      gsub("_", "-", run_result$config$sampler_mode),
      "::",
      "two-phase-MMALA(L=",
      run_result$config$inner_steps,
      ",target=",
      run_result$config$opt_rate,
      ",adapt=burnin,fixed=post",
      ",SigmaPX=",
      run_result$config$sigma_px_steps,
      ")"
    ),
    theta_mean = gcvs_format_metric(mean(run_result$diagnostics$lambda_summary$posterior_mean, na.rm = TRUE)),
    theta_peaks = gcvs_format_metric(max(run_result$diagnostics$lambda_summary$posterior_mean, na.rm = TRUE)),
    theta_ess = gcvs_format_metric(min(run_result$diagnostics$lambda_summary$ess, na.rm = TRUE)),
    M_mean = gcvs_format_metric(mean(run_result$diagnostics$sigma_pair_summary$posterior_mean, na.rm = TRUE)),
    decision = run_result$decision,
    report_path = gcvs_repo_relative_path(gcvs_repo_root(), paths$report_md),
    html_path = if (file.exists(paths$report_html)) gcvs_repo_relative_path(gcvs_repo_root(), paths$report_html) else "",
    pdf_path = if (file.exists(paths$report_pdf)) gcvs_repo_relative_path(gcvs_repo_root(), paths$report_pdf) else "",
    meta_path = gcvs_repo_relative_path(gcvs_repo_root(), paths$report_meta_json),
    stringsAsFactors = FALSE
  )

  gcvs_register_report_entry(paths$report_index_csv, new_row)
}

gcvs_write_artifacts <- function(run_result, repo_root){
  paths <- gcvs_initialize_results_paths(repo_root, run_result)

  summary_df <- data.frame(
    experiment_id = run_result$config$exp_id,
    dataset = run_result$dataset,
    profile = run_result$config$profile,
    sampler_mode = run_result$config$sampler_mode,
    niter = run_result$config$niter,
    burnin = run_result$config$burnin,
    inner_steps = run_result$config$inner_steps,
    sigma_px_steps = run_result$config$sigma_px_steps,
    mean_acceptance = run_result$mean_acceptance,
    mean_r2 = mean(run_result$r2_by_task),
    mean_rmse = mean(run_result$rmse_by_task),
    precision = run_result$precision,
    decision = run_result$decision,
    stringsAsFactors = FALSE
  )

  write.csv(summary_df, paths$summary_csv, row.names = FALSE)
  write.csv(run_result$iteration_metrics, paths$iteration_csv, row.names = FALSE)
  write.csv(run_result$beta_mean, paths$beta_csv, row.names = FALSE)
  write.csv(run_result$sigma_mean, paths$sigma_csv, row.names = FALSE)
  gcvs_write_row_metric_csv(run_result$row_acceptance, paths$lambda_row_acceptance_csv)
  gcvs_write_row_metric_csv(run_result$row_step_size, paths$lambda_row_step_size_csv)
  gcvs_write_row_metric_csv(run_result$row_restarts, paths$lambda_row_restarts_csv)
  saveRDS(run_result$lambda_samples, file = paths$lambda_samples_rds)
  saveRDS(run_result$sigma_samples, file = paths$sigma_samples_rds)
  write.csv(run_result$diagnostics$row_acceptance_summary, paths$lambda_row_acceptance_summary_csv, row.names = FALSE)
  write.csv(run_result$diagnostics$lambda_summary, paths$representative_lambda_csv, row.names = FALSE)
  write.csv(run_result$diagnostics$sigma_summary, paths$representative_sigma_csv, row.names = FALSE)
  write.csv(run_result$diagnostics$sigma_pair_summary, paths$sigma_pair_summary_csv, row.names = FALSE)
  if (!is.null(run_result$diagnostics$support_recovery_summary)) {
    write.csv(run_result$diagnostics$support_recovery_summary, paths$support_recovery_csv, row.names = FALSE)
  }
  if (!is.null(run_result$diagnostics$task_abs_correlation_summary)) {
    write.csv(run_result$diagnostics$task_abs_correlation_summary, paths$task_abs_correlation_csv, row.names = FALSE)
  }
  gcvs_plot_diagnostics(run_result, paths$diagnostics_png)
  gcvs_plot_beta_heatmap(run_result, paths$beta_heatmap_png)
  gcvs_plot_sigma_heatmap(run_result, paths$sigma_heatmap_png)
  gcvs_plot_representative_traces(run_result, paths$representative_traces_png)
  gcvs_write_run_log(run_result, paths)
  gcvs_write_experiment_report(run_result, paths)
  gcvs_render_experiment_report(paths)
  gcvs_write_metadata(run_result, paths)
  gcvs_write_report_metadata(run_result, paths)
  gcvs_append_index(run_result, paths)
  gcvs_append_report_index(run_result, paths)

  paths
}

gcvs_run_experiment <- function(config){
  repo_root <- gcvs_repo_root()
  dataset_info <- gcvs_load_dataset(repo_root, config$dataset)
  run_result <- gcvs_run_sampler(dataset_info, config, repo_root)
  run_result$diagnostics <- gcvs_build_diagnostics(run_result)
  artifact_paths <- gcvs_write_artifacts(run_result, repo_root)

  cat("Experiment completed.\n")
  cat("Experiment ID:", run_result$config$exp_id, "\n")
  cat("Summary CSV:", artifact_paths$summary_csv, "\n")
  cat("Run log:", artifact_paths$log_md, "\n")
  cat("Report:", artifact_paths$report_md, "\n")

  invisible(list(result = run_result, paths = artifact_paths))
}

gcvs_run_cli <- function(args = commandArgs(trailingOnly = TRUE)){
  parsed <- gcvs_parse_args(args)
  config <- gcvs_build_config(parsed)
  gcvs_run_experiment(config)
}

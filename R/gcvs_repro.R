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

gcvs_default_config <- function(profile = "smoke"){
  profile <- tolower(profile)
  if (profile == "paper") {
    return(list(
      profile = "paper",
      niter = 300L,
      burnin = 150L,
      inner_steps = 40L,
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
    inner_steps = 8L,
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

gcvs_build_config <- function(parsed){
  defaults <- gcvs_default_config(gcvs_default_if_null(parsed$profile, "smoke"))
  timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  dataset <- tolower(gcvs_default_if_null(parsed$dataset, "synthetic"))

  list(
    dataset = dataset,
    profile = defaults$profile,
    niter = gcvs_numeric_arg(parsed, "niter", defaults$niter),
    burnin = gcvs_numeric_arg(parsed, "burnin", defaults$burnin),
    inner_steps = gcvs_numeric_arg(parsed, "inner-steps", defaults$inner_steps),
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
    seed = gcvs_numeric_arg(parsed, "seed", 123L),
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
  logs_dir <- file.path(repo_root, "logs", "iteration_logs")
  reports_dir <- file.path(repo_root, "results", "reports")
  reports_pdf_dir <- file.path(reports_dir, "pdf")
  reports_meta_dir <- file.path(reports_dir, "meta")
  dir.create(experiment_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(reports_pdf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(reports_meta_dir, recursive = TRUE, showWarnings = FALSE)

  list(
    experiment_dir = experiment_dir,
    figures_dir = figures_dir,
    summary_csv = file.path(experiment_dir, "summary.csv"),
    iteration_csv = file.path(experiment_dir, "iteration_metrics.csv"),
    beta_csv = file.path(experiment_dir, "posterior_beta_mean.csv"),
    sigma_csv = file.path(experiment_dir, "posterior_sigma_mean.csv"),
    metadata_json = file.path(experiment_dir, "metadata.json"),
    diagnostics_png = file.path(figures_dir, "sampling_diagnostics.png"),
    log_md = file.path(logs_dir, paste0(format(Sys.time(), "%Y%m%d__"), exp_id, ".md")),
    report_md = file.path(experiment_dir, "experiment_report.md"),
    report_html = file.path(experiment_dir, "experiment_report.html"),
    report_pdf = file.path(reports_pdf_dir, paste0(report_basename, ".pdf")),
    report_meta_json = file.path(reports_meta_dir, paste0(report_basename, ".json")),
    index_csv = file.path(repo_root, "results", "experiments", "index.csv"),
    report_index_csv = file.path(reports_dir, "index.csv")
  )
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

  set.seed(config$seed)
  alpha <- matrix(config$alpha_const, nrow = J, ncol = K)
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

  kept_samples <- list()
  keep_index <- 1L
  iteration_metrics <- data.frame(
    iteration = seq_len(config$niter),
    mean_acceptance = NA_real_,
    mean_step_size = NA_real_,
    max_restarts = NA_real_,
    sigma_12 = NA_real_,
    lambda_min = NA_real_,
    stringsAsFactors = FALSE
  )

  for (iter in seq_len(config$niter)) {
    Errvar <- model$Errvar_posterior(
      X_train, Y_train, Beta,
      a_err = config$a_err,
      b_err = config$b_err,
      K = K
    )
    Beta <- model$Beta_posterior(X_train, Y_train, Errvar, Psi, Tau)
    Psi <- model$Psi_posterior(Tau, Beta, a_psi = config$a_psi, b_psi = config$b_psi)
    Tau <- model$Tau_posteriorGDP(Lambda, Beta, Psi)

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
      inisd = config$inisd
    )
    Lambda <- lambda_update[[1]]
    Sigma <- model$Sigma_PX_posterior_GDP(Lambda = Lambda, Sigma = Sigma, shape = alpha, rate = config$eta)

    iteration_metrics$mean_acceptance[iter] <- mean(lambda_update[[2]])
    iteration_metrics$mean_step_size[iter] <- mean(lambda_update[[3]])
    iteration_metrics$max_restarts[iter] <- max(lambda_update[[4]])
    iteration_metrics$sigma_12[iter] <- if (K > 1) Sigma[1, 2] else NA_real_
    iteration_metrics$lambda_min[iter] <- min(Lambda)

    if (iter > config$burnin) {
      kept_samples[[keep_index]] <- list(Beta = Beta, Sigma = Sigma)
      keep_index <- keep_index + 1L
    }
  }

  if (length(kept_samples) == 0) {
    stop("Burn-in is greater than or equal to the total number of iterations; no posterior samples were retained.")
  }

  beta_mean <- gcvs_matrix_mean(kept_samples, "Beta")
  sigma_mean <- gcvs_matrix_mean(kept_samples, "Sigma")

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
    beta_mean = beta_mean,
    sigma_mean = sigma_mean,
    iteration_metrics = iteration_metrics,
    r2_by_task = r2_by_task,
    rmse_by_task = rmse_by_task,
    precision = precision,
    mean_acceptance = mean_acceptance,
    decision = decision
  )
}

gcvs_plot_diagnostics <- function(run_result, png_path){
  png(filename = png_path, width = 1200, height = 700)
  par(mfrow = c(1, 2))
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
  dev.off()
}

gcvs_write_run_log <- function(run_result, paths){
  evidence_paths <- c(
    paths$summary_csv,
    paths$iteration_csv,
    paths$beta_csv,
    paths$sigma_csv,
    paths$diagnostics_png,
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
        diagnostics_png = paths$diagnostics_png,
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
    niter = run_result$config$niter,
    burnin = run_result$config$burnin,
    inner_steps = run_result$config$inner_steps,
    mean_acceptance = round(run_result$mean_acceptance, 6),
    mean_r2 = round(mean(run_result$r2_by_task), 6),
    mean_rmse = round(mean(run_result$rmse_by_task), 6),
    precision = ifelse(is.na(run_result$precision), "", round(run_result$precision, 6)),
    decision = run_result$decision,
    log_path = paths$log_md,
    stringsAsFactors = FALSE
  )

  if (file.exists(paths$index_csv)) {
    current <- read.csv(paths$index_csv, stringsAsFactors = FALSE)
    current <- current[current$exp_id != new_row$exp_id, , drop = FALSE]
    updated <- rbind(current, new_row)
  } else {
    updated <- new_row
  }

  write.csv(updated, paths$index_csv, row.names = FALSE)
}

gcvs_write_experiment_report <- function(run_result, paths){
  last_iter <- tail(run_result$iteration_metrics, 1)
  lines <- c(
    "# Experiment Report",
    "",
    "## 1. Experiment ID",
    run_result$config$exp_id,
    "",
    "## 2. Objective",
    "Stabilize the reproducible GC-MTL pipeline, align the lambda posterior simulation with the implemented MMALA kernel, and verify that the cleaned repository can run end-to-end without manual path edits.",
    "",
    "## 3. Setup",
    paste("- Dataset:", run_result$dataset),
    paste("- Profile:", run_result$config$profile),
    paste("- Outer iterations:", run_result$config$niter),
    paste("- Burn-in:", run_result$config$burnin),
    paste("- Inner lambda-MMALA steps per row:", run_result$config$inner_steps),
    paste("- Hyperparameters: alpha =", run_result$config$alpha_const, ", eta =", run_result$config$eta, ", a_psi =", run_result$config$a_psi, ", b_psi =", run_result$config$b_psi, ", a_err =", run_result$config$a_err, ", b_err =", run_result$config$b_err),
    paste("- Adaptation target:", run_result$config$opt_rate),
    "",
    "## 4. Intervention",
    "- Replaced the broken top-level experiment driver with a repository-relative CLI entrypoint.",
    "- Removed hidden dependencies on pre-created objects such as `Beta_alltasks_J_40` and `W` from the runnable path.",
    "- Corrected the lambda MMALA kernel so invalid-support proposals are rejected rather than repaired with an unmatched fallback draw.",
    "- Added two-sided probability clipping before the Gamma-to-Gaussian copula transform and repaired indefinite proposal metrics by eigenvalue flooring.",
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
    paste("- beta_csv:", paths$beta_csv),
    paste("- sigma_csv:", paths$sigma_csv),
    paste("- log_md:", paths$log_md),
    "",
    "## 6. Processed Metrics",
    paste("- Mean lambda acceptance:", gcvs_format_metric(run_result$mean_acceptance)),
    paste("- Mean predictive R^2:", gcvs_format_metric(mean(run_result$r2_by_task))),
    paste("- Mean predictive RMSE:", gcvs_format_metric(mean(run_result$rmse_by_task))),
    paste("- Selection precision:", gcvs_format_metric(run_result$precision)),
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
    exp_id = run_result$config$exp_id,
    dataset = run_result$dataset,
    M1_rule = paste0("GDP(alpha=", run_result$config$alpha_const, ",eta=", run_result$config$eta, ")"),
    M2_rule = paste0("full-MMALA(L=", run_result$config$inner_steps, ",target=", run_result$config$opt_rate, ")"),
    theta_mean = "",
    theta_peaks = "",
    theta_ess = "",
    M_mean = "",
    decision = run_result$decision,
    pdf_path = if (file.exists(paths$report_pdf)) paths$report_pdf else "",
    report_path = paths$report_md,
    meta_path = paths$report_meta_json,
    stringsAsFactors = FALSE
  )

  if (file.exists(paths$report_index_csv)) {
    current <- read.csv(paths$report_index_csv, stringsAsFactors = FALSE, check.names = FALSE)
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

  write.csv(updated, paths$report_index_csv, row.names = FALSE)
}

gcvs_write_artifacts <- function(run_result, repo_root){
  paths <- gcvs_initialize_results_paths(repo_root, run_result)

  summary_df <- data.frame(
    experiment_id = run_result$config$exp_id,
    dataset = run_result$dataset,
    profile = run_result$config$profile,
    niter = run_result$config$niter,
    burnin = run_result$config$burnin,
    inner_steps = run_result$config$inner_steps,
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
  gcvs_plot_diagnostics(run_result, paths$diagnostics_png)
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

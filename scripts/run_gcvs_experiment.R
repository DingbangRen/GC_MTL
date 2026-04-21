#!/usr/bin/env Rscript

script_file <- if (!is.null(sys.frames()[[1]]$ofile)) sys.frames()[[1]]$ofile else "scripts/run_gcvs_experiment.R"
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)

source(file.path(repo_root, "R", "gcvs_repro.R"))
gcvs_run_cli()

#!/usr/bin/env Rscript

script_file <- if (!is.null(sys.frames()[[1]]$ofile)) sys.frames()[[1]]$ofile else "Our_MTL_GCVS.R"
repo_root <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))

source(file.path(repo_root, "R", "gcvs_repro.R"))
gcvs_run_cli()

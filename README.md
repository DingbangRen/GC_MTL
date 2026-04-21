# GC_MTL

Gaussian-copula structured shrinkage for sparse multi-task regression, with a cleaned reproducible entrypoint and a manuscript-aligned posterior simulation workflow.

## What is in this repository

This repository contains the research code for the GC-MTL model together with a cleaned execution path that can be run from the repository root without hand-editing absolute paths or pre-loading hidden workspace objects.

The cleaned path focuses on the core Bayesian sampler:

- coefficient updates for each task
- GDP local-scale updates
- Gaussian-copula dependence updates
- adaptive full MMALA updates for `lambda`
- experiment logging, diagnostics, and report generation

## What was fixed

Compared with the original committed scripts, the cleaned pipeline now:

- provides a reproducible command-line entrypoint
- removes the broken dependency on undefined objects in `Our_MTL_GCVS.R`
- fixes the `lambda` MMALA kernel so proposals leaving the positive support are rejected instead of silently redrawn from a different law
- clips the Gamma CDF on both tails before the probit transform, avoiding `qnorm(0)` and `qnorm(1)` failures
- repairs near-singular proposal metrics numerically before sampling
- writes structured experiment artifacts following the local AGENTS workflow

## Quick start

Run a short smoke test from the repository root:

```bash
Rscript scripts/run_gcvs_experiment.R --dataset synthetic --profile smoke
```

Run the longer paper-style profile:

```bash
Rscript scripts/run_gcvs_experiment.R --dataset synthetic --profile paper
```

The legacy root-level entrypoint is still available and forwards to the same cleaned pipeline:

```bash
Rscript Our_MTL_GCVS.R --dataset synthetic --profile smoke
```

## R dependencies

The cleaned experiment path relies on the packages used by the original sampler plus a small reporting layer. In particular, make sure the following are installed in your R environment:

- `dplyr`
- `LaplacesDemon`
- `psych`
- `truncnorm`
- `BayesLogit`
- `extras`
- `matrixcalc`
- `Matrix`
- `jsonlite`
- `CholWishart`

The report HTML export additionally uses a local `pandoc` installation when available.

## Supported datasets

- `synthetic`
- `sarcos`
- `isolet`

The cleaned driver loads committed `.RData` files from the repository for the real-data runs. The synthetic run is regenerated from `simulation data setting.R`.

## Repository layout

- `scripts/run_gcvs_experiment.R`: command-line experiment entrypoint
- `R/gcvs_repro.R`: cleaned orchestration layer, dataset loader, artifact writer, and report generator
- `GCVS_posteriors.R`: posterior updates and MMALA sampler
- `Our_MTL_GCVS.R`: compatibility wrapper to the cleaned pipeline
- `results/README.md`: artifact conventions
- `results/experiments/index.csv`: machine-readable experiment index
- `results/reports/index.csv`: human-readable report archive index

## Generated artifacts

Each experiment run writes:

- `results/experiments/<EXP_ID>/summary.csv`
- `results/experiments/<EXP_ID>/iteration_metrics.csv`
- `results/experiments/<EXP_ID>/posterior_beta_mean.csv`
- `results/experiments/<EXP_ID>/posterior_sigma_mean.csv`
- `results/experiments/<EXP_ID>/figures/sampling_diagnostics.png`
- `results/experiments/<EXP_ID>/metadata.json`
- `results/experiments/<EXP_ID>/experiment_report.md`
- `results/experiments/<EXP_ID>/experiment_report.html` when `pandoc` is available
- `logs/iteration_logs/YYYYMMDD__<EXP_ID>.md`

The repository tracks the artifact structure and indexes, while transient experiment outputs are ignored by `.gitignore`.

## Legacy helper scripts

The following files are still present because they document the broader research workflow, but they are not part of the cleaned core pipeline:

- `preprocessing real data.R`
- `RMTL_and_DMFS.R`
- `ARWUL_R_running.R`

Their hard-coded Windows paths have been replaced by environment-variable based configuration:

- `GCVS_RAW_DATA_DIR`
- `GCVS_DMFS_ROOT`
- `GCVS_NORMT3_ROOT`
- `GCVS_BKMTL_ROOT`
- `GCVS_ARMUL_ROOT`

If those variables are not defined, the legacy comparison blocks now fail explicitly instead of silently depending on a machine-specific layout.

## Reproducibility status

The cleaned core experiment path is reproducible from the repository root and has been validated with smoke runs. Full paper-level replication still requires running the longer `paper` profile and comparing the resulting figures against the manuscript outputs.

## Notes on the manuscript

The manuscript source used for the current revision is maintained outside this code repository. The posterior computation section has been rewritten locally to match the implemented sampler, and the theoretical analysis section has been revised to separate the idealized adaptive-MAP argument from the actual fixed-hyperparameter posterior sampler.

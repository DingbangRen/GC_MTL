# Results Convention

This repository stores reproducible experiment artifacts under `results/experiments/` and structured run logs under `logs/iteration_logs/`.

## Minimum evidence per run

Each experiment run should write:

- one machine-readable summary CSV
- one sampling diagnostic PNG
- one structured markdown log
- one experiment report in the AGENTS section order

## Reserved report archive directories

The following directories are reserved for future human-readable report exports:

- `results/reports/pdf/`
- `results/reports/meta/`

Current runs also write:

- `results/experiments/<EXP_ID>/experiment_report.md`
- `results/experiments/<EXP_ID>/experiment_report.html` when `pandoc` is available
- `results/reports/index.csv`
- `results/reports/meta/<REPORT_BASENAME>.json`

If a local PDF engine is installed later, the same report basename is reserved under `results/reports/pdf/`.

## Git tracking policy

The repository keeps the artifact folder structure and index files under version control, but ignores experiment-specific output folders, iteration logs, PDF exports, and external exchange caches. This keeps GitHub readable while preserving a stable on-disk layout for local runs.

# Results Convention

This repository stores reproducible experiment artifacts under `results/experiments/`, unified human-readable reports under `results/reports/`, and structured run logs under `logs/iteration_logs/`.

## Minimum evidence per run

Each experiment run should write:

- one machine-readable summary CSV
- one sampling diagnostic PNG
- one structured markdown log
- one experiment report in the AGENTS section order

## Unified report directories

All reviewable reports should be discoverable from the single index:

- `results/reports/index.csv`

Report storage is split by role, not by ad hoc naming:

- `results/reports/daily/` for daily synthesis reports that summarize experiment progress, code changes, algorithm changes, manuscript corrections, and the next plan
- `results/reports/artifacts/` for manuscript, theory, shutdown, or other special-purpose review bundles
- `results/reports/pdf/`
- `results/reports/meta/`

Current experiment runs also write:

- `results/experiments/<EXP_ID>/experiment_report.md`
- `results/experiments/<EXP_ID>/experiment_report.html` when `pandoc` is available
- `results/reports/meta/<REPORT_BASENAME>.json`

If a local PDF engine is installed later, the same report basename is reserved under `results/reports/pdf/`.

The legacy top-level folder `results/daily_reports/` is deprecated. Daily reports should now be placed under `results/reports/daily/` and registered through the unified report index.

## Git tracking policy

The repository keeps the artifact folder structure and index files under version control, but ignores experiment-specific output folders, iteration logs, PDF exports, and external exchange caches. This keeps GitHub readable while preserving a stable on-disk layout for local runs.

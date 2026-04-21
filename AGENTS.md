# GC_MTL Agent Workflow

This file is the local operating guide for reproducibility work in `GC_MTL`.
It is adapted from a stronger experiment-governance template, but rewritten for
this repository's actual posterior sampler, artifact layout, and manuscript
workflow.

## A) Research Iteration Protocol

### A1. Mandatory State Machine

Every non-trivial experiment run must follow this order:

1. `Problem`: define the concrete question, scope, and success/failure criteria.
2. `Hypothesis`: state a falsifiable mechanism or expectation.
3. `Plan`: specify the minimal run, variables, metrics, and stop rule.
4. `Run`: execute and capture artifacts.
5. `Evaluate`: compare outcomes to the declared criteria.
6. `Decide`: set `PASS` | `NEEDS-DEBUG` | `INCONCLUSIVE`.
7. `Log`: persist a structured Markdown run log.

If a step is missing, the run is not complete.

### A2. Required Run Fields

Each run log and experiment report should record:

- `problem_id`
- `hypothesis_id`
- `experiment_id`
- `dataset`
- `profile`
- `assumption_list`
- `changes_made`
- `metrics_checked`
- `evidence_paths`
- `result_summary`
- `decision`
- `next_action`
- `risk_or_blocker`

### A3. Hard Gates

- No falsifiable hypothesis: do not run.
- No declared metric or criterion: do not run.
- No saved evidence path: do not make a final claim.
- Conflicting evidence: set `decision = INCONCLUSIVE`.
- If code behavior and manuscript wording disagree, treat the run as
  `INCONCLUSIVE` until one side is corrected.

### A4. Minimum Evidence Bundle

Every meaningful experiment should save at least:

- one machine-readable summary file, usually `summary.csv`
- one iteration-level diagnostic file, usually `iteration_metrics.csv`
- one visual diagnostic, usually `figures/sampling_diagnostics.png`
- one structured Markdown run log in `logs/iteration_logs/`

### A5. GC_MTL Posterior Diagnostics Focus

For posterior-sampling runs in this repository, explicitly track:

- `lambda` acceptance behavior
- `lambda` step-size adaptation behavior
- `lambda` lower-tail stability, especially very small positive values
- `Sigma` trace behavior and whether the PX update stays numerically stable
- retained-draw ESS when available, especially for `lambda` and `Sigma`
- posterior mean beta heatmap or support-recovery summaries on synthetic data
- task-correlation structure implied by posterior `Sigma`
- predictive metrics on real data
- whether the manuscript's Chapter 3 posterior-computation story matches the
  implemented kernel exactly

### A6. Cadence

- Per experiment: one pre-run plan and one post-run evaluation log.
- Per major milestone: one concise repository-level report.
- When touching manuscript-facing methodology: one explicit code-to-paper sync
  note.

## B) Artifact Storage Rules

### B1. Experiment Directory

Each run should write into:

- `results/experiments/<EXPERIMENT_ID>/`

Required core artifacts for the cleaned pipeline:

- `summary.csv`
- `iteration_metrics.csv`
- `posterior_beta_mean.csv`
- `posterior_sigma_mean.csv`
- `metadata.json`
- `experiment_report.md`

Recommended additional artifacts for longer or debugging-heavy runs:

- `figures/sampling_diagnostics.png`
- retained `lambda` draws
- retained `Sigma` draws
- ESS summary tables
- heatmap or task-correlation comparison figures
- `progress.json`
- `checkpoints/checkpoint_latest.rds`

### B2. Structured Log Location

Every run should append or create one Markdown log in:

- `logs/iteration_logs/YYYYMMDD__<EXPERIMENT_ID>.md`

The log should be short, checkable, and written after the run decision is made.

### B3. Report Archive Location

All reviewable reports should be discoverable from a single index:

- `results/reports/index.csv`

Use these subdirectories by report role:

- `results/reports/daily/` for daily synthesis reports that summarize
  experiment progress, code changes, key algorithm updates, manuscript
  corrections, and the next plan
- `results/reports/artifacts/` for manuscript, theory, shutdown, or other
  special-purpose review bundles

Archived PDFs, when generated, belong in:

- `results/reports/pdf/`

Associated metadata JSON belongs in:

- `results/reports/meta/`

### B4. PDF Filename Pattern

When a readable PDF is archived, use:

- `YYYYMMDD__EXPID__DATASET__PRIOR__SAMPLER__STATUS.pdf`

Example:

- `20260421__EXP-GCVS-SYNTHETIC-PAPER-001__SYNTHETIC__GDP-a3-eta1__FULL-MMALA__INCONCLUSIVE.pdf`

### B5. Index Files

Maintain these repository-level indexes:

- `results/experiments/index.csv`
- `results/reports/index.csv`

The experiment index should point back to the run log and key metrics.
The report index is the single catalog for experiment, daily, manuscript, and
handoff reports.

If a PDF is archived, create a same-basename metadata JSON in
`results/reports/meta/` linking back to:

- `summary.csv`
- `iteration_metrics.csv`
- `experiment_report.md`
- archived PDF path

Do not create a parallel top-level `results/daily_reports/` workflow. If a
daily report is generated, place it under `results/reports/daily/` and register
it in `results/reports/index.csv`.

### B6. Data and Path Rules

- Prefer repository-relative paths in runnable scripts.
- Avoid hard-coded user-specific or OS-specific paths in experiment code.
- If a manuscript file outside the repository is updated, record its absolute
  path in the run log and note that it is local-only until separately versioned.

### B7. Checkpoint Rules

- Long runs must write `results/experiments/<EXPERIMENT_ID>/progress.json`.
- Long runs must keep a latest resumable checkpoint at:
  `results/experiments/<EXPERIMENT_ID>/checkpoints/checkpoint_latest.rds`
- Archived checkpoints may be kept at regular iteration milestones in the same
  `checkpoints/` directory.
- `checkpoint_latest.rds` is for direct resume; archived checkpoints are for
  fallback and debugging.
- Checkpoints are local-first artifacts. Sync them to GitHub only when
  explicitly needed, typically with `git add -f`.

## C) Experiment Report Template

All experiment-readable reports should follow this section order:

1. `Experiment ID`
2. `Objective`
3. `Setup`
4. `Intervention`
5. `Raw Evidence`
6. `Processed Metrics`
7. `Key Findings`
8. `Mechanism Hypothesis`
9. `Alternative Explanations`
10. `Comparison`
11. `Failure / Limitation`
12. `Next Step`
13. `Confidence`

### C1. Language Rule

Experiment reports should be written in English only.
This is mainly to reduce PDF and cross-platform rendering problems and to keep
the report format stable across environments.

### C2. Raw Evidence Constraint

Section `[5] Raw Evidence` should contain checkable fragments only:

- raw numeric metrics
- direct artifact paths
- short log snippets when genuinely useful
- data samples when relevant

Do not interpret results in this section.

### C3. Processed Metrics for GC_MTL

For this repository, `[6] Processed Metrics` should normally include the subset
that applies:

- mean `lambda` acceptance
- acceptance range across rows or iterations
- mean step size
- `lambda_min` or another lower-tail stability summary
- `Sigma` trace summaries
- ESS summaries for retained `lambda` and `Sigma` draws
- support-recovery metrics on synthetic data
- predictive `R^2` / RMSE on real data
- final run decision

## C4. Daily Synthesis Report Template

When the report is a daily or milestone synthesis rather than a single run,
reuse the experiment-template discipline but include these sections explicitly:

1. `Date / Scope`
2. `Experiments Advanced Today`
3. `Code Changes`
4. `Algorithm Changes`
5. `Manuscript Revisions`
6. `Code-to-Paper Alignment`
7. `Key Findings`
8. `Open Risks`
9. `Next Plan`

The daily report should still cite concrete evidence paths and should register
into the same `results/reports/index.csv` catalog as all other report types.

## D) Decision Convention

Use one of:

- `PASS`: evidence supports the planned claim and diagnostics are acceptable for
  that claim.
- `NEEDS-DEBUG`: a major convergence, reproducibility, or logic criterion fails.
- `INCONCLUSIVE`: evidence is partial, conflicting, or too short to justify a
  stronger statement.

Default to `INCONCLUSIVE` if:

- the chain is too short for the claim being made
- ESS is not available for the parameters under discussion
- figure-level agreement with the manuscript is still qualitative only
- the run is numerically stable but methodological alignment is still under
  review

## E) Manuscript Sync Rules

When the code and manuscript interact, follow these rules:

- Do not describe a sampler step in the manuscript unless it matches the code.
- If the implemented kernel changes, update the Chapter 3 posterior-computation
  section and log the change.
- If synthetic or real-data claims are tightened or weakened, cite the exact
  artifact path that justifies the change.
- Record manuscript edits as local-only when the manuscript lives outside the
  repository.

## F) Current Runnable Entry Points

Current cleaned entry points include:

- `Rscript scripts/run_gcvs_experiment.R --dataset synthetic --profile smoke`
- `Rscript scripts/run_gcvs_experiment.R --dataset synthetic --profile paper`
- `Rscript scripts/run_gcvs_experiment.R --dataset synthetic --profile paper --exp-id <EXP_ID> --resume`
- `Rscript Our_MTL_GCVS.R`
- `scripts/start_codex_tmux_supervisor.sh`
- `scripts/prepare_gcvs_handoff.sh <EXP_ID> [TMUX_SESSION]`

These are the preferred starting points for reproducibility work. Legacy helper
scripts may still exist for data preparation or external baselines, but they
should not be treated as the primary reproducible entrypoint unless their path
and dependency assumptions are first cleaned and logged.

# Reproducibility Summary — 2026-04-21

## Synthetic paper profile

- Final run: `EXP-GCVS-SYNTHETIC-PAPER-FIXED-20260421`
- Configuration: `niter = 300`, `burnin = 150`, `inner_steps = 40`, `sigma_px_steps = 1`
- Key code-path correction: the lambda MMALA block now warm-starts from the previous outer-iteration `Lambda`, rather than repeatedly reinitializing rows from fresh copula draws.

Observed behavior:

- Posterior mean coefficient recovery is qualitatively strong.
- Task-wise absolute correlations between posterior mean `|Beta|` and truth are high: roughly `0.87` to `0.99`.
- Row-wise lambda acceptance is stable but somewhat below the nominal target: mean row acceptance ranges about `0.506` to `0.560`.
- Representative lambda ESS is mixed: approximately `5` to `30`.
- Representative Sigma ESS remains weak: approximately `3` to `8`.

Alignment call:

- The coefficient heatmap is broadly aligned with the manuscript description.
- The Sigma correlation summary is only partially aligned.
- Tasks `1`, `2`, and `3` are strongly positively correlated, and task `7` is an outlier with negative correlations to the main block.
- Task `8` does not match the manuscript-era synthetic correlation story: in the cleaned paper-profile run it stays positively correlated with tasks `1` to `6` and negatively correlated mainly with task `7`.

Conclusion:

- The lambda-mixing path is materially improved relative to the earlier cleaned baseline.
- The synthetic paper-profile run is qualitatively closer to the manuscript, but the Sigma posterior is not converged strongly enough to justify tightening manuscript wording yet.

## Long-chain check

- A resumable long-chain checkpoint exists at `EXP-GCVS-SYNTHETIC-PAPER-LONG-20260421`.
- The early long-chain warmup moved into the same improved qualitative regime as the fixed paper-profile run.
- A full long-chain retained-draw assessment is still pending.

## Real-data reproducibility

- Completed run: `EXP-GCVS-SARCOS-SMOKE-20260421`
- Result: the cleaned SARCOS path runs end-to-end and writes the expected artifacts, checkpoints, report, and index entries.
- Smoke-length posterior structure is directionally plausible but not yet aligned tightly enough with the manuscript’s two-block task-correlation description to count as a full reproduction.

## Manuscript status

- `main.tex` was not updated in this pass.
- Reason: the synthetic coefficient story is acceptable, but the Sigma posterior summary and ESS do not yet support stronger manuscript claims.

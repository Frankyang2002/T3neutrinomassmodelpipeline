# Reduced-scalar frozen UV: existing interactive RGE and sensitivity oval reports

This is the exact report family described by `docs/FERMION_START_REUSE_EXISTING_REPORTS.md`, **not** the general-scalar 48-case HTML/PDF atlas and **not** the lambda-style summary/print adapter.

## Established HTML generators reused unchanged

- RGE: `Numerical.plotting.InteractiveModelComparison.write_dashboard()`.
- Ovals: `Numerical.plotting.BuildFittedModelEnvelopes.build()`.

`BuildReducedScalarExistingInteractiveReports.py` prepares input data for those two established renderers. No source renderer changes.

## What is new relative to the saved final report

The completed `reduced_frozen_overnight_v2` campaign records only low-energy endpoints for 2 scenarios x 2 Yukawa modes x 16 models. **It does not contain gamma scans or 41-point Weinberg RGE trajectories.** These must be evaluated with the *same* frozen-UV reduced-scalar prescription used by the Sobol runner before the HTML can be rendered. This script computes those **diagnostics**, but never optimises Yukawas or reruns symbolic matching.

For T1/T2/T3, the sensitivity scan scales the input lambdaT3 and the corresponding Yukawa matrices as in `Numerical.fitting.BenchmarkSensitivity`, and **re-evolves the scalar EFT1 when lambdaT3 changes**. Gamma=1 is checked against the actual saved model/mode neutrino endpoint; numerical failures are marked failed and not interpolated. The original `BuildFittedModelEnvelopes` HTML defines descriptive empirical ovals, not confidence regions.

RGE views use the sampled *physical SM+Weinberg* evolution from `muS=1e3 GeV` to `muLow=1e2 GeV`, matched using the original fixed-order C5 evaluator. The scalar EFT1 segment from `muF=1e5 GeV` to `muS` is included in the matched endpoint, **but no intermediate Wilson C5 trajectory is invented**. UV running is disabled; there is no UV RGE curve. Four legacy C5 source schemes remain nonuniform. Additional scalar quadratics are zero at input and not evolved.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -m pytest -q tests\python\test_reduced_scalar_existing_interactive_reports.py

# First test one model in one case (separate pilot output).
python -u -m scripts.BuildReducedScalarExistingInteractiveReports `
  --case small_uv_frozen_5e-9 --mode optimised --model A_-2 --kind all

# Complete: both frozen scenarios x both modes x all 16 models.
python -u -m scripts.BuildReducedScalarExistingInteractiveReports --kind all `
  2>&1 | Tee-Object -FilePath reduced_interactive_reports.log
```

To run a faster first pass, add `--kind rge` (no sensitivity evaluations). To build the ovals later, run with `--kind sensitivity`. The full scan is computational work, about 2 x 2 x 16 x 3 x 11 = 2,112 sensitivity evaluations, plus RGE samples; **runtime has not been benchmarked**. Rerun the same command to reuse completed saved diagnostics; `--force` discards that reuse. `--gamma-min 0.1 --gamma-max 10 --points 11` is default, as in the original documented pilot.

## Outputs

Full results:

```text
Reports/output/display_reduced_scalar/reduced_frozen_overnight_v2/
  small_uv_frozen_5e-9/reused/optimised/interactive_rge.html
  small_uv_frozen_5e-9/reused/optimised/sensitivity_ovals.html
  small_uv_frozen_5e-9/reused/shared_a_m2/interactive_rge.html
  small_uv_frozen_5e-9/reused/shared_a_m2/sensitivity_ovals.html
  large_uv_frozen_0p5/reused/optimised/interactive_rge.html
  large_uv_frozen_0p5/reused/optimised/sensitivity_ovals.html
  large_uv_frozen_0p5/reused/shared_a_m2/interactive_rge.html
  large_uv_frozen_0p5/reused/shared_a_m2/sensitivity_ovals.html
```

Saved, reproducible diagnostic sources:

```text
output/studies/reduced_scalar_interactive_diagnostics/reduced_frozen_overnight_v2/
  <case>/<mode>/models/<model>/running_diagnostics.json
  <case>/<mode>/sensitivity/<model>/sensitivity_diagnostics.json
  <case>/<mode>/sensitivity/comparison_summary.json
```

A single-model pilot writes to `.../pilots/<model>/<case>/<mode>/` in **both** roots so it does not impersonate an all-model report. All original lambda-study report directories remain untouched.

This adapter is local-only. It should be run from the same checkout containing the completed Sobol fits and native RGBeta/C5 artifacts. It has not been physically executed on the user's Windows repository; unit tests do not certify the upstream C5 conventions.

## RGE integration failures at large gamma (29 October 2026 correction)

The initial report adapter allowed the known high-coupling EFT1 integration
failures during `ReducedEvaluator.__init__` to abort the *entire* model scan,
even though `ReducedEvaluator.evaluate()` already returns `NumericalFailure`
for the same failures occurring later in a trial. This happens for the
`large_uv_frozen_0p5` sensitivity scan (e.g. `Required step size is less than
spacing between numbers`). This is a report sensitivity gamma-point issue, not
a failure of the previously accepted gamma=1 fitted benchmark.

The corrected adapter records these **specific** numerical failures and
nonpositive running scalar thresholds as `NumericalFailure` points, with their
original exception messages. It reuses a cached failed EFT1 boundary for the
same lambdaT3 where applicable. It does **not** interpolate or replace failed
points, modify the solver, or widen perturbativity bounds. Metadata, beta
inventory and other unrecognised errors still stop the run. All three gamma=1
centres must pass their original neutrino endpoint checks.

The unchanged oval renderer requires at least three successful gamma points per
T1/T2/T3 direction, including gamma=1. Where there are too few valid points,
that model is **excluded from the oval dashboard**, while its source JSON
retains failed points and reasons. Such an exclusion is not model exclusion.
The original wide scan may therefore show fewer than 16 models (in the first
small-coupling `optimised` run it showed 14, excluded 2).

To resume from the partially completed all-model run, reinstall this package
and execute the original command **without `--force`**. Existing source-hash-
verified RGE trajectories and sensitivity scans are reused; the twelve large-
coupling optimised models whose scans were aborted are computed afresh. No
frozen fits are modified.

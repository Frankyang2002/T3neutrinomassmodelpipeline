# T3 high scalar (0.5) — fermion-threshold-start reproduction

**Install:** `Numerical/orchestration/RunHighScalarHalfFermionStart.py` (one new module only).  **Use locally; do not push to GitHub.**

Reference source inspected: GitHub `main` commit `9abe36f583e3cb2d7a3d1ce709032af7bb357a7f` (24 Oct 2026).  This is an *adapter* to the repository's already implemented `RunFermionStartLambdaStudy.py`, `BuildFermionStartRunningReports.py`, and `BuildFermionStartExistingReports.py`; it does not replace those files.

## Physics convention — UV FROZEN

The benchmark parameter values are imposed **at the fermion matching scale** `mu_F`. There is no UV RGE integration. The verified evaluation sequence in `RunFermionStartLambdaStudy.FermionStartEvaluator` is:

1. `frozen_result(initial, mu_f, uv=True)`: copies the UV input state to `mu_F` without beta-function evaluations;
2. `build_intermediate_scalar_boundary`: integrates out the heavy fermions at `mu_F`;
3. `run_intermediate_scalar_segment`: evolves the scalar-only EFT1 from `mu_F` to `mu_S`;
4. `FinalC5TrajectoryEvaluator`: matches the Weinberg coefficient at the scalar threshold;
5. `continue_with_weinberg_running`: evolves the SM+Weinberg EFT to `mu_low`.

Numerical scalar coupling fixed during fitting: `lambdaT3=0.5` in classes A–D, **0.25** in class E (preserves the original special Class E convention). Yukawa entries `y1` and `y2` are refitted; fermion/scalar masses and threshold scales are not reoptimised. The rescaled initial Yukawas derive from existing **5e-9 individually fitted** input files in `output/studies/lambda_5e-9_refit/fitted/`; these are *seeds*, not accepted results at 0.5.

The script passes through the existing 18-real-parameter Sobol and ranked 5/8/12-active-parameter least-squares fit. It uses a smaller Sobol neighbourhood: `--sobol-fraction 0.6`, `--sobol-floor 1e-5`, `--bound 0.01`. The old 0.1 floor is too large for the O(1e-4) scaled Yukawa seeds. The chi-squared acceptance cutoff remains **5.86** as in the previous study.

`T3-A-m2` is the donor for `shared_a_m2`; importantly, it is the **new 0.5 fit** and not the saved low-scalar or UV-start A(-2) matrices.

These standard reports use the existing fixed-order EFT implementation. They do **not** automatically evaluate the exploratory *coupled* C12–C5 solver, and fitted oscillation observables are calibration targets, not independent physics predictions.

## Local installation and complete PowerShell commands

Run everything from the repository root, with your usual project's Python environment enabled:

```powershell
# Save T3_HIGH_SCALAR_HALF_FERMION_START.zip somewhere accessible.
# Verify working tree before extracting a new module.
git status --short
Expand-Archive -LiteralPath .\T3_HIGH_SCALAR_HALF_FERMION_START.zip -DestinationPath . -Force
python -m Numerical.orchestration.RunHighScalarHalfFermionStart --help

# Inspect old fitted seeds and per-model symbolic matching/RGE artifacts.
python -m Numerical.orchestration.RunHighScalarHalfFermionStart preflight --check-artifacts

# 1. Fermion-start, UV-frozen inputs (all 16 models).
python -m Numerical.orchestration.RunHighScalarHalfFermionStart prepare

# 2. Fit the donor first (256 Sobol candidates, rescalings + local refinement).
python -m Numerical.orchestration.RunHighScalarHalfFermionStart fit --model T3-A-m2 --sobol-points 256 --resume

# 3. Fit the other 15 models and retain the successful donor fit.
python -m Numerical.orchestration.RunHighScalarHalfFermionStart fit --sobol-points 256 --resume

# 4. Compare own fitted Yukawas vs the new A(-2) shared Yukawas.
python -m Numerical.orchestration.RunHighScalarHalfFermionStart compare --mode both

# 5. Endpoint HTML, PNG comparison graphs, and central overview page.
python -m Numerical.orchestration.RunHighScalarHalfFermionStart display

# 6. Actual EFT1 / Weinberg trajectories, CSVs and interactive SVG HTML.
python -m Numerical.orchestration.RunHighScalarHalfFermionStart running --mode both --samples 15

# 7. Existing interactive-RGE and sensitivity-oval HTML generators.
python -m Numerical.orchestration.RunHighScalarHalfFermionStart reports --mode both --kind all

# View results.
Start-Process 'Reports\output\display_fermion_start\0p5\overview.html'
```

**Single-run alternative:** `python -m Numerical.orchestration.RunHighScalarHalfFermionStart all --resume --check-artifacts` runs the same stages sequentially. This is a significant optimisation + trajectory study; it is not the single fixed-Yukawa `RunT3BHighScalarHalf.py` analysis.

**Checkpoint note:** Keep `--sobol-points` and all Sobol range arguments identical when using `--resume`. A change in those parameters invalidates the saved checkpoint signature. If you pilot at 64 points, move the previous model's `sobol/<model>.json` out of the study folder before repeating at 256. Do not overwrite the original `5e-9` or `5e-10` study paths.

If your local artifact directory layout differs, pass `--artifact-map path\to\artifact_map.json` at preflight and at every numerical stage (`fit`, `compare`, `running`, `reports`). Its entries use the pre-existing keys `uv_rgbeta`, `eft1_rgbeta` and `final_weinberg`. UV beta expressions may be *located* for compatibility, but **UV beta functions are never numerically integrated** by this evaluator.

## Intended output locations (relative to project root)

| Data | Folder/file |
|---|---|
| Prepared 0.5 inputs | `output/studies/fermion_start_lambda_comparison/0p5/prepared/` |
| Sobol checkpoints | `output/studies/fermion_start_lambda_comparison/0p5/sobol/` |
| Fitted 0.5 Yukawas and fit status | `output/studies/fermion_start_lambda_comparison/0p5/fitted/` and `fits/` |
| Optimised and shared endpoints | `output/studies/fermion_start_lambda_comparison/0p5/comparison/{optimised,shared_a_m2}/` |
| Existing report-ready numeric diagnostics | `output/studies/fermion_start_lambda_comparison/0p5/reports_existing/` |
| Perturbativity magnitude screen | `output/studies/fermion_start_lambda_comparison/0p5/eft1_perturbativity_screen.json` |
| Case overview HTML | `Reports/output/display_fermion_start/0p5/overview.html` |
| Endpoint figures | `Reports/output/display_fermion_start/0p5/figures/` |
| True EFT1 + Weinberg trajectories and CSV | `Reports/output/display_fermion_start/0p5/running/{optimised,shared_a_m2}/` |
| Pre-existing renderer HTML outputs | `Reports/output/display_fermion_start/0p5/reused/{optimised,shared_a_m2}/` |

Perturbativity is screened by the largest absolute scalar quartic captured over sampled EFT1 running, against `4*pi`. This is a **necessary practical magnitude check only**; it does not prove potential stability, unitarity, convergence or validity of large-gamma variations. The existing sensitivity directions remain T1/T2/T3 from `BenchmarkSensitivity.py`; default gamma range here is `0.5..2` to keep the first large-coupling inspection contained. To reproduce the original wider display range, add `--gamma-min 0.1 --gamma-max 10 --force-reports` to the *reports* command; review failures and perturbativity before interpreting it.

## Validation scope

The new module has been syntax-compiled in isolation. It has **not** been numerically executed against your local matching artifacts in this environment; no 0.5 fit values or plots are claimed to exist until you run the commands above and inspect the outputs. The unchanged report generators themselves were verified against current GitHub source, not against a local modified checkout. Do not silently replace locally modified versions with the GitHub files.

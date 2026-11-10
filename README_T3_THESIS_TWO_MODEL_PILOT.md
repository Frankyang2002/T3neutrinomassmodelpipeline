# A(-4) and D(-2) thesis-benchmark pilot — v3 Sobol + least squares

> **Current source audit — 10 November 2026:** This pilot fits a **restricted conditional** one-loop EFT evaluator, **not** the complete general-Lagrangian one-loop two-threshold calculation. Its runtime enters `scripts/RunT3FHardTriangleBatch1.py` with `builder.Options(reduced_potential=True)` and includes only the analytic F-hard triangle plus the specified scalar bubble at MS. It does **not** invoke `scripts/RunT3CompleteGeneralScalarSymbolic.py`, `scripts/RunT3FullScalarContinuous.py`, or the separate general Matchete one-loop branches. A low chi-square would be a fit within the restricted approximation, not the full-general result requested. For the source-level feasibility determination and model-specific A(-4)/D(-2) limitations, see `docs/IMPLEMENTATION_STATUS.md` and `docs/T3_A4_D2_GENERAL_PIPELINE_FEASIBILITY_2026_11_10.md`. **Do not start a large optimisation under the assumption that it uses the complete general Lagrangian.**


## What changed

This is the **complete replacement** of the v2 pilot driver, its tests, and this README. The thesis benchmark configuration and all physics modules are unchanged.

Two problems in v2 were corrected:

1. **Incorrect exit status:** the script included a previously blocked *other* model in the shell exit-code decision. A successful `--baseline --model D_-2` could print `COMPLETE_CONDITIONAL` and nevertheless exit with code 2 because a historic A(-4) record was blocked. In v3, the overview still preserves both model rows, but the exit status depends only on the model(s) requested in this invocation.
2. **Not actually least-squares fitting:** v2 `--fit` ran a few Sobol points followed by random-coordinate descent. V3 replaces that fit routine with scrambled **Sobol samples across all 18 real Yukawa entries**, sensitivity checks, and **bounded `scipy.optimize.least_squares`** refinement of the most sensitive entries. No other physics or matching formulas are changed.

The source-model evaluator continues to recompute UV RG, all-fermion matching, the analytic hard triangle, strict seven-channel EFT1 running, simultaneous scalar matching, and low-energy SMEFT running at 100 GeV for **every new trial**. Successful candidate evaluations from v2 are reusable because their physics inputs and evaluator were not changed. Previously blocked candidates are retried. Each new candidate has its own isolated scratch data. The historical all-16 Batch 8 outputs remain unchanged.

## Physics and benchmark

- **Models:** A(-4) and D(-2) only.
- **Thesis benchmark:** defined in the existing `configs/thesis_numerics/t3_thesis_A4_Dm2_pilot.json`; `mu_UV=1e7`, `M_F=1e5`, `M_S=1e3`, `mu_low=100` GeV, `lambdaT3=0.1`, other source scalar couplings `0.01`, both original real 3x3 Yukawa matrices as starting points, no class-E half multiplier.
- **Prescription fitted:** `analytic_triangle` with the current strict EFT1 and scalar-matching calculation. `native_formal` and `native_as_UV_minus_soft_hypothesis` are **not** refitted by this pilot.
- **Five Gaussianised targets:** `configs/nufit_6_1_ic24_no.json`. An original-thesis Yukawa seed with high chi-square is a **valid numerical prediction but a poor fit**, not an execution error.
- **Scope:** conditional incomplete one-loop matching; the existing reduced source tensor omits unmapped additional scalar quartics. Do not describe these results as a fully certified physical MSbar matching calculation. The original Workstream B SM `yu/yd/ye` defaults are still used because Table 4.1 alone did not define a full substitute.

## Installation

Extract `t3_thesis_A4_Dm2_pilot_sobol_ls_v3.zip` **into the repository root**, replacing the three files below:

- `scripts/RunT3ThesisTwoModelPilot.py`
- `tests/test_t3_thesis_two_model_pilot.py`
- `README_T3_THESIS_TWO_MODEL_PILOT.md`

All paths are relative to `C:\thesistemp\T3neutrinomassmodelpipeline`. No GitHub changes.

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive `
  -LiteralPath "$HOME\Downloads\t3_thesis_A4_Dm2_pilot_sobol_ls_v3.zip" `
  -DestinationPath . -Force

.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_thesis_two_model_pilot.py
```

## Run the two baselines

The **D(-2) numeric calculation already completed** in your previous log, despite the shell throwing an error. Its conditional strict-one-loop masses were

`[0.0008061911, 0.0150880588, 0.0869668958] eV`, with an approximate five-observable `chi2=85161.4086`.

After installing v3, the same input is reused from its successful evaluation cache, and the exit status should be zero even if the previous A(-4) result is blocked:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3ThesisTwoModelPilot.py --baseline --model D_-2
if ($LASTEXITCODE -ne 0) { throw "D_-2 baseline failed" }

.\.venv\Scripts\python.exe scripts\RunT3ThesisTwoModelPilot.py --baseline --model A_-4
if ($LASTEXITCODE -ne 0) { throw "A_-4 baseline failed" }
```

## Run the actual Sobol + least-squares fits

```powershell
.\.venv\Scripts\python.exe scripts\RunT3ThesisTwoModelPilot.py `
  --fit --model D_-2 --max-evals 80 --sobol-points 16 `
  --screen-points 18 --active-params 5 --starts 1 --local-max-nfev 12

if ($LASTEXITCODE -ne 0) { throw "D_-2 Sobol+least-squares pilot incomplete" }

.\.venv\Scripts\python.exe scripts\RunT3ThesisTwoModelPilot.py `
  --fit --model A_-4 --max-evals 80 --sobol-points 16 `
  --screen-points 18 --active-params 5 --starts 1 --local-max-nfev 12

if ($LASTEXITCODE -ne 0) { throw "A_-4 Sobol+least-squares pilot incomplete" }
```

The budget counts only **new** full EFT evaluations; successful cached points cost no new runs. The scrambled Sobol points search the full 18-dimensional bounded region around the thesis seed. The sensitivity stage perturbs up to 18 coordinates at the best sampled candidate. The default local least-squares search changes the five highest-sensitivity entries; the other 13 are held at their selected Sobol start values during that local solve. This is the first computationally manageable numerical fit, not a global search proof. To fit a larger subset later, increase `--active-params`, `--max-evals`, and possibly `--starts`.

The old `--max-evals 12` settings are intentionally rejected in v3 for the default Sobol + screening + least-squares configuration: they do not permit one sufficiently evaluated local gradient. The script reports the minimum required budget. `--perturbation` controls the bounded trial region around the old fitted Yukawa entries (default `1.0` times `max(abs(seed entry), 2e-5)` for each coordinate).

## Outputs

Per model (`A_m4` or `D_m2`):

- `output/parallel_t3/thesis_A4_Dm2_pilot_v1/<tag>/pilot_summary.json` — best/seed, Sobol evaluations, sensitivity ranking, least-squares run metadata, conditional status
- `output/parallel_t3/thesis_A4_Dm2_pilot_v1/<tag>/best_yukawas.json` — best real `Y1_UV`, `Y2_UV` matrices and chi-square
- `output/parallel_t3/thesis_A4_Dm2_pilot_v1/<tag>/evaluations/<SHA256>.json` — saved individual candidate, prediction, masses, residuals, source code classification
- `output/parallel_t3/thesis_A4_Dm2_pilot_v1/two_model_overview.csv` — one row for each previously attempted model, even when run on separate commands

All fitted results are `CONDITIONAL_NOT_PHYSICAL_MSbar`. `PILOT_BEST_FOUND` means **best verified candidate among those evaluated**, not proof of a unique/global optimum; check `fit_acceptable_project_cut`, the per-observable residuals, and `least_squares_runs` before using a fit in a thesis table.

## Tests/status

The eight isolated tests passed for the benchmark, UV source-key handling, failed/cached candidate behavior, selected-model exit code, and synthetic Sobol-to-least-squares improvement. These do **not** constitute an end-to-end test of your Windows physics repository. No new real A(-4) or fitted D(-2) numerical results were generated during preparation of this bundle.

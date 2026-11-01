# T3 thesis: per-model Sobol + least-squares numerical fit

## Physics scope

This runs the repository's **existing restricted/ordinary-T3** UV -> F -> EFT1 -> scalar -> Weinberg pipeline using the saved `output/full/comparison/<scenario>/.../data` JSON artifacts. It does not implement complete full-general-quartic threshold matching. The scalar parameter under study is `ordinary.lambdaT3` (not the full scalar-quartic basis).

Each of the 16 configurations has its **own** Yukawa fit. The objective is the existing Gaussian approximation to five NuFIT 6.1 NO observables. Fitted values are calibration targets, not independent predictions. The scalar quartics other than `lambdaT3`, gauge inputs, masses, and scales remain fixed to the existing per-model reference config.

## Starting inputs

Source: `Numerical/benchmarks/fitted_reference/<model>.json` and `benchmark_manifest.json` in the GitHub repository. The defaults use `lambdaT3=0.1` in A-D, `lambdaT3=0.05` in E, and the original full 3x3 real `y1`, `y2` matrices as the seed. The class-E half coupling corresponds to the original comparison normalisation convention, not a newly established full-theory equivalence.

For a requested alternative `--lambda-t3 0.5`, models A-D use `0.5` and E uses `0.25` by default (disable that with `--no-class-e-half`). Initial Yukawas are multiplied by `sqrt(lambda_reference/lambda_new)` as a *seed heuristic* only; the actual fit and diagnostic rerun the full existing numerical evaluator.

## Algorithm

1. Find matching saved `uv_rgbeta_rge.json`, `eft1_rgbeta_rge.json` and `final_weinberg_coefficient.json` from the restricted comparison archive.
2. Fix `lambdaT3` at the requested value (unless `--fit-lambda` is explicitly used).
3. Generate **128 scrambled Sobol points** in the 18 real entries of `y1`, `y2`. Each entry has a model-specific range centred on the rescaled fitted reference.
4. Retain the previous fitted-reference seed as an explicit candidate; never discard it in favour of a worse Sobol point.
5. Run the existing reduced, bounded least-squares refinement (sensitivity-ranked 5 -> 8 -> 12 free parameters). If the fit target has already been met, the runner still attempts one 5-parameter local confirmation by default; it never replaces a better fit with a worse confirmation.
6. Save best full input JSON, score, predictions, numerical running diagnostics and gamma-direction sensitivity data. Generate HTML via **existing** `InteractiveModelComparison` and `BuildFittedModelEnvelopes` renderers.

The search is local to the supplied bounds; it does not prove a unique global optimum. A lowest `chi2` near zero is a fitted calibration point, not evidence the model predicted five oscillation observables a priori.

## Install

Extract the ZIP at the root of `C:\thesistemp\T3neutrinomassmodelpipeline`.

- `scripts/RunThesisSobolModels.py` — new orchestration (full file).
- `scripts/BuildThesisNumericHTML.py` — same existing-report adapter as the earlier thesis numerical pilot; included for stand-alone installation.
- `tests/python/test_thesis_sobol_models.py` — pure config tests.

## Windows PowerShell commands

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive "$HOME\Downloads\T3_Thesis_Sobol_16_Model_Fits.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_thesis_sobol_models.py
python -X utf8 -m scripts.RunThesisSobolModels --list
```

**Start with one model:**

```powershell
python -X utf8 -m scripts.RunThesisSobolModels --model A_-2 --sobol-points 128 --run
```

**Then all remaining/all 16 (completed A(-2) will be reused):**

```powershell
python -X utf8 -m scripts.RunThesisSobolModels --sobol-points 128 --run
```

**Second scalar setting, new output directory:**

```powershell
python -X utf8 -m scripts.RunThesisSobolModels --lambda-t3 0.5 --sobol-points 256 --run
```

**Optional: optimize the scalar itself alongside the Yukawas:**

```powershell
python -X utf8 -m scripts.RunThesisSobolModels --model A_-2 --lambda-t3 0.1 --fit-lambda --tag fit_scalar_A_m2 --run
```

**If there is no saved `output/full/comparison/smallY_smallL` case for a model**, first generate the standard restricted comparison artifacts. This script deliberately fails with `Missing saved matching/RGE artifact` rather than silently claiming to use the full-general matching data.

## Outputs (example with default reference scalar)

- `configs/thesis_sobol/reference_scalar/<long_model_key>.json` — best fitted model input, including full matrices.
- `output/thesis_sobol/reference_scalar/<long_model_key>/fit_summary.json` — per-model `chi2`, scalar value, parameter map, prediction, Sobol and local refinement metadata.
- `output/thesis_sobol/reference_scalar/<long_model_key>/sobol_phase.json` — retained candidates and search metadata; separate from existing native `output/benchmarks/...` caches.
- `output/thesis_sobol/reference_scalar/<long_model_key>/running_diagnostics.json` — original trajectory and observable format.
- `Reports/output/thesis_sobol/reference_scalar/fit_summary.csv` — scalar values and fit statistics across the models.
- `Reports/output/thesis_sobol/reference_scalar/observables.csv` and `.json` — low-energy masses, mass splittings and mixing angles.
- `Reports/output/thesis_sobol/reference_scalar/interactive_RGE.html` — existing interactive running/observable renderer.
- `Reports/output/thesis_sobol/reference_scalar/sensitivity_oval.html` — existing covariance-oval renderer if successful gamma trajectories were generated.

By default, T1/T2/T3 use 11 gamma points over 0.8-1.2 including gamma=1, with the original sensitivity definitions and fit target. The ovals are descriptive covariance of the chosen parameter response, **not experimental confidence contours**.

## Practical notes

- This is an expensive numerical calculation: each Sobol point and each local/sensitivity step evaluates the production running and matching pipeline. Start with one model first.
- The fitter can skip deep least squares when its target has been met. The extra confirmation LS attempt may be disabled with `--no-force-ls`.
- `--no-sensitivity` runs fitting and trajectory evaluation without all gamma scans.
- `--force` replaces results only within this explicitly selected `--tag` study. To retain multiple settings, give each experiment a distinct `--tag`.
- Never treat a fitted-reference chi-squared as a new independent neutrino-mass prediction. Keep a common-input (unfitted) comparison alongside fitted per-model results when evaluating representation-dependent patterns.
- Code was syntax-checked and the 4 local pure-parameter tests passed; full physical numerical runs require your local saved artifacts and dependencies.

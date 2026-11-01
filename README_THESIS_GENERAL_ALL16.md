# Full-general RGE numerical pilot — 16 models

This is a **local-only exploratory campaign**. It reuses the repository's existing
canonical full-general one-loop UV/EFT1 RGEs, existing numerical solvers,
model-specific matching kernels, NuFIT scoring, Sobol and staged least squares.

**Scope:** B(-1)/C(-1) reuse the working general-quartic + mixed-bilinear
one-scale Matchete C5 calculation (`RunT3BCGeneralFromReference.py`). The
other 14 models evolve all canonical general renormalisable coupling beta keys
(including allowed additional scalar operators and conjugate Yukawas), but use
an existing original-topology leading-C5 matching kernel. Additional general
quadratic/Yukawa direct matching contributions are not universally included;
this is **not a full general sequential Wilson/EFT calculation**.

## Files and intended project paths

- `scripts/RunT3GeneralAll16FromReference.py` — full new campaign controller.
- `scripts/RunT3BCGeneralFromReference.py` — *identical complete script* previously supplied for the successful B(-1) run, included for installation reproducibility.
- `tests/python/test_thesis_general16.py` — local utility/CLI smoke tests.
- `README_THESIS_GENERAL_ALL16.md` — this guide, at the repository root.

None of these overwrite existing physics modules or archived result files.

## Installation

From PowerShell in `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Expand-Archive "$HOME\Downloads\T3_General_All16_Optimisation.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_thesis_general16.py
```

## 1. Inspect data availability before fitting

```powershell
python -X utf8 -m scripts.RunT3GeneralAll16FromReference `
  --run-id thesis16_v1 --preflight `
  --reuse-b-run thesis_b_m1_128 --export-bc-betas
```

The preflight writes
`output/studies/general_scalar_extended/thesis_general_all16/thesis16_v1/preflight.json`.
The selected `--reuse-b-run` must refer to the successful B(-1) direct pilot.
`--export-bc-betas` marks C(-1) beta exports as generatable, and the subsequent
fit actually exports them using Wolfram/RGBeta.

Non-mixed models require two things locally:

1. Saved canonical full-general symbolic UV/EFT1 beta exports under
   `output/studies/general_scalar_extended/complete_scalar_symbolic/full16_general_v1/`.
2. A saved original-topology C5 kernel under one of the repository's original
   matching artifact locations (`output/full/hypercharge/.../data/final_weinberg_coefficient.json`
   or `output/hypercharge/.../data/c5_coefficient.txt`).

If a matching source is absent, the model is marked **Blocked**. The runner
never invents an equivalent C5 expression or silently uses an unrelated model.
It does not depend on `output/studies/small_lambda_1e-7/fitted` or
`output/studies/lambda_5e-9_refit/fitted`.

## 2. Run first small campaign A(-2), B(-1), C(-1)

```powershell
python -X utf8 -m scripts.RunT3GeneralAll16FromReference `
  --run-id thesis16_pilot --models A_-2,B_-1,C_-1 `
  --reuse-b-run thesis_b_m1_128 --export-bc-betas `
  --sobol-points 8 --max-nfev 20
```

## 3. Run all 16 models

```powershell
python -X utf8 -m scripts.RunT3GeneralAll16FromReference `
  --run-id thesis16_v1 --reuse-b-run thesis_b_m1_128 `
  --export-bc-betas --sobol-points 128 --max-nfev 60
```

If some models fail because an artifact is missing, supply it locally, then
resume without losing completed work:

```powershell
python -X utf8 -m scripts.RunT3GeneralAll16FromReference `
  --run-id thesis16_v1 --reuse-b-run thesis_b_m1_128 `
  --export-bc-betas --sobol-points 128 --max-nfev 60 --resume
```

`--models A_-2` runs only that model. The first all-16 attempt may require
substantial computation (128 Sobol evaluations plus a least-squares stage
per model); fitting is performed on the local computer, **not in this ZIP**.

## Inputs and fitting

- Common scales: `1e7 -> 1e5 -> 1e3 -> 1e2 GeV`.
- Fixed `lambdaT3`: tracked per-model values: `0.1` in A–D and `0.05` in E.
- Additional quartic boundary: `--quartic-fraction 0.1` times `|lambdaT3|`.
- Allowed holomorphic bilinear: `--mu2-fraction 0.01` times the smaller scalar mass square.
- When allowed, conjugate Yukawa boundary fixed to `--extra-yukawa-fraction 0.05` times its ordinary Yukawa matrix. They run with the canonical RGEs, but are **not** independently fitted or fully inserted in C5.
- Sobol samples 18 real `y1`, `y2` entries per model, then uses the existing staged sensitivity-ranked bounded least-squares 5 -> 8 -> 12 parameters.
- Five normal-ordering NuFIT oscillation targets. `chi2 <= 5.86` is an exploratory acceptance label, not full model exclusion/statistical confidence.
- `--sobol-points` must be a power of two. Seed values and scalar values are frozen for each run ID.

## Outputs

`output/studies/general_scalar_extended/thesis_general_all16/<run_id>/`:

- `summary.csv` / `summary.json` — rows for every selected model, including blocked models.
- `settings.json` — immutable resumption signature.
- `<model>/fitted_config.json` — full fitted reference-style config.
- `<model>/result.json` — best chi2, five observables, fitted Yukawas, provenance and matching limitations.
- `<model>/observables.csv` — prediction vs NuFIT target.
- `<model>/general_uv_overrides.json` — for nonmixed models, exact full-general UV coupling values used (including extra allowed interactions).
- `<model>/error.json` — clear traceback for blocked models; rerun with `--resume` when inputs become available.

The original interactive RGE and sensitivity-oval HTML tools are **not** automatically
run by this first 16-model fitter. They require trajectories/sensitivity scans in
their expected saved-data schemas. The purpose of this stage is to obtain
reproducible best-fit parameters and observable tables first.

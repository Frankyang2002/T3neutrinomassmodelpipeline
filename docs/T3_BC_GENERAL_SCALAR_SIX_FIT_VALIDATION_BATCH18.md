# Batch 18 — Six-fit B/C general-scalar validation and consolidation

## Purpose

Consolidate the **six accepted conditional exploratory** B(-1)/C(-1) neutrino fits from Batches 11–13 and 17, check their provenance and saved observables, then independently replay their existing **candidate RGBeta** running and conditional **one-scale Matchete C5**. The intent is a numerical-physics benchmark audit, **not** a new EFT matching derivation.

The six input directories are **read only**:

| Study | Model | Input run ID | Format |
|---|---|---|---|
| `small_uv_running_5e-9` | `B_-1` | `b_uv_constrained_pilot_v1` | Batch 17 `fit.json` |
| `small_uv_running_5e-9` | `C_-1` | `c_uv_constrained_pilot_v1` | Batch 17 `fit.json` |
| `small_uv_frozen_5e-9` | `B_-1` | `b_small_fit_v2` | Batch 12 continuation |
| `small_uv_frozen_5e-9` | `C_-1` | `c_small_fit_v1` | Batch 11 fit |
| `large_uv_frozen_0p5` | `B_-1` | `b_large_fit_v1` | Batch 11 fit |
| `large_uv_frozen_0p5` | `C_-1` | `c_large_fit_v2` | Batch 12 continuation |

The UV-running solutions fit `lambdaT3` alongside the 18 real Yukawas; they are not fixed-`5e-9` repetitions of the original UV-running study. Frozen cases retain their original fixed scalar boundary.

## What is checked

**Saved-input checks** (both modes):

- Six exact identities, status/acceptance flags, `chi2 <= 5.86`.
- Fitted Yukawas and UV `lambdaT3` agree between the saved fit and its complete fitted configuration.
- Source input SHA-256 fingerprints match the files currently on disk. Any mismatch stops the run.
- Saved predictions reproduce the five-variable NuFIT chi-square with the currently saved `configs/nufit_6_1_ic24_no.json` target.
- Saved complex symmetric C5 matrices at 1 TeV and 100 GeV are finite and nonzero.
- Scalar 2x2 quadratic mass eigenvalues are positive; UV-running cases satisfy `rhoF <= 0.5` and 0.5–5 TeV at 1 TeV.
- Frozen-study `rhoF` estimates are formed from the **100 TeV scalar starting masses** and 100 TeV heavy masses, never mixing values evaluated at different scales.

**Default full replay** additionally:

- Recreate the fitted UV/EFT1 boundary using the *saved best-fit Yukawas*, not fresh seeds.
- UV-running B/C: candidate RGBeta 10^7 -> 10^5 GeV; frozen B/C: no UV RGE; all six: candidate EFT1 10^5 -> 10^3 GeV.
- Independently evaluate the same conditional saved one-scale C5 and run SM+Weinberg 10^3 -> 10^2 GeV, cross-checking the saved five observables and chi-square.
- Sample running scalar 2x2 eigenmasses, heavy-fermion ratios, and maximum raw scalar quartic-coordinate magnitude. For UV, also sample largest singular value of the running `y1` and `y2` matrices.
- Flag coordinate-level perturbativity when `max |lambda| > 4*pi` or `max singular(y1,y2) > sqrt(4*pi)`. **These are numerical screening conventions**, not a complete physical perturbative-unitarity or potential-boundedness analysis.

No production RGE, Matchete coefficient, fitted input, or historical result is edited. The runner fails **before writing output** when any saved-file check or replay comparison fails.

## Project paths (new files)

- `scripts/ValidateT3BCGeneralScalarSixFits.py`
- `tests/python/test_t3_bc_general_scalar_six_validation.py`
- `docs/T3_BC_GENERAL_SCALAR_SIX_FIT_VALIDATION_BATCH18.md`

## PowerShell commands

From `C:\thesistemp\T3neutrinomassmodelpipeline`, with `.venv` active:

```powershell
Expand-Archive -LiteralPath "C:\path\to\T3_BC_General_Scalar_Six_Fit_Validation_Batch18.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_bc_general_scalar_six_validation.py

# Full RG and C5 replay, 10 saved samples per running interval, new-only output:
python -m scripts.ValidateT3BCGeneralScalarSixFits --run-id bc_six_validation_v1

Get-Content -Encoding UTF8 output\studies\general_scalar_extended\validation\bc_six_validation_v1\summary.md
```

For a fast *saved-file-only* source/fit check (not full trajectory validation), with a different output ID:

```powershell
python -m scripts.ValidateT3BCGeneralScalarSixFits --run-id bc_six_saved_v1 --saved-only
```

The outputs for the full run are:

- `output/studies/general_scalar_extended/validation/bc_six_validation_v1/report.json`
- `output/studies/general_scalar_extended/validation/bc_six_validation_v1/summary.md`
- `output/studies/general_scalar_extended/validation/bc_six_validation_v1/trajectory_samples.csv`

Use a **new** `--run-id` if an output directory already exists. The report never claims a full 16-model general-scalar fit. No new numerical fit is performed.

## Physical limits

The candidate signed Matchete/RGBeta quartic identification is **not** certified; the saved C5 is a one-scale insertion, not complete sequential 100 TeV fermion and 1 TeV scalar Wilson-operator matching. Scalar masses from 2x2 RG quadratic matrices are not pole masses. Sampling positive masses does not guarantee vacuum stability for the complete general potential. The 100 TeV/1 TeV common-scale approximation persists, particularly when a running scalar mass is several TeV. Extremely small `chi2` values are existence checks with many adjustable Yukawa parameters, not evidence of unique predictions.

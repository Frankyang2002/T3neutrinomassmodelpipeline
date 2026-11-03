# Continuing the poor/blocked full-general T3 fits

## Scope and inputs

This is a **targeted continuation of independent 16-model fits**, not the fixed-A(-2)-Yukawa comparison. The untouched source is

`output/studies/general_scalar_extended/thesis_general_all16/thesis16_general_c5_full_v1/`.

The checked campaign snapshot has 11 exploratory accepted fits, four exploratory poor fits, and one blocked case:

| Model | Previous conditional chi2 | Previous status |
|---|---:|---|
| A(-2) | 34.5294821 | Poor fit |
| B(+1) | 54.3916572 | Poor fit |
| C(-1) | 35.1605246 | Poor fit, mixed-scalar conditional matching |
| D(+0) | — | Blocked: negative real self-bilinear at mass-basis matching |
| D(+2) | 14.7397265 | Poor fit |

The source's local refinement stopped after **12 active entries**. This continuation starts from each saved best point and refines **12, then all 18 real entries** of `y1`, `y2` using the existing `refine_sobol_candidates` optimiser. It **does not refit** the other 11 models and does not run a fresh wide Sobol search. Fixed general scalar parameters, initial scales, and matching conventions are inherited from the saved source campaign. C(-1) uses the saved B/C general scalar + mixed-bilinear evaluator; the others use canonical full-general renormalisable UV and EFT1 beta exports and general Matchete one-scale C5.

The screening criterion is `chi2 <= 5.86`. This is the saved campaign's exploratory acceptance screen, not a calibrated global exclusion threshold. These calculations are **conditional one-scale C5 predictions**, not certified complete sequential EFT matching (heavy fermions -> scalar EFT1 -> SMEFT). In particular, coupling dictionary normalization and full Wilson running remain uncertified.

## Install and run (PowerShell)

Unzip at the local repository root (do not extract into `scripts` directly):

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\T3_general_fit_continuation.zip" -DestinationPath . -Force
python -X utf8 -m scripts.ContinueT3GeneralPoorFits --self-test
python -X utf8 -m scripts.ContinueT3GeneralPoorFits --run-id continue_general16_v1 --preflight
python -u -X utf8 -m scripts.ContinueT3GeneralPoorFits --run-id continue_general16_v1
```

The last command runs local RGE integrations, C5 matching and optimisation. It **does not touch GitHub or overwrite** `thesis16_general_c5_full_v1`.

If interrupted after some models finish:

```powershell
python -u -X utf8 -m scripts.ContinueT3GeneralPoorFits --run-id continue_general16_v1 --resume
```

For a deeper *second* continuation, use a **different run ID** and adjust budgets, e.g. `--restarts 2 --stage18-nfev 300`. The script continues from the original saved study each time; it does not yet chain from a prior continuation. Existing accepted models remain untouched.

## D(+0) mass-basis blocker

The full-general self-bilinear Matchete C5 kernel currently requires a **positive real** bilinear at the scalar matching scale. The original D(+0) attempt produced no valid starting point, so there is no prior fitted minimum to continue. The script first tries the model's saved general reference texture, then the fitted D(-2) general texture under D(+0)'s own RGEs.

If D(+0) remains blocked, that is a **matching/benchmark issue**, not evidence that the representation cannot fit oscillation data. One optional diagnostic changes D(+0)'s scalar-bilinear UV boundary while keeping the rest of the settings:

```powershell
python -u -X utf8 -m scripts.ContinueT3GeneralPoorFits `
  --run-id continue_general16_d0_mu2_05_v1 --models D_+0 `
  --d0-mu2-fraction 0.05
```

**This changes a non-Yukawa input.** Results from that diagnostic are marked `modified_scalar_boundary=true` and are **excluded from the default merged all-16 table**. Do not label this an improvement at the original fixed benchmark.

## Outputs

`output/studies/general_scalar_extended/thesis_general_all16/<run-id>/`

- `summary.csv`, `summary.json`: the five targeted models, before/after chi2 and blockers.
- `merged_all16.csv`, `merged_all16.json`: original 16-model table, replacing rows only when the continuation improved chi2 **without modifying fixed scalar boundaries**.
- `<model-tag>/result.json`: best parameters, chi2, observables, optimiser stages, C5 and provenance.
- `<model-tag>/fitted_config.json`: complete new fitted model configuration (if evaluated).
- `<model-tag>/general_uv_overrides.json`: the continued full-general UV boundary dictionary (nonmixed models only).
- `<model-tag>/c5_conditional.json`: full symmetric C5 matrices at the scalar threshold and low scale.
- `<model-tag>/error.json`: explicit blocker where the calculation could not run.

The previous `A_-2` benchmark and the separate **shared A(-2) all-16 comparison** remain unchanged until a successful new A(-2) point is selected as that study's common Yukawa reference. Never silently replace the reference while comparing representations.

## Source verification

Read-only GitHub snapshot: `main` commit `d1eb7bc8c9fd4b6473752d09ba2106f011a845ae` (3 November 2026). Actual local files take precedence if they differ. The runner reuses the current repository's `scripts/RunT3GeneralAll16FromReference.py`, `scripts/RunT3BCGeneralFromReference.py`, `scripts/RunT3GeneralScalarBCExploratoryFit.py`, and `Numerical/fitting/SobolBenchmarkSearch.py`.

# Shared fitted A(-2) Yukawas across 16 general-T3 models

## Scope

Additive, **local-only** runner `scripts/RunT3SharedA2General16Reports.py`.

- Reads the completed A(-2) fit from `output/studies/general_scalar_extended/thesis_general_all16/thesis16_general_c5_full_v1/A_m2/` and validates the saved `result.json`, `fitted_config.json`, `settings.json` and matching-file SHA256 using the installed `scripts.RunT3Aminus2YukawaBenchmark.load_fit`.
- Copies the **same real 3x3 fitted UV `Y1` and `Y2`** into each model; the other UV parameters remain **model-specific**. No refitting.
- Uses `RunT3GeneralAll16FromReference.NonMixedEvaluator` plus model-specific canonical full-general RGBeta and Matchete one-scale C5 for the **14 nonmixed** models. The six A/D/E self-bilinear models retain their split scalar masses through the matched kernel.
- Uses the original `RunT3GeneralScalarBCExploratoryFit.CandidateEvaluator` for **B(-1)/C(-1)**, whose mixed-scalar matching is **conditional**, not independently certified complete-general.
- Additional symmetry-allowed conjugate Yukawas are included **only where present in the model's UV beta export**. By default, their boundary matrix is `extra_yukawa_fraction` times the corresponding **A(-2) reference** matrix, exactly reproducing the A(-2) donor campaign's `y2Conj` boundary value. `--conj-source model-reference` instead reproduces each model's own original reference-matrix convention. This avoids falsely claiming that all 16 independently fit additional Yukawas.
- Re-evaluates the RGEs and `C5` at every gamma value using the existing numerical implementation. Three standard `BenchmarkSensitivity` directions:
  - T1: λT3 → γλT3; Y1,Y2 → γ^-1/2 Y1,Y2
  - T2: λT3 → γλT3; Y2 → γ^-1 Y2
  - T3: Y1 → γY1; Y2 → γ^-1 Y2
- When λT3 changes, other **new quartics seeded as fractions of λT3** also change proportionally, as in the existing `fixed_overrides`/B-C constructors. Conjugate Yukawa boundaries remain fixed through each gamma trajectory. T1/T2 are therefore parameter sensitivities of the implemented general benchmark, not exact scaling symmetries.
- Saves the real numerical observable predictions, chi2, C5/UV/EFT threshold details for each gamma=1 centre; failed gamma values are explicitly tagged.
- Feeds saved diagnostics to **existing** `Numerical.plotting.BuildFittedModelEnvelopes.build` for an interactive `sensitivity_oval.html` with model selectors, observable selectors, gamma bounds and descriptive covariance ellipses.
- Saves `index.html`, `comparison.csv`, `gamma_points.csv` and manifest. **No RGE trajectory plots are fabricated**: the current general fitter returns true UV/EFT1 *endpoints*, not scale-resolved samples compatible with `InteractiveModelComparison`.

## Installation

Extract the ZIP from the project root using PowerShell:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive "$HOME\Downloads\T3_Shared_A2_General16_Interactive_Reports.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_thesis_shared_a2_general16_reports.py
```

## Full 16-model interactive sensitivity study

```powershell
python -X utf8 -m scripts.RunT3SharedA2General16Reports `
  --donor-run-id thesis16_general_c5_full_v1 `
  --run-id shared_A2_general16_ovals_v1 `
  --gamma-min 0.7 `
  --gamma-max 1.4 `
  --gamma-points 7
```

The original A(-2) full16 result must exist and the 14 model-specific Matchete C5 expressions and general beta exports must exist. To troubleshoot a short three-model scan first, use a different run ID and:

```powershell
python -X utf8 -m scripts.RunT3SharedA2General16Reports `
  --donor-run-id thesis16_general_c5_full_v1 `
  --run-id shared_A2_three_smoke_v1 `
  --models A_-2,B_-1,C_-1 `
  --gamma-points 3
```

If interrupted, rerun the exact command with `--resume`; successful completed model scans are reused. Model-level failures are written to `A_m2/error.json` etc. Run returns `2` if any selected model is blocked or incomplete; no partial result is falsely marked complete.

## Result paths

```
output/studies/general_scalar_extended/shared_a2_general16_reports/shared_A2_general16_ovals_v1/
  index.html                    cross-model HTML comparison
  sensitivity_oval.html         interactive T1/T2/T3 prediction-space ovals
  comparison.csv                gamma=1 observables and chi2
  gamma_points.csv              all actual computed gamma points
  comparison_summary.json       interface for existing oval generator
  manifest.json                 provenance and limitations
  settings.json                 resumable settings and donor fingerprints
  A_m2/centre.json              physical endpoints and C5 for shared benchmark
  A_m2/sensitivity_diagnostics.json
  ... per-model data for each of the 16 models ...
```

View the report using:

```powershell
Start-Process "output\studies\general_scalar_extended\shared_a2_general16_reports\shared_A2_general16_ovals_v1\index.html"
```

### Verification and limitations

The A(-2) gamma=1 evaluation is automatically checked against the numerical chi2 originally saved in the donor 16-model run, so an accidental change of its extra Yukawa boundary is reported.

This is a **conditional exploratory comparison**, not an independently validated full sequential EFT result. The descriptive ovals are sample covariance ellipses, **not** physical uncertainty/confidence regions; shared Yukawa inputs are deliberately *not* separate fits for each model. General Matchete coupling-map normalisations and the complex phases of self-bilinears remain unverified.

This package creates no plots without the user's local numerical exports and Matchete coefficient files. All offline tests test the program adapter, not a live Wolfram integration.

# Shared A(-2) general 16-model: 21-point sensitivity plus real UV/EFT1 RGEs

**Full replacement project files:**

- `scripts/RunT3SharedA2General16Reports.py` — existing shared-A2 scan and oval orchestration, enhanced with RGE generation, `--rge-only`, and HTML links.
- `scripts/RunT3SharedA2General16RGE.py` — added module that uses **the existing** `Numerical.running.UVRunner.run_uv_segment` and `Numerical.running.IntermediateScalarRunner.run_intermediate_scalar_segment`, with `save_scales_gev`, for genuine trajectory samples at the shared A(-2) gamma=1 point. It also renders an offline interactive UV/EFT1 HTML dashboard. It does not implement separate beta functions or matching.
- `tests/python/test_thesis_shared_a2_general16_reports.py` — full replacement original shared-A2 tests.
- `tests/python/test_thesis_shared_a2_rge_html.py` — extra tests for RGE samples, endpoint parity, HTML coverage and the new gamma grid.

All are project-relative paths. The ZIP does not contain or overwrite the general fitting or matching implementation. It does not modify GitHub.

## PowerShell install

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive "$HOME\Downloads\T3_Shared_A2_Gamma21_Interactive_RGE.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_thesis_shared_a2_general16_reports.py tests/python/test_thesis_shared_a2_rge_html.py
```

## New 21-point sensitivity scan: gamma=0.8..1.2

```powershell
python -X utf8 -m scripts.RunT3SharedA2General16Reports `
  --donor-run-id thesis16_general_c5_full_v1 `
  --run-id shared_A2_general16_gamma08_12_v1 `
  --gamma-min 0.8 --gamma-max 1.2 --gamma-points 21 `
  --rge-points 41
```

The established oval renderer requires the exact gamma=1 centre. A 20-point grid evenly spanning 0.8..1.2 does not guarantee gamma=1, so this uses **21** values, with ten below and ten above 1.0. Each model can have 3*21=63 sensitivity evaluations, requiring up to 1,008 model evaluations across all 16 models. Saving RGE curves requires one additional UV/EFT1 integration for each of the 14 nonmixed models; no additional Matchete coefficient evaluation is required for these dense curves.

Resume an interrupted run with the identical command plus `--resume`.

## Add RGE HTML to the existing seven-point run without repeating the scans

```powershell
python -X utf8 -m scripts.RunT3SharedA2General16Reports `
  --donor-run-id thesis16_general_c5_full_v1 `
  --run-id shared_A2_general16_ovals_v1 `
  --rge-only --rge-points 41
```

This option reads the saved donor and gamma=1 centre files, re-integrates only the UV and EFT1 ODEs with 41 saved scales per stage, verifies endpoint equality with the saved gamma=1 result, writes `rge_trajectories.json` per available model, and regenerates report HTML. It does not rescan gamma or redo Matchete matching.

To rebuild just HTML from saved scans/trajectories, use `--render-only` with the saved run ID.

## Report outputs

`output/studies/general_scalar_extended/shared_a2_general16_reports/<run-id>/`

- `index.html`: existing model overview with links to both HTML views;
- `sensitivity_oval.html`: existing interactive T1/T2/T3 oval visualisation, including non-mixed models and conditional B/C where the numerical scans succeed;
- `interactive_RGE.html`: actual sampled scale-dependent running from gamma=1 across 14 nonmixed models, with model and coupling selectors, UV/EFT1 toggles and logarithmic-value option;
- `rge_status.json`: explicit per-model RGE coverage;
- `<model_tag>/rge_trajectories.json`: actual numerical ODE points with solver information and source fingerprint (for nonmixed models);
- `<model_tag>/rge_error.json`: reason for unavailable trajectories;
- `gamma_points.csv`, `comparison.csv`: original numerical scan outputs.

**Scope:** RGEs here are only the renormalisable UV and EFT1 beta systems. The old `Numerical.plotting.InteractiveModelComparison` expects an additional full SMEFT/Weinberg scale-by-scale diagnostic that the current general fitter does not produce. This package does not invent such a trajectory or relabel the older restricted-matching numerical outputs. The existing `Numerical.plotting.BuildFittedModelEnvelopes` is reused for sensitivity ovals.

**Mixed B(-1)/C(-1):** Their conditional mixed-scalar numerical implementation is used for gamma sensitivity, but it does not expose comparable dense UV/EFT1 samples. They are explicitly omitted from the RGE curves (and named in `rge_status.json`), not replaced by 2-point lines. They remain shown in the oval dashboard if valid.

**Physics status:** all results conditional/provisional one-scale matching; no sequential general EFT Wilson matching certification. All additional symmetry-allowed couplings continue to be managed by each model's existing general evaluator. This package does not change any Lagrangian, Yukawa or matching expression, or the previous signed-real bilinear fix.

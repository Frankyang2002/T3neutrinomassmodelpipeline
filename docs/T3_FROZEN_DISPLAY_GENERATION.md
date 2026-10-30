# Frozen general-scalar display generation

This package adds the **data-generation step** missing from the previous frozen-model display work.

## Implemented

### 1. `scripts/GenerateT3FrozenRGEData.py`
Builds saved JSON trajectories for the frozen report.

What is actually implemented:
- reads accepted frozen fit states from
  `output/studies/general_scalar_extended/bulk_48/<run-id>/<case>/<model>/state.json`
- writes trajectory JSON files under
  `output/studies/general_scalar_extended/frozen_display/<run-id>/trajectories_v1/`
- stores:
  - flat UV Yukawa curves: `||Y1||_F`, `||Y2||_F`, `max |Y|`
  - a saved low-energy endpoint proxy over the EFT interval:
    `sqrt(|Δm²3ℓ|)` plus accepted low-energy observables

What is verified:
- frozen UV inputs are not RG-run in the accepted frozen fits;
- the accepted states contain Yukawa parameters and low-energy oscillation observables.

What is inferred:
- the displayed EFT segment is an endpoint proxy, not a reconstructed continuous Weinberg-operator evolution.

### 2. `scripts/GenerateT3FrozenSensitivityScans.py`
Builds the genuine saved T1/T2/T3 scan CSV used by the oval report.

Implemented directions:
- `T1`: scale first row of `Y1` by `γ`
- `T2`: scale second row of `Y2` by `γ`
- `T3`: scale third row of both Yukawas by `γ`

This uses the existing frozen evaluator from `scripts.RunT3GeneralScalar48Campaign` and writes one CSV row per `(model, case, direction, γ)`.

### 3. `Numerical/plotting/BuildFrozenGeneralScalarLambdaReports.py`
Builds:
- `index.html`
- `interactive_rge.html`
- `sensitivity_report.html`
- `all_results.html`
- `all_results.csv`

The pages are self-contained normal HTML/JS reports and do not require external libraries.

## PowerShell commands

Project root:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

python -X utf8 -m pytest -q tests\test_frozen_display_pipeline.py

python -X utf8 -m scripts.BuildT3FrozenDisplaySummary `
    --run-ids frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1

python -u -X utf8 -m scripts.GenerateT3FrozenRGEData `
    --run-id rge_v1 `
    --run-ids frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1

python -u -X utf8 -m scripts.GenerateT3FrozenSensitivityScans `
    --run-id sensitivity_v1 `
    --run-ids frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1 `
    --source-campaign overnight_six_v1_bulk48

python -X utf8 -m Numerical.plotting.BuildFrozenGeneralScalarLambdaReports `
    --run-ids frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1 `
    --trajectory-dir .\output\studies\general_scalar_extended\frozen_display\rge_v1\trajectories_v1 `
    --sensitivity-csv .\output\studies\general_scalar_extended\frozen_display\sensitivity_v1\frozen_gamma_points.csv `
    --output-dir frozen_massnorm_reports_v1

Start-Process .\Reports\output\display_frozen_general_scalar\frozen_massnorm_reports_v1\index.html
```

## Intended project paths

- `Numerical/plotting/BuildFrozenGeneralScalarLambdaReports.py`
- `scripts/GenerateT3FrozenRGEData.py`
- `scripts/GenerateT3FrozenSensitivityScans.py`
- `scripts/BuildT3FrozenDisplaySummary.py`
- `tests/test_frozen_display_pipeline.py`
- `docs/T3_FROZEN_DISPLAY_GENERATION.md`

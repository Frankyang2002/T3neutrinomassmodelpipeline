# Frozen general-scalar lambda-style reports

This package restores the **old Lambda-style HTML reporting layout** for the **accepted frozen-UV T3 general-scalar fits**.

## What is implemented

- `Numerical/plotting/BuildFrozenGeneralScalarLambdaReports.py`
  - reads accepted `state.json` files from frozen search runs;
  - builds three HTML pages and one index page;
  - reuses the old Lambda report templates where those files exist locally;
  - never fabricates sensitivity ovals or RGE curves.
- `scripts/BuildT3FrozenDisplaySummary.py`
  - collects the accepted frozen fits into `summary.csv`/`summary.json`.

## What is verified from the available sources

Verified in the previous work products:

- accepted frozen fits are saved as `output/studies/general_scalar_extended/bulk_48/<run-id>/<case>/<model>/state.json`;
- the old HTML template sources are `Numerical/plotting/InteractiveModelComparison.py` and `Numerical/plotting/BuildFittedModelEnvelopes.py`.

## What is inferred

The exact frozen display-data filenames for saved trajectories and saved T1/T2/T3 gamma scans were not verified in the available material.
So this builder accepts them as explicit inputs:

- `--trajectory-dir <directory containing JSON trajectories>`
- `--sensitivity-csv <csv containing genuine frozen T1/T2/T3 scan points>`

Required sensitivity CSV columns:

`model,case,direction,gamma,mass_hi_ev,mass_mid_ev,mass_lo_ev,delta_m21_sq_ev2,delta_m3l_sq_ev2,sin2_theta12,sin2_theta13,sin2_theta23`

Each accepted `(model, case)` should contain all three directions `T1`, `T2`, `T3`, at least three gamma values per direction, and include `gamma = 1`.

## PowerShell

From the repository root:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

python -X utf8 -m pytest -q tests\test_frozen_general_scalar_lambda_reports.py

python -X utf8 -m scripts.BuildT3FrozenDisplaySummary `
    --run-ids frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1

python -X utf8 -m Numerical.plotting.BuildFrozenGeneralScalarLambdaReports `
    --run-ids frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1 `
    --output-dir frozen_massnorm_reports_v1
```

If genuine saved data already exist, include them:

```powershell
python -X utf8 -m Numerical.plotting.BuildFrozenGeneralScalarLambdaReports `
    --run-ids frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1 `
    --trajectory-dir .\output\studies\general_scalar_extended\frozen_display\trajectories_v1 `
    --sensitivity-csv .\output\studies\general_scalar_extended\frozen_display\frozen_gamma_points_v1.csv `
    --output-dir frozen_massnorm_reports_v1
```

Open:

```powershell
Start-Process .\Reports\output\display_frozen_general_scalar\frozen_massnorm_reports_v1\index.html
```

## Scope

This package only changes **reporting/summary generation**. It does **not** rerun the neutrino fits, regenerate RGEs, or synthesize historical sensitivity envelopes.

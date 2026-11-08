# Rebuilt frozen-64 optimised_low sensitivity from accepted refits

This is a new report, not a modification of the old 64-case dashboard. It reads all 16 `ConditionalFitAccepted` Yukawa fit JSON files from `output/studies/general_scalar_extended/full_scalar_frozen_refits/<run-id>/`, verifies each recorded scalar-background hash and target hash, and recomputes model-specific conditional C5 + SM-Weinberg predictions for the fitted Yukawas and gamma-scaled Yukawas.

## Run (PowerShell in repository root)

```powershell
python -X utf8 -m pytest -q tests/test_rebuild_t3_frozen_full_scalar_refit_display.py
python -m scripts.RebuildT3FrozenFullScalarRefitDisplay --refit-run-id frozen_full_scalar_refit_low_pathfix_v2 --plan
python -u -m scripts.RebuildT3FrozenFullScalarRefitDisplay --refit-run-id frozen_full_scalar_refit_low_pathfix_v2 --display-id frozen_low_refitted_sensitivity_v1
Start-Process .\Reports\output\display_full_scalar_refitted\frozen_low_refitted_sensitivity_v1\index.html
```

Dependencies: numpy, scipy, matplotlib, plus the repository runtime and saved source data. Each model gets a multi-page sensitivity PDF. The HTML index links to the PDFs, CSV scan points, JSON fit diagnostics. The gamma directions are exactly those in `GenerateT3FrozenSensitivityScans.py`: T1 scales electron-flavour Y1 row; T2 scales muon-flavour Y2 row; T3 scales tau-flavour rows of both Yukawas.

The central gamma=1 evaluation must reproduce the saved refit chi2. Any model that fails this check is blocked, not silently included. New display data are **conditional fixed-background** results: scalar Yukawa-feedback is *not* rerun, A/D/E self-bilinear corrections to C5 remain absent and B/C(-1) use a conditional direct matching. Do not claim these are complete physical general-potential predictions. The sensitivity curves are not experimental confidence contours. No original outputs are overwritten.

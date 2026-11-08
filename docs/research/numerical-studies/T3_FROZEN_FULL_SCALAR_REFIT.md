# Frozen full-scalar fixed-background refit: source-path and diagnostics correction

## Verified bug
The previous default source was `.../frozen_general64_v1` but the trajectory files are nested under `.../frozen_general64_v1/optimised_low/<model-tag>/continuous_trajectory.json`. The new default includes the `optimised_low` segment.

The console now prints the exception that causes any `Blocked` result instead of printing `None`. Older `Blocked` files are not repaired automatically. Use a **new run ID**; do not use `--resume` on the previous failed campaigns.

## Important scope
This is an **approximate fixed-background Yukawa refit**, NOT a consistent full coupled scalar-Yukawa RGE refit: EFT1 scalar running and the scalar endpoint are read from the original frozen 64 calculation and are not recomputed when Yukawas change. Six self-bilinear C5 corrections are not included; mixed B/C mass-eigenstate one-scale matching and general coupling normalisation remain conditional. Fit quality cannot certify the full extended model.

## Commands (PowerShell, project root)
```powershell
python -m pytest -q tests\test_refit_t3_frozen_full_scalar_low.py
python -m scripts.RefitT3FrozenFullScalarLow --plan
python -u -m scripts.RefitT3FrozenFullScalarLow --run-id smoke_full_scalar_pathfix_v2 --models A_+0,B_-1 --max-nfev 5
python -u -m scripts.RefitT3FrozenFullScalarLow --run-id frozen_full_scalar_refit_low_pathfix_v2 --models all --max-nfev 30
```

New output: `output/studies/general_scalar_extended/full_scalar_frozen_refits/<run-id>/`. No original project files or historical outputs are changed.

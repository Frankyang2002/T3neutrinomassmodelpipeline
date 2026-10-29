# Frozen-UV persistent mass-normalized search

## Files
- `scripts/SearchT3FrozenPersistent.py`: full replacement for the previous diverse-restart runner.
- `tests/test_t3_frozen_persistent.py`: full replacement offline tests.

## Algorithm
For each Sobol candidate satisfying physical screening:

1. Read `delta_m3l_sq_ev2` from the authoritative evaluation.
2. Compute the common Yukawa scale factor `a=(abs(target)/abs(prediction))**0.25` only for finite, compatible mass-splitting signs.
3. Scale all 18 real Yukawa entries by `a`. Skip if outside the existing Sobol box.
4. Run the same authoritative frozen-UV evaluator on the scaled point (including constraints).
5. Require agreement with the `a**4` predicted atmospheric splitting to `--mass-scaling-rtol` (default 3%). If this fails, retain the unscaled point.
6. Use verified normalized points to select diverse local starts for standard 18D least-squares. Keep the best saved point and checkpoint every cycle.

This does not modify neutrino fitting targets, matching, scalar RGEs or screening. It is a proposal/seed transformation, **not** a substitute for the authoritative model evaluation. The quartic relation assumes a bilinear Yukawa dependence at fixed scalar input; it is checked for each transformed seed.

## Run (PowerShell)
Extract ZIP to `C:\thesistemp\T3neutrinomassmodelpipeline`, preserving folders.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -X utf8 -m pytest -q tests\test_t3_frozen_persistent.py
python -u -X utf8 -m scripts.SearchT3FrozenPersistent `
    --run-id frozen_massnorm_a0_v1 `
    --model A_+0 `
    --case small_uv_frozen_5e-9 `
    --source-campaign overnight_six_v1_bulk48 `
    --from-run frozen_persistent_diverse_a0_v1 `
    --from-run frozen_persistent_a0_v1 `
    --points 256 `
    --seeds 3 `
    --nfev 300 `
    --radii 0.02,0.05,0.10 `
    --min-seed-distance 0.08 `
    *> .\frozen_massnorm_a0_v1.log
```

To resume, use the same command with `--resume`, preserving all original options. A new run ID protects previous results. The checkpoint is:

`output/studies/general_scalar_extended/bulk_48/frozen_massnorm_a0_v1/small_uv_frozen_5e-9/A_p0/state.json`

`cycles.jsonl` adds `normalization_counts` with `normalized`, `unscalable`, `outside_search_box`, `scaling_mismatch`, etc. `normalized` means a candidate passed the quartic consistency check, **not** that its full neutrino chi-squared is acceptable. `CHECKPOINT` prints status counts. `Ctrl+C` stops the process; only completed cycles are guaranteed fully checkpointed.

## Interpretation
A mass-rescaled point need not have correct mixing angles or the solar-to-atmospheric mass-squared ratio. Least-squares still optimizes the original five whitened residuals; no residual is down-weighted. `chi2 <= 5.86` remains the acceptance condition.

# Frozen-UV 18-Yukawa radius sweep

The new `scripts/SweepT3Frozen18Radius.py` orchestrates the existing
`scripts/RefineT3Frozen18Stable.py` with 0.005, 0.01, then 0.02
local half-widths in unit-coordinate Yukawa space.

The `RefineT3Frozen18Stable.py` complete replacement is bundled to ensure the
runner has the exact expected source. The new sweep does not change model physics.
It uses each preceding stage's saved best fit as an additional input to the next
stage. It stops if chi2 <= 5.86 or after all three radii. Previous results are
preserved under new stage-specific run IDs.

Run from the repository root:

```powershell
python -X utf8 -m pytest -q tests\test_t3_frozen_stable18.py tests\test_t3_frozen_radius_sweep.py
python -u -X utf8 -m scripts.SweepT3Frozen18Radius `
  --run-id frozen_radius_sweep_a0_v1 `
  --model A_+0 `
  --case small_uv_frozen_5e-9 `
  --source-campaign overnight_six_v1_bulk48 `
  --from-run frozen_stable18_a0_v1 `
  --from-run frozen_full18_v1 `
  --seeds 3 --nfev 600 --diff-step 1e-6
```

To resume after a completed stage, repeat with `--resume`. An interrupted
incomplete inner stage is not automatically resumed; it is left intact for
inspection rather than silently overwritten.

Outputs: `output/studies/general_scalar_extended/bulk_48/frozen_radius_sweep_a0_v1/radius_sweep.json`,
and the corresponding `_r01`, `_r02`, `_r03` fit directories.

Exit code 0 means accepted, 2 means no fit at target after the sweep.
Check the reported bounds before assuming increasing radius will improve chi2.
No full sequential EFT matching is certified by this conditional fitter.

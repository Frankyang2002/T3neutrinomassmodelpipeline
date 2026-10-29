# Persistent diverse-restart frozen-UV search

Complete replacement for `scripts/SearchT3FrozenPersistent.py`, retaining the same
physics evaluator, acceptance test, checkpoint file names, and CLI arguments.

New behaviour:
- In each cycle, retain the lowest-chi2 candidate plus geometrically distinct
  starting points selected from the Sobol candidate pool. Separation is RMS
  distance in 18-dimensional, normalized Yukawa space (default `0.08`).
- Rotate local least-squares window radii `0.02,0.05,0.10` across local fits
  and cycles. These are normalized half-widths, clipped to the established
  physical search box.
- Save per-local `start_chi2`, `start_parameters`, `start_whitened_residual`,
  `final_chi2`, `final_parameters`, `final_whitened_residual`,
  `normalized_rms_displacement`, SciPy `nfev` and termination information
  into `cycles.jsonl`.
- The acceptance criterion remains `chi2 <= 5.86` and no future model is started.

## IMPORTANT: checkpoint migration

The provenance now includes the radii schedule and the minimum seed distance.
**Do not resume `frozen_persistent_a0_v1` using this replacement**, since it was
created by the old strategy and has a different provenance. Preserve that run
and use a NEW run-id. Import its saved best with `--from-run`.

From project root, with Python environment activated:

```powershell
python -X utf8 -m pytest -q tests\test_t3_frozen_persistent.py

python -u -X utf8 -m scripts.SearchT3FrozenPersistent `
  --run-id frozen_persistent_diverse_a0_v1 `
  --model A_+0 `
  --case small_uv_frozen_5e-9 `
  --source-campaign overnight_six_v1_bulk48 `
  --from-run frozen_persistent_a0_v1 `
  --from-run frozen_radius_sweep_a0_v1_r03 `
  --points 256 `
  --seeds 3 `
  --nfev 300 `
  --radii 0.02,0.05,0.10 `
  --min-seed-distance 0.08 `
  *> .\frozen_persistent_diverse_a0_v1.log
```

The imported best fit is re-evaluated before the new run begins. To resume
completed cycles of the **new** run, add `--resume` and use the identical
options (the `--from-run` inputs are then ignored). The run is unbounded unless
`--max-cycles` is specified. Stop at a checkpoint using Ctrl+C; the currently
unfinished cycle may be repeated on resume.

Output directory:
`output/studies/general_scalar_extended/bulk_48/frozen_persistent_diverse_a0_v1/small_uv_frozen_5e-9/A_p0/`

The production numerical physics calculation was not run during packaging.

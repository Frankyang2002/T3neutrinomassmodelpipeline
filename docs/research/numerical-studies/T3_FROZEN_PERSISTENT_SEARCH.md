# Persistent Sobol + full-18 refinement (frozen UV)

Project path: `scripts/SearchT3FrozenPersistent.py`. The existing `scripts/RefineT3Frozen18Stable.py` is required, unchanged.

## Operation

* Default target: chi-square <= 5.86 (mass/perturbativity screening unchanged).
* Each cycle samples a fresh scrambled, 18-dimensional Sobol sequence and performs up to 3 local full-18 trust-region refinements.
* `xtol` ends only the local attempt; the next cycle runs with a fresh Sobol seed.
* No preset number of cycles (`--max-cycles 0`). Interrupt with Ctrl-C; only **completed** cycles are checkpointed, so an interrupted cycle is repeated.
* Previous fits are imported and evaluated against current inputs; the existing files are never overwritten.
* Experimental conditional general-quartic approximation; not certified full sequential matching.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -X utf8 -m pytest -q tests\test_t3_frozen_persistent.py
python -u -X utf8 -m scripts.SearchT3FrozenPersistent `
  --run-id frozen_persistent_a0_v1 `
  --model A_+0 --case small_uv_frozen_5e-9 `
  --source-campaign overnight_six_v1_bulk48 `
  --from-run frozen_radius_sweep_a0_v1_r03 `
  --from-run frozen_full18_v1 `
  --points 256 --seeds 3 --nfev 300 --radius 0.02 `
  *> .\frozen_persistent_a0_v1.log
```

Resume with the same arguments and `--resume`. Do not alter settings when resuming; start a new run id instead. If you prefer a finite pilot, add `--max-cycles 4`, but use that same limit on resume or a new run id.

Outputs: `output/studies/general_scalar_extended/bulk_48/frozen_persistent_a0_v1/small_uv_frozen_5e-9/A_p0/state.json`, `cycles.jsonl`, and `fitted_config.json`.

**Important:** There is no guarantee of convergence. This script will run indefinitely if a suitable fit is not found; the user must stop it manually. It does not proceed to other models. Extreme hypercharge legacy C5 parity failures still need independent resolution.

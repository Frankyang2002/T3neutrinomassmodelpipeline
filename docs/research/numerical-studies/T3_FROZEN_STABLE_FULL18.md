# Frozen-UV local 18D fit: controlled numerical refinement

**Project path**: `scripts/RefineT3Frozen18Stable.py`

This runner re-evaluates saved candidates in the original candidate general-quartic model, then uses SciPy `least_squares` with all 18 Yukawa entries active. The fit is restricted to a narrow unit-coordinate neighborhood of each seed, with explicit 3-point numerical derivatives, Jacobian scaling, and a large budget. It saves evaluations, the best validated point, and SciPy termination diagnostics.

It **does not** alter scalar RGEs, mass screens, topology coupling, matching, or the physical model. It **does not** certify sequential EFT matching. It will not mark a poor fit accepted.

## PowerShell (project root)

```powershell
python -X utf8 -m pytest -q tests\test_t3_frozen_stable18.py
python -u -X utf8 -m scripts.RefineT3Frozen18Stable `
  --run-id frozen_stable18_a0_v1 `
  --model A_+0 --case small_uv_frozen_5e-9 `
  --source-campaign overnight_six_v1_bulk48 `
  --from-run frozen_escalating_v1 `
  --from-run frozen_full18_v1 `
  --seeds 3 --nfev 600 --radius 0.002 --diff-step 1e-6 `
  *> .\frozen_stable18_a0_v1.log
```

Result directory: `output/studies/general_scalar_extended/bulk_48/frozen_stable18_a0_v1/small_uv_frozen_5e-9/A_p0/`.

- `state.json`: best fit, complete SciPy termination diagnostics and input hashes; updated on improvements and after each seed.
- `fitted_config.json`: historical configuration with best new Yukawas, if one was evaluated.
- `evaluations.jsonl`: one line per numerical evaluation with status and chi-square.

Exit code `0` means chi-square <= 5.86; exit code `2` means a complete but unsuccessful fit. **Do not reuse the same run-id**: start a new run-id to preserve the outputs.

**Status**: offline syntax and structural tests only. No local physics replay has been run in this environment.

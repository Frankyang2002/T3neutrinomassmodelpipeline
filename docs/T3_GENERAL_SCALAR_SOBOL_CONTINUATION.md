# Bounded continued Sobol search — 16 × 3 conditional general-quartic campaign

This continues **unresolved** model/scenario cells from the existing `overnight_six_v1_bulk48` results, preserving all 20 previously accepted cells. It calls `scripts.RunT3GeneralScalar48Campaign.run_case` unchanged for each new random seed: 18 real Yukawa Sobol samples, existing mass/perturbativity screens, then the same local SciPy least-squares refinement. Accepted cells are not re-evaluated. Every round gets an immutable `fit.json` and (when available) `fitted_config.json` and is checkpointed individually. B(-1)/C(-1) historical accepted fits are preserved and never rewritten.

**Not an unbounded search:** `--max-rounds` caps rounds per unresolved cell, so computation terminates. A mass-screen rejection is not a poor oscillation fit and is not guaranteed fixable by more Yukawa points. All calculations are **conditional native general-quartic** studies, excluding extra holomorphic scalar quadratic matching and complete sequential EFT Wilson matching.

## Files to install (new files only)

- `scripts/ContinueT3GeneralScalarSobol.py`
- `tests/test_t3_general_scalar_sobol_continuation.py`
- `docs/T3_GENERAL_SCALAR_SOBOL_CONTINUATION.md`

Existing campaign, matching, and plot files stay unchanged.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -X utf8 -m pytest -q tests\test_t3_general_scalar_sobol_continuation.py

# Focus first on the A(-2) missing low-scalar donor solutions.
python -u -X utf8 -m scripts.ContinueT3GeneralScalarSobol `
    --source-campaign overnight_six_v1_bulk48 `
    --run-id sobol_recovery_a2_v1 `
    --models A_-2 --sobol-points 256 --max-nfev 150 --max-rounds 8 `
    *> .\sobol_recovery_a2_v1.log

# Or continue every unresolved cell using a new output ID:
python -u -X utf8 -m scripts.ContinueT3GeneralScalarSobol `
    --source-campaign overnight_six_v1_bulk48 `
    --run-id sobol_recovery_all_v1 `
    --sobol-points 128 --max-nfev 100 --max-rounds 8 `
    *> .\sobol_recovery_all_v1.log

# On interruption, append --resume to the SAME command.
```

Outputs: `output/studies/general_scalar_extended/bulk_48/<run-id>/campaign.csv`, `campaign.json`, `continuation_progress.json`, `summary.md`, and individual immutable `attempts/<case>/<model>/round_NNN/fit.json`. The run's `campaign.json` remains compatible with the existing `CompareT3GeneralScalarSharedA2` generator. A new shared comparison ID is required:

```powershell
python -u -X utf8 -m scripts.CompareT3GeneralScalarSharedA2 `
  --campaign sobol_recovery_all_v1 --run-id shared_after_recovery_v1 --allow-poor-donor
```

The shared report will still flag missing A(-2) donors and numerical failures. No automatic physics claim can be made about missing cells.

# Scalar-threshold EOMEFT Hermiticity diagnostic and conservative rewrite

## Live finding

After the canonical EOMEFT handoff, local Matchete passed canonical normalization,
heavy mass basis, gauge invariance and object definition in both B_-1 and C_-1.
Only `Hermiticity` was `False`. Consequently, no second-stage Match ran.

## Changes

1. Parse the saved EOMEFT after removing integrated F while X,Y are still registered Light.
   Record the **light-field** Hermiticity result in a compact JSON file.
2. Re-register X,Y Heavy. Require the original `CheckLagrangian` as before.
3. Only if Hermiticity is the sole failed check, attempt the documented identity-preserving
   Matchete simplifiers `HcSimplify` and `GreensSimplify` (with `RelabelIndices`).
4. A rewrite is accepted only if the full ordinary check becomes True, and detailed
   Hermiticity, canonical normalization, heavy mass basis, gauge invariance and
   object definition all become True. No additional HC terms are invented.
5. Export a compact `_hermiticity_diagnostics.json` even on failure. The CLI
   requires a successful, complete Hermiticity certificate for any claimed success.

If no equivalence-preserving simplification resolves the failure, the stage **still fails**
with a diagnostic JSON. In that case inspect `light_EFT1_hermiticity` and the rewrite attempts
before attempting further EFT1 integration. This patch is not a verified live repair.

## Full replacement files

- `Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl`
- `scripts/RunGeneralT3EFT1ScalarThresholdProbe.py`
- `tests/python/test_general_t3_eft1_scalar_threshold_probe.py`
- `docs/T3_GENERAL_EFT1_HERMITICITY_SECOND_THRESHOLD.md`

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_General_EFT1_Hermiticity_Recovery.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_eft1_scalar_threshold_probe.py
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model B_-1
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model C_-1
```

No need to rerun either fermion-first threshold; this batch consumes existing
canonical tree EOMEFT exports. No files in the production RGE path are changed.

# EFT1 scalar threshold: heavy-field free-Lagrangian rebase and Match audit

## Scope

Experimental **one-loop LLHH-sector probe** for `B_-1`, `C_-1`. This does not repair the complete scalar potential, evolve EFT1, subtract hard UV poles, compute physical sequential C5, or reproduce the oscillation fit.

## Observed failure before this update

Both restricted-sector (`--c5-relevant-only`) inputs passed `CheckLagrangian=True`, with 20/24 saved LLSS terms, 52 Higgs-containing quartics, and 21 background terms. Matchete's subsequent `Match` emitted diagnostic messages that `NewScalarX` and `NewScalarY` were defined **Heavy** but were **Light** in the supplied Lagrangian. Although `RunT3Matching` returned `Status=Success` through its processing stages, the final tree/loop difference did not pass the explicit `hbar`/heavy-field guards and correctly failed with `ScalarThresholdLoopCorrectionInvalid`.

Matchete stage success means the stage returned a non-`$Failed` value; it does **not** by itself establish that heavy fields were integrated out.

## New handling

1. Before changing the field registry, capture the actual `FreeLag[NewScalarX,NewScalarY]` under `Mass->{Light,MSX/MSY}`.
2. After re-registering both fields with `Mass->{Heavy,MSX/MSY}`, capture their new `FreeLag` under the heavy registry.
3. In the opt-in `--c5-relevant-only` matching subtheory **only**, substitute the heavy free sector for the light one: `Expand[scalarEFT-lightEigenFree+heavyEigenFree] // RelabelIndices`. Preserve all matched LLSS and Higgs-bearing interaction terms. Do **not** alter the saved EFT1 source file.
4. Require full `CheckLagrangian=True` for this restricted subtheory. Reject any genuinely invalid Lagrangian.
5. Export a JSON audit of the Matchete `RawEFT`, `EOMEFT` and `MatchedEFT` at **both** loop orders before validating the correction. Record expression head, character count, `hbar`, scalar-field presence, `GammaCC` and a short prefix.
6. Continue only if the difference is nonzero with explicit `hbar` and the final Matchete EFT contains no `NewScalarX`, `NewScalarY` or `NewFermion`. **There is no fallback that invents a contribution if Matchete matches zero.**

This is a power-counting-aware mass/free-Lagrangian handoff attempt. It is **not an algebraic-equivalence certification** between light and heavy field power-counting conventions. The live Matchete result is still unknown until the user runs the scripts.

## Paths (relative to repository root)

**Full replacements**
- `Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl`
- `scripts/RunGeneralT3EFT1ScalarThresholdProbe.py`

**New files**
- `tests/python/test_general_t3_eft1_heavy_free_rebase.py`
- `docs/T3_GENERAL_EFT1_SCALAR_MATCHING_HEAVY_FREE_REBASE.md`

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_General_EFT1_Scalar_Heavy_Free_Rebase.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_eft1_heavy_free_rebase.py tests\python\test_general_t3_eft1_c5_relevant_sector.py tests\python\test_general_t3_eft1_scalar_threshold_probe.py
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model B_-1 --c5-relevant-only
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model C_-1 --c5-relevant-only
```

Outputs remain under `output/general_t3_mass_eigen_matchete/eft1_scalar_threshold_probe/`, with independent `_c5_relevant` stems. New diagnostics per model:

- `<stem>_free_lagrangian_rebase.json`
- `<stem>_matching_stage_audit.json`

If `ScalarThresholdLoopCorrectionInvalid` remains, upload the **two matching-stage audit JSONs** and the **two free-Lagrangian rebase JSONs**. Those identify whether the issue is the free mass power counting or Matchete's handling of dimension-five LLSS insertions. Do not rerun earlier fermion-threshold one-loop matching.

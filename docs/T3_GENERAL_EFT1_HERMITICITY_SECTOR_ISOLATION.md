# T3 EFT1 Hermiticity: exact-sector isolation

## Verified live failure

In both `B_-1` and `C_-1`, the canonical first-threshold `EOMEFT` failed
`CheckLagrangian[...,DetailedOutput->True]["Hermiticity"]` while **X and Y
were still Light**. Reclassifying to Heavy did not cause the failure.
`HcSimplify` and `GreensSimplify` did not restore the check. Gauge invariance,
canonical normalization, heavy mass basis and object definitions passed in
the heavy-scalar check.

## This batch

The scalar-threshold probe first tests mathematically identity-preserving
rewrites of the EFT (`RelabelIndices`, `Expand` and `Expand`+`RelabelIndices`)
without changing coefficients. It accepts one only after ordinary and detailed
Matchete checks pass. Existing `HcSimplify` and `GreensSimplify` are retained,
with nonfatal simplifier messages no longer converted into `$Failed`.

If the full check still fails, the script **does not run scalar matching**.
Instead, it expands the *saved* EFT1 into disjoint operator sectors and
performs bounded detailed Matchete checks on each:

- Background (all terms other than LLSS and BSM scalar quartics)
- LLSS (all; XX, XY, YY individually)
- Heavy scalar quartics (all; purely heavy and Higgs portals individually)
- The full expanded EFT and EFT with each of the two operator classes excluded

These tests distinguish Hermiticity of each **complete operator class**;
checking a single chirality term would not be meaningful. They preserve
all original terms in the actual matching input and never add/remove HC
operators. The diagnostic output is not a certificate of pole cancellation,
physical `C5`, or a complete EFT1 RGE.

## Full replacement paths

- `Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl`
- `scripts/RunGeneralT3EFT1ScalarThresholdProbe.py`
- `tests/python/test_general_t3_eft1_scalar_threshold_probe.py`

New document:

- `docs/T3_GENERAL_EFT1_HERMITICITY_SECTOR_ISOLATION.md`

## Windows PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_General_EFT1_Hermiticity_Sector_Isolation.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_eft1_scalar_threshold_probe.py
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model B_-1
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model C_-1
```

On a genuine Hermiticity failure inspect:

`output/general_t3_mass_eigen_matchete/eft1_scalar_threshold_probe/B_m1_final_eft1_scalar_threshold_probe_hermiticity_sector_audit.json`

and its `C_m1` counterpart. In particular compare the `hermiticity` keys in
`sector_checks` for `background_no_LLSS_or_BSM_quartics`,
`background_plus_all_LLSS`, and `background_plus_all_BSM_quartics`.
The log contains a concise failure-sector summary. Upload those two small
JSON diagnostics together with the failed model's existing
`_hermiticity_diagnostics.json` for the subsequent exact physics fix.

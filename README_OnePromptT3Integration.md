# T3 — one-prompt local general-Yukawa RGE integration

**Scope:** workstream A, local repository only; no GitHub writes or changes to the physical matching pipeline.

## Files (install ZIP at repository root)

| ZIP path | Function |
|---|---|
| `scripts/IntegrateT3ExtraYukawas.py` | Guarded full replacement of the user's **current** `RGE/running/rgbeta/T3RGBetaModel.wl`, with automatic timestamped backup; opt-in `AddYukawa[y1Conj/y2Conj]` and `BetaTerm` export; preserves previously corrected `MF` doublet epsilon. |
| `scripts/RunT3ExtraYukawaStack.py` | Production **general-scalar** RGBeta UV runner for 5 smoke or 16 full models; reads existing `RGE/general/generated/*_extra_quartics.wl`; validates output beta keys and model metadata. |
| `Numerical/experimental/RunT3ExtraYukawaLeadingLog.py` | Frozen-background RGBeta-only one-step numerical test for the six new flavour matrices; computes 3×3 matrix beta shifts; keeps separate from the old fixed-layout coupled ODE state. |

No existing matching, neutrino-mass, UV RGE report, numerical fitter, or fitted reference is overwritten.

## Exact PowerShell commands

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$HOME\Downloads\T3OnePromptIntegrationStack_Local.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

# Pure Python structural tests before changing production source
python scripts\IntegrateT3ExtraYukawas.py --self-test
python scripts\RunT3ExtraYukawaStack.py --self-test
python Numerical\experimental\RunT3ExtraYukawaLeadingLog.py --self-test

# Guarded local full-file replacement and backup
python scripts\IntegrateT3ExtraYukawas.py --install
python scripts\IntegrateT3ExtraYukawas.py --verify

# Reverify the corrected charged doublet Dirac mass
python Lagrangian\RunT3A0DiracMassWard.py --all

# Five model smoke tests, then all sixteen general-scalar RGBeta UV models
python scripts\RunT3ExtraYukawaStack.py --smoke
python scripts\RunT3ExtraYukawaStack.py --all

# All six extra-Yukawa RGBeta-only numerical one-step benchmarks
python Numerical\experimental\RunT3ExtraYukawaLeadingLog.py --all
```

`--smoke` checks A(0), A(-2), B(-1), D(-2), E(0). `--all` tests 16 models with six optional independent Yukawa couplings. Each native Wolfram process sets `T3_EXPERIMENTAL_GENERAL_FERMIONS=1` **only for that process**. Existing restricted runs are not changed by default.

## Output paths

- `Reports/output/workstream_a/integrated_general_yukawa_uv/native/<model>.json` — 16 model-wise full general quartic UV beta exports.
- `Reports/output/workstream_a/integrated_general_yukawa_uv/logs/<model>.log` — native Wolfram execution diagnostics.
- `Reports/output/workstream_a/integrated_general_yukawa_uv/summary.json` — model-wise registry / full beta-key validation.
- `Reports/output/workstream_a/general_yukawa_leading_log/summary.json` — six test-point one-step matrix beta shifts.
- `RGE/running/rgbeta/T3RGBetaModel.before_general_yukawa_<timestamp>.wl` — the automatic backup of the *entire* original model source.

## Exact implementation and certification boundaries

- Source basis: `T3RGBetaModel.wl`, `RunGeneralT3RGBeta.wl`, original Matchete T3 interactions, and the isolated six-extra-Yukawa registration bridge in the public project main branch as inspected 2026-10-31.
- This package **requires** the user's previously verified A(0) `eps[SU2L[fund]]` charged-doublet mass fix to be present. It refuses to install if that correction is missing.
- It adds `y1Conj` for alpha=0 and `y2Conj` for alpha=-2 **only when** `T3_EXPERIMENTAL_GENERAL_FERMIONS=1`, matching the six independent Matchete Yukawa operators. No special self-conjugate scalar assumptions are added.
- It reuses the repository's existing general-scalar generator and full quartic projection, not the restricted 16-model quartic basis.
- The Python smoke tests and an **entirely simulated** 16-model Wolfram-export pipeline passed in the builder environment. Actual RGBeta and numerical evaluation against real exported beta formulas require the user's local run.
- **Not certified**: native RGBeta symbolic Ward checker, absolute Matchete–RGBeta triplet component axes and coupling phases, the physical `MF/y1/y2/yConj` action dictionary, full coupled ODE evolution with extended state, EFT matching, C5 neutrino-mass predictions and complete UV operator closure.
- The numerical runner deliberately seeds a test-only independent complex 3×3 Yukawa matrix in RGBeta coordinates and uses existing fitted-reference magnitudes only as **numerical placeholders**. **Do not interpret its shifts as physical neutrino-mass predictions.**

## Undoing the local source change

The installer prints the exact backup path. With all Wolfram programs closed, restore the same file, e.g.:

```powershell
Copy-Item -LiteralPath 'RGE\running\rgbeta\T3RGBetaModel.before_general_yukawa_YYYYMMDD_HHMMSS.wl' -Destination 'RGE\running\rgbeta\T3RGBetaModel.wl' -Force
```

Replace `YYYYMMDD_HHMMSS` with the timestamp printed by the installer. The three added Python files can remain installed; without the switch the production model continues to use its restricted Yukawa basis.

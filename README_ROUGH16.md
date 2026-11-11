# Batch 4 — Finish rough T3 16-model numerical campaign

This package **does not change solvers or physics parameters**. It consolidates the existing, completed 16-model coupled-IVP cases, selecting the successful Batch 3C A(-4), D(-2) cached non-HHSS calculations and the other 14 existing coupled-IVP calculations (including mixed B(-1), C(-1)). It uses **all 16 stage-freeze cases per model**.

Project files:

- `scripts/ConsolidateT3Rough16Results.py` — new full file; standard library only.
- `README_ROUGH16.md` — instructions.

Input campaign explicitly used:

- A(-4), D(-2): `output/studies/general_scalar_extended/experimental_coupled_cached_scalar_a4d2/scalar_cached_verified_source_20261111_235757/<model>/cases.json`
- Remaining 12 nonmixed: `output/studies/general_scalar_extended/experimental_simultaneous_eft1_general_nonmixed/simultaneous_nonmixed13_v1/<model>/cases.json`
- Mixed B(-1), C(-1): `output/studies/general_scalar_extended/experimental_simultaneous_eft1_mixed_BC/simultaneous_mixed_BC_v1/<model>/cases.json`
- Previous A(-4) reference is also checked for before/after comparison.

All input files are read-only. Output directory must not already exist.

## Run from PowerShell, in repository root

```powershell
$Project = "C:\thesistemp\T3neutrinomassmodelpipeline"
$Downloads = (New-Object -ComObject Shell.Application).Namespace('shell:Downloads').Self.Path
$Zip = Get-ChildItem -LiteralPath $Downloads -Filter "T3_Batch4_Rough16Completion*.zip" -File | Sort-Object LastWriteTime -Descending | Select-Object -First 1
if ($null -eq $Zip) { throw "Download the Batch 4 ZIP first" }
Expand-Archive -LiteralPath $Zip.FullName -DestinationPath $Project -Force
Set-Location -LiteralPath $Project
$Out = Join-Path $Project ("output\studies\general_scalar_extended\rough16_consolidated\rough16_" + (Get-Date -Format 'yyyyMMdd_HHmmss'))
python -X utf8 -m scripts.ConsolidateT3Rough16Results --project "$Project" --experimental-run scalar_cached_verified_source_20261111_235757 --output "$Out"
if ($LASTEXITCODE -ne 0) { throw "Consolidation blocked; see error above" }
Import-Csv -LiteralPath (Join-Path $Out 'models_1111.csv') | Format-Table model,scalar_context,chi2_1111,m1_eV,m2_eV,m3_eV -AutoSize
Get-Content -LiteralPath (Join-Path $Out 'summary.json') -Raw
Write-Host "Final results folder: $Out"
```

Output:

- `models_1111.csv` — 16-model all-stages-on observables, neutrino masses, chi².
- `all_256_stage_freeze_cases.csv` — all four-bit U/E/W/S cases.
- `single_stage_freeze_effects.csv` — 64 UV/EFT1/Wilson/SMEFT one-at-a-time freezes, relative to `1111`.
- `A4_D2_cached_vs_previous.csv` — 32 before/after case comparisons.
- `stage_freeze_delta_chi2.svg` — 16 × 16 heatmap of `chi2(case)-chi2(1111)` (colour saturates at ±1000).
- `summary.json` — file hashes, model-to-source mapping, complete caveats, 16 baseline results.

## Status of the physics

**Conditional rough numerical results only.** Normalisation *does* matter for physical mass scales and the resulting chi². Existing `spinor_factor=1` is a provisional convention. The two scalar approximations are mixed only at the report level, **never inside one model's EFT calculation**. The program never claims these are physical one-loop predictions. Explicit unresolved items include the signed Matchete/GammaCC→Weyl map, one-loop UV-hard minus EFT1-soft matching, scalar native-to-RGBeta quartic map, D(-2) quadratic mass direction, and complete joint scalar threshold.

No Wolfram run or new optimisation is needed: the numerical calculations already exist. If all inputs validate, the rough-results campaign is finished at this approximation level.

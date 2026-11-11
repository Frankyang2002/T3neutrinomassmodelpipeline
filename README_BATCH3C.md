# Batch 3C — Experimental full non-HHSS quartic in coupled EFT1 IVP

## Files (project-relative paths)

- `Numerical/running/T3CoupledQuarticContext.py` — new cached tensor loader, source-HHSS splice and context diagnostic collector.
- `scripts/RunT3CoupledScalarA4D2Experimental.py` — new full A(-4)/D(-2) 16-stage-freeze-case experimental runner (per model), reuses the existing coupled DOP853 UV/EFT1/SMEFT workflow.
- `scripts/TestT3CoupledScalarCached.py` — standalone tensor and coupled toy-RGE tests.
- `README_BATCH3C.md` — this file.

**No existing implementation files are replaced.** The new runner does not modify GitHub or historical results. It refuses to reuse any run ID and writes under:

`output/studies/general_scalar_extended/experimental_coupled_cached_scalar_a4d2/<run-id>/`

## Physics scope

1. The existing simultaneous `RunT3SimultaneousEFT1All16` runner and the source-resolved `HHSSLinearBasis` are retained. At each DOP853 RHS evaluation, the cached full real-scalar quartic tensor is evaluated with the instantaneous RGBeta scalar couplings.
2. The cached tensor's `HHSS` component is **replaced by the original source-resolved HHSS block** (all six index placements). Thus initial fermion matching, scalar-pair calibration, and finite two-scalar matching are not silently redefined.
3. The other quartic blocks (`H^4`, `S^4`, `S1^2 S2^2`) come from the earlier Batch 3B potential, including the D(-2) scalar-only directions. This is a **conditional replacement of the old non-HHSS ansatz**; physical quartic coordinate normalisation is still unverified.
4. The experimental runner preserves **UV -> integrate out all fermions -> EFT1 with both scalars -> integrate out both scalars at once -> SMEFT**. Neither scalar is integrated out at the fermion threshold.
5. A(-4) and D(-2) are supported; the mixed B(-1)/C(-1) are not. The D(-2) fixed-mass rotation is taken from the archived Batch3B tensor; the missing quadratic direction is unresolved. Do not regard the result as a certified running mass basis.
6. Existing scalar finite bubble and one-scale matching approximations remain unchanged, and no independent hard one-loop fermion matching is inserted.

## Prerequisite

A successful Batch3B run producing `report.json` and `quartic_response.npz` in **each** of `<Batch3B-root>/A_-4/` and `<Batch3B-root>/D_-2/`. The runner hashes/checks the cached NPZ and registry/seed, and reads the existing `load_snapshot()` reference state used to construct it.

## Usage (PowerShell)

```powershell
$ErrorActionPreference = "Stop"
$Project = "C:\thesistemp\T3neutrinomassmodelpipeline"
$Downloads = (New-Object -ComObject Shell.Application).Namespace('shell:Downloads').Self.Path
$Zip = Join-Path $Downloads "T3_Batch3C_ExperimentalCoupledScalar.zip"
if (-not (Test-Path -LiteralPath $Zip)) { throw "ZIP not found: $Zip" }
Expand-Archive -LiteralPath $Zip -DestinationPath $Project -Force
Set-Location -LiteralPath $Project

python -X utf8 scripts\TestT3CoupledScalarCached.py
if ($LASTEXITCODE -ne 0) { throw "Coupled scalar test failed" }
python -X utf8 -m scripts.RunT3CoupledScalarA4D2Experimental --self-test
if ($LASTEXITCODE -ne 0) { throw "Experimental runner self-test failed" }

$Previous = Get-ChildItem -LiteralPath (Join-Path $Project 'output\audits') -Directory |
    Where-Object { $_.Name -like 'batch3b_quartic_fixed_*' -or $_.Name -like 'batch3b_quartic_*' } |
    Where-Object { (Test-Path (Join-Path $_.FullName 'A_-4\report.json')) -and (Test-Path (Join-Path $_.FullName 'D_-2\report.json')) } |
    Sort-Object LastWriteTime -Descending | Select-Object -First 1
if ($null -eq $Previous) { throw 'Batch3B report directory missing. Run Batch3B first.' }
$ResponseRoot = $Previous.FullName
Write-Host "Using cached tensor response: $ResponseRoot"

$RunId = 'scalar_cached_A4_D2_' + (Get-Date -Format 'yyyyMMdd_HHmmss')
python -X utf8 -m scripts.RunT3CoupledScalarA4D2Experimental `
    --project "$Project" `
    --response-root "$ResponseRoot" `
    --source-run thesis16_merged_v1 `
    --run-id "$RunId" `
    --models A_-4,D_-2 `
    --check-only
if ($LASTEXITCODE -ne 0) { throw 'General-fit preflight failed' }

python -X utf8 -m scripts.RunT3CoupledScalarA4D2Experimental `
    --project "$Project" `
    --response-root "$ResponseRoot" `
    --source-run thesis16_merged_v1 `
    --run-id "$RunId" `
    --models A_-4,D_-2
if ($LASTEXITCODE -ne 0) { throw 'Coupled experiment blocked; inspect error.json in the new run directory' }

$RunRoot = Join-Path $Project "output\studies\general_scalar_extended\experimental_coupled_cached_scalar_a4d2\$RunId"
Get-Content -LiteralPath (Join-Path $RunRoot 'quartic_context_diagnostics.json')
Get-Content -LiteralPath (Join-Path $RunRoot 'summary.json')
```

The experiment runs four UV/EFT1 combinations for each model and all four W/S switch settings: 16 cases/model, 32 cases total. It is more computationally demanding than the synthetic tests.

## Verification boundary

- Source inspected on public `main` at SHA `6fa0b6426e1c3fddb4f0b31a7ae14d30b9bf7eae`.
- Synthetic cached tensor / Bose symmetry / toy coupled IVP tests run before delivery.
- The actual A(-4)/D(-2) end-to-end outputs **must be executed on your local Windows repository**; we do not have your fitted model state or Mathematica/Wolfram on this machine.
- The old experimental pipeline's result fields say `Computed` for numerically completed conditional calculations; the new output carries explicit `full_sequential_general_C5_certified=false` and `complete_physical_C5=false` markers.
- Historical baseline comparisons are diagnostics only, not a physical correction or refit.

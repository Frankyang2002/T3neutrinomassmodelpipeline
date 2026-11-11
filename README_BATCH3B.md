# Batch 3B: full real-scalar quartic linear-response construction

**New local project files**

- `Numerical/running/T3ScalarQuarticLinearResponse.py`
- `scripts/BuildT3ScalarQuarticResponseA4D2.py`
- `README_BATCH3B.md`

**Scope:** A(-4) and D(-2) at the existing fitted fermion threshold, using the repository's own source-aware quartic-tensor constructor, native SU(2) CG basis, scalar rotation, and D(-2) extra HHSS/scalar-only reference directions.

The module constructs a cache of real/imaginary response columns of the full symmetric tensor `Q[12,12,12,12] = d^4 V4 / d phi^4` and exposes `response.tensor(running_couplings)` for a fixed heavy-real-scalar basis. It runs a multi-coupling direct-reconstruction test plus two independent random complex perturbations. It does **not** change the production ODE. Basis rotation must stay fixed during EFT1 running and be handled at the scalar threshold.

**Physical status:** These columns encode the **current candidate** RGBeta/reference potential conventions. The repository explicitly does not independently certify the native-to-RGBeta scalar normalisations, D(-2) missing quadratic operator, gauge Ward consistency, full running mass Hessian, or full joint one-loop scalar matching. The script does not claim to remedy any of these. Do not label output physical or use it as a certified replacement matching result.

## PowerShell (from Downloads)

```powershell
$Project = "C:\thesistemp\T3neutrinomassmodelpipeline"
$Downloads = (New-Object -ComObject Shell.Application).Namespace('shell:Downloads').Self.Path
$Zip = Join-Path $Downloads "T3_Batch3B_QuarticLinearResponse.zip"
if (!(Test-Path -LiteralPath $Zip)) { throw "Missing ZIP: $Zip" }
Expand-Archive -LiteralPath $Zip -DestinationPath $Project -Force
Set-Location -LiteralPath $Project
python -X utf8 -m scripts.BuildT3ScalarQuarticResponseA4D2 --self-test
$Report = Join-Path $Project ("output\audits\batch3b_quartic_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
python -X utf8 -m scripts.BuildT3ScalarQuarticResponseA4D2 --project "$Project" --model both --output "$Report"
Get-Content -LiteralPath (Join-Path $Report "A_-4\report.json")
Get-Content -LiteralPath (Join-Path $Report "D_-2\report.json")
```

Expected new files:
- `$Report\A_-4\quartic_response.npz` and `$Report\A_-4\report.json`
- `$Report\D_-2\quartic_response.npz` and `$Report\D_-2\report.json`

If any required existing Batch9 scalar snapshot/Matchete seed is missing, the run intentionally fails rather than inventing data. The output directory is never overwritten; use a fresh timestamp. Upload both `report.json` files, or the traceback if either model fails. NPZ arrays need not be uploaded unless numerical tensor diagnostics are requested.

## Source basis and validation semantics

The script imports the repository's existing `two_model_master_wilson_batch13.scalar_quartic_from_matchete`, `add_d2_extra_hhss`, and `add_d2_extra_scalar_only`, with `native_triplet_basis_batch15.heavy_source_from_reference`. It uses `load_snapshot` to take a consistent coupling/mass snapshot and `build_special_basis` for D(-2) at the matching reference point. Its full tensor response is

`Q(c) = Q(c_ref) + sum_{k, Re/Im} [c_k - c_ref,k]_(component) Q_k`.

This is an implementation of efficient *evaluation* of the already defined source/reference potential, not an independent field-theory derivation of the potential. The second task (certify source/RGBeta map) is outstanding; no silent physical certification is implied by passing the interpolation tests.

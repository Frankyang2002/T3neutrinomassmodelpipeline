# T3 coupled-IVP Batch 1: read-only baseline audit

Install `scripts/AuditT3CoupledIVPBaseline.py` into the project. This script performs no RG evolution, fits, data replacement, or GitHub writes. It checks the saved 1+13+2 case locations, required bits (`1111`, `0111`, `1001`, `1110`), JSON numeric finiteness, reported status, and source hashes. It also records static A(-4) boundary indicators from the **local** runner. The CSV and JSON retain problems rather than fabricating missing results.

The code returns exit status 2 if any static check fails; inspect the reports even if it fails. By design it refuses to reuse an existing output directory. If local uncommitted versions differ from the read-only public `main`, the local files are audited and hashed; no changes are imposed.

The source read in the public repository was `scripts/RunT3SimultaneousEFT1A4.py` at `main` commit `6fa0b6426e1c3fddb4f0b31a7ae14d30b9bf7eae` (11 Nov 2026). This file passes `tree.raw * spinor_factor` as the coupled EFT1 boundary and labels complete independent fermion one-loop matching absent. This is a static observation, not physical verification.

## Project-relative files in ZIP

- `scripts/AuditT3CoupledIVPBaseline.py` (new, standalone)
- `README_BATCH1.md` (installation documentation, project root)

## PowerShell (run from any directory)

```powershell
$Project = 'C:\thesistemp\T3neutrinomassmodelpipeline'
$Downloads = (New-Object -ComObject Shell.Application).Namespace('shell:Downloads').Self.Path
$Zip = Join-Path $Downloads 'T3_Batch1_CoupledIVP_Audit.zip'
if (!(Test-Path -LiteralPath $Zip)) { throw "Missing ZIP: $Zip" }
if (!(Test-Path -LiteralPath $Project)) { throw "Missing project: $Project" }
Expand-Archive -LiteralPath $Zip -DestinationPath $Project -Force
Set-Location -LiteralPath $Project
$Report = Join-Path $Project ('output\audits\coupled_ivp_baseline_' + (Get-Date -Format 'yyyyMMdd_HHmmss'))
python -X utf8 scripts\AuditT3CoupledIVPBaseline.py --project "$Project" --output "$Report"
Get-Content -LiteralPath (Join-Path $Report 'SUMMARY.txt')
Import-Csv -LiteralPath (Join-Path $Report 'models.csv') | Format-Table -AutoSize
```

To share the results for the next implementation, upload the three generated files (`audit.json`, `models.csv`, `SUMMARY.txt`) and, if needed, the A(-4) source/export files identified in the local JSON provenance. Do not run the all-16 fit campaign yet.

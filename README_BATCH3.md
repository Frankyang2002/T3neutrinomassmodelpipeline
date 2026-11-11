# Batch 3A — A(-4)/D(-2) scalar source coverage gate

**Intended project path:** `scripts/AuditT3ScalarSectorA4D2.py`

This is a read-only source-contract implementation, not the final scalar tensor / joint one-loop threshold. It checks the actual local `general_rgbeta_quartic_registry.json`, verifies the current core source files exist, and reports unresolved quartic and quadratic directions. It optionally inspects the original Batch 2 Matchete evidence ZIP. The downstream `eft1_after_F_S1_S2_scalar_quartic_seed.json` has already integrated out both scalars; **it must not be used to reconstruct the EFT1 scalar potential**.

Expected on currently inspected public main: A(-4) zero additional quartics; D(-2) five additional quartics (two HHSS and three scalar-only), with one blocked quadratic. No correction to the coupled IVP is justified by this audit alone.

## PowerShell

```powershell
$Project = 'C:\thesistemp\T3neutrinomassmodelpipeline'
$Downloads = (New-Object -ComObject Shell.Application).Namespace('shell:Downloads').Self.Path
$Zip = Join-Path $Downloads 'T3_Batch3A_ScalarCoverage.zip'
$Evidence = Join-Path $Downloads 'T3_Batch2_Matchete_Evidence.zip'
if (!(Test-Path -LiteralPath $Zip)) { throw "ZIP missing: $Zip" }
if (!(Test-Path -LiteralPath $Evidence)) { throw "Evidence missing: $Evidence" }
Expand-Archive -LiteralPath $Zip -DestinationPath $Project -Force
Set-Location -LiteralPath $Project
$Report = Join-Path $Project ('output\audits\batch3_scalar_coverage_' + (Get-Date -Format 'yyyyMMdd_HHmmss') + '.json')
python -X utf8 scripts\AuditT3ScalarSectorA4D2.py --project "$Project" --evidence-archive "$Evidence" --output "$Report"
Get-Content -LiteralPath $Report
```

To preserve local uncommitted changes, the ZIP **adds** a new script rather than replacing the existing runners. Review potential filename collisions before `Expand-Archive -Force`. Output refuses to overwrite an existing report.

Current public main inspected: `6fa0b6426e1c3fddb4f0b31a7ae14d30b9bf7eae`, 11 November 2026.

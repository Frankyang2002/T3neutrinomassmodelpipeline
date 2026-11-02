# B1 Windows source-guard fix (2026-11-03)

This ZIP replaces **only** the B1-owned Python audit runner, and adds
a new isolated test and handoff. It does not modify the Matchete driver,
physics source, benchmark inputs, core.py, RGE, EFT, B2, or A files.

## Installation from PowerShell (repository root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Get-ChildItem "$env:USERPROFILE\Downloads" -Filter "B1_Windows_EOL_Guard_Fix*.zip" -File | Sort-Object LastWriteTime -Descending | Select-Object -First 1
if (-not $zip) { throw "Download the B1 guard fix ZIP first." }
Expand-Archive -LiteralPath $zip.FullName -DestinationPath (Get-Location).Path -Force
.\.venv\Scripts\python.exe -m unittest discover -s tests -p "test_physical_b1_*.py" -v
$native = "output\general_observable\physical_normalisation\b1_02\native"
$run = "output\general_observable\physical_normalisation\b1_eolrun_$(Get-Date -Format yyyyMMdd_HHmmss)"
.\.venv\Scripts\python.exe scripts\RunPhysicalB1WeylAudit.py --native-dir $native --output-dir $run
Get-Content (Join-Path $run "all16_component_check.json") | ConvertFrom-Json | Select-Object status,count_conditional_comparison_passed,physical_source_certified,error
```

The native JSON evidence from `b1_02/native` is reused; no need to rerun Wolfram
for this file-hashing correction. The new output directory is never overwritten.

Expected when local benchmarks and audits match: \`All16ConditionalCandidateWeylTensorsReconstructed_PhysicalDictionaryPending\` with 16 conditional comparisons, still **NOT** physical certification. Any other blockage remains a real follow-up task and must not be hidden.

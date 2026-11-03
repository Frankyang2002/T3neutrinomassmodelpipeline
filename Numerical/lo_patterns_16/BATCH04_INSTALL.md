# Batch 4 installation — local Windows PowerShell

Place `Chat_A_Batch04_Flavour_Geometry.zip` in `Downloads`, then run from a PowerShell terminal:

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip = Join-Path $HOME 'Downloads\Chat_A_Batch04_Flavour_Geometry.zip'

Expand-Archive -LiteralPath $zip -DestinationPath $repo
Set-Location $repo

.\.venv\Scripts\python.exe -m pytest -q tests\test_lo_flavour_geometry_16.py
.\.venv\Scripts\python.exe scripts\ExplainLOFlavour16.py
```

The default output is `output\general_observable\lo_pattern_study16\analysis_batch04_01\`.

For an intentional second report (the first directory is overwrite-protected):

```powershell
.\.venv\Scripts\python.exe scripts\ExplainLOFlavour16.py --output-dir output/general_observable/lo_pattern_study16/analysis_batch04_02
```

For a full local integration test (optional, normally skipped):

```powershell
$env:RUN_LO_PATTERN_INTEGRATION = '1'
.\.venv\Scripts\python.exe -m pytest -q tests\test_lo_flavour_geometry_16.py
Remove-Item Env:\RUN_LO_PATTERN_INTEGRATION
```

**Expected ordinary tests:** `17 passed, 1 skipped`, if the same source and dependencies are present. **This is not a claim that the actual local batch-4 run has already completed.**

The important file to review is `BATCH04_FLAVOUR_REPORT.md`; optional machine-readable evidence is in `batch04_flavour_analysis.json` and `batch04_mode_projections.csv`. This batch is source conditional, not B1/B2 certified.

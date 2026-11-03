# Batch 7 — install, run and inspect

New files only; unpack ZIP at repository root. Depends on already installed Batches 1, 3, 5 and 6 (and the original source audit and frozen inputs). Does not require a fresh UV run.

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip  = Join-Path $HOME 'Downloads\Chat_A_Batch07_Factorisation_Theorem.zip'
Expand-Archive -LiteralPath $zip -DestinationPath $repo
Set-Location $repo
.\.venv\Scripts\python.exe -m pytest -q tests\test_lo_factorization_theorem_16.py
.\.venv\Scripts\python.exe scripts\DeriveLOFactorisationTheorem16.py
```

Output (must be absent or empty):

`output/general_observable/lo_pattern_study16/analysis_batch07_01/`

The script generates `BATCH07_FACTORISATION_THEOREM.md`, `batch07_theorem.json`, `batch07_cg_weight_table.csv`, `batch07_first_corrections.csv`, `batch07_scenario_residuals.csv`, `SOURCE_MANIFEST.json` and `HANDOFF.json`. All these outputs are generated; no precomputed real benchmark outputs are shipped in the ZIP.

Optional rerun with a **fresh** sibling directory:

```powershell
.\.venv\Scripts\python.exe scripts\DeriveLOFactorisationTheorem16.py --output output\general_observable\lo_pattern_study16\analysis_batch07_02
```

Unit tests are self-contained synthetic checks. The skipped test is real-checkout integration by explicit opt-in only; the standalone script performs the genuine integration directly.

To inspect dominant correction errors:

```powershell
$dir = 'output\general_observable\lo_pattern_study16\analysis_batch07_01'
Import-Csv "$dir\batch07_first_corrections.csv" |
    Sort-Object { [double]$_.charge_split_linearized_relative_error } -Descending |
    Select-Object -First 8 | Format-Table -AutoSize
```

**If the script fails**, do not rerun earlier batches or modify production files. Supply the full traceback, the first few lines of the named input, and the current `git status --short`; the error is intentionally fail-closed.

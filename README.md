# Automatic analytical-report cleanup; I3 investigation

## Project files
- `scripts/BuildCombinedAnalyticalReports.py` — full replacement of saved 16-model combined report entrypoint; archives disposable files automatically after successful generation.
- `scripts/AnalyticalReportHousekeeping.py` — safe archive policy; keeps finished `.pdf`/`.tex` outputs and original raw data.
- `scripts/FindI3Definition.py` — read-only source scan for a mathematical implementation of the matching master integral `I3`.
- `tests/test_analytical_housekeeping_auto.py` — tests.

## Commands (PowerShell, repository root)

```powershell
python -m pytest -q tests/test_analytical_housekeeping_auto.py
python scripts/BuildCombinedAnalyticalReports.py --check
python scripts/BuildCombinedAnalyticalReports.py
python scripts/FindI3Definition.py .
```

The report generator now archives known LaTeX byproducts/backups and clearly superseded standalone C5 presentation assets to `Reports/report_cleanup_archive/<timestamp>/`, *automatically after a successful report-generation run*. It preserves `Reports/output/analytical/**/*.tex` and `/**/*.pdf`, aside from the specifically named superseded `C5_original_matchete.*` and `C5_original.*`. It never touches the main `output/` physics calculations.

If the combined report generator fails, it does not clean, so source/diagnostic output remains available. To inspect what would be archived before generation:

```powershell
python -c "from pathlib import Path; from scripts.AnalyticalReportHousekeeping import clean; print(clean(Path('Reports/output/analytical'), apply=False))"
```

## Critical remaining I3 issue

The checked repository report writer copies `i_3` from the saved Matchete expression and renames it to `I_3` for display; it does not implement the integral. The upstream publication verifies the general one-loop method but does not uniquely identify `i_3` in this saved T3 expression. No mathematically unjustified substitution is made here. Use `FindI3Definition.py` against your locally installed Matchete/associated Mathematica packages if needed; provide the matching definition for exact insertion in the C5 generator.

This ZIP intentionally does not claim to resolve `I3`, and does not change the C5 physics formula.

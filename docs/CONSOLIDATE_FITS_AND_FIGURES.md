# One authoritative fitted benchmark per T3 model

Source snapshot verified: GitHub main `166ed4c57b8c32eb0bb55e5b0cc9e30b67e7bb64` (2026-10-20). This package does **not** change GitHub or any physics input values.

## Project paths (full replacements / additions)

- `tools/ConsolidateFittedBenchmarks.py`: new, preview and checksum-verify 16 winning fits; optionally archive three historical stages.
- `Numerical/orchestration/studies/fitted_rge/RunFittedSensitivityComparison.py`: full replacement, authoritative `best_candidate` reads `fitted_benchmarks` and verifies SHA256 against the manifest; original `fit_stage` is retained for provenance.
- `tools/RebuildCurrentReports.py`: full replacement, sends plotted figures into `Reports/output/display/full_studies/`, uses saved numerical data, leaves original HTML destinations.
- `tools/RedirectStudyFigures.py`: new display-only `Figure.savefig` redirection for current report rebuilding.
- `tools/MoveStudyFigures.py`: new preview-first relocation of existing PNGs from `full/` to `display/full_studies/`.

## Procedure (PowerShell at repository root)

```
cd C:\thesistemp\T3neutrinomassmodelpipeline

# Preview selections BEFORE updating code. Requires the ZIP files already extracted.
python -m tools.ConsolidateFittedBenchmarks
python -m tools.ConsolidateFittedBenchmarks --apply

# Confirm 16 copied configs + one manifest:
(Get-ChildItem .\Reports\output\full\fitted_benchmarks\smallY_smallL -Directory).Count
python -m tools.RebuildCurrentReports --only displays

# Inspect current and historical figures before moving:
python -m tools.MoveStudyFigures
python -m tools.MoveStudyFigures --apply

# Rebuild plots under display/; this does not rerun any RGEs:
python -m tools.RebuildCurrentReports --only plots

# Only when the above succeeds and there are no consumer scripts relying on old stage paths:
python -m tools.ConsolidateFittedBenchmarks --apply --archive-old

git status --short
```

**Important:** do not archive old stages before running the initial copy and inspecting the manifest. The archive is under `.cleanup_quarantine/fitted_stages_YYYYMMDDTHHMMSS/`, outside `Reports/output/full/` but inside the project. Keep the quarantine until future numerical runs have been checked.

The consolidation keeps bytes unchanged for `feasibility_summary.json` and `fitted_candidate_config.json`. Historical scan diagnostics still contain old absolute fit paths; the updated fitted sensitivity reader explicitly accepts those old path strings when reusing unchanged saved scans. They are provenance, not a reason to rewrite JSON.

The plot move only handles PNGs already inside `Reports/output/full/`. Afterward, the current `tools.RebuildCurrentReports` uses `display/full_studies/<study>/<scenario>/...` for its newly created figures. The direct numerical `--run` commands in individual study modules may still emit their original plot paths inside `full/`; rerun `tools.MoveStudyFigures --apply` after such runs, or use the report rebuild.

`Reports/output/analytical/` TeX/PDF reports and all unrelated numerical outputs are untouched.

## Validation

The migration tool and PNG move were tested with 16 synthetic model directories across staged fits, with checksum-matching copies and archived original directories. Python files were syntax-checked. No full user-specific Windows checkout was available, and the end-to-end report rebuild was not run against your saved data in this environment.

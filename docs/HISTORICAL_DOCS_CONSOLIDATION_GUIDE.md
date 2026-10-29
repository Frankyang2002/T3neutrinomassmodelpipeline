# Historical T3 documentation consolidation

## Purpose

The large root `docs/` directory contains standalone research-step reports and repair notes. This package consolidates **eligible clean historical Markdown files** into up to seven topic-based volumes. It does **not** rewrite scientific claims, claim completion of conditional calculations, regenerate Mathematica artifacts, or remove unverified results.

**Important distinction:** This is a *lossless-preservation, navigational consolidation* of the original research notes, not a new synthesized physics manuscript. The full historical Markdown is retained in each volume, with local links to moved notes updated. The original raw bytes are archived externally with SHA-256 hashes, so the scientific chronology remains recoverable.

The tool reads only the user's **current local checkout**. Remote GitHub was used read-only to establish the broad source inventory, not as the authoritative source for file replacements. No GitHub writes occur.

## Source files installed by ZIP

- `tools/ConsolidateHistoricalDocs.py` — complete Python source for planner and fail-closed archive-first local migration.
- `tests/python/test_historical_docs_consolidation.py` — synthetic regression coverage for link classification, transitive reference blockers, dirty-file protection, original-byte preservation, and saved-data non-modification.
- `docs/HISTORICAL_DOCS_CONSOLIDATION_GUIDE.md` — these instructions; it is a permanent guide, not a retirement candidate.

## Protected documents and data

The following primary docs remain individual files: `PHYSICS.md`, `SCALAR_SECTOR.md`, `NUMERICAL.md`, `RESULTS.md`, `T3_THESIS_PLAN.md`, `REGENERATE_CURRENT_REPORTS.md`, `THESIS_FIGURE_PLACEMENT.md`, `THESIS_FIGURE_CAPTIONS.md`, `THESIS_EXISTING_COMPARISONS_PLACEMENT.md`, `GENERAL_SCALAR_48_FINAL_REPORTING.md`, `T3_GENERAL_SCALAR_STUDY_SEPARATION.md`, `FERMION_START_STUDY.md`, and `T3B_YUKAWA_BASIS.md`.

**No code or experimental data is changed:** `scripts/`, `Numerical/`, `Lagrangian/`, `RGE/`, `tests/`, `configs/`, `output/`, `Reports/`, scientific JSON/CSV/TEX, fitted inputs, provenance JSON and Matchete/RGBeta outputs are untouched. Untracked/locally modified documentation is not moved. Old instructions embedded in the consolidated volumes may be superseded; use current master documents and source as the authority.

## Thematic volumes generated *on the local checkout*

Depending on the eligible documents and references, the tool generates a subset of:

- `docs/RESEARCH_GENERAL_T3_TENSORS_AND_QUARTICS.md`
- `docs/RESEARCH_GENERAL_T3_RUNNING_AND_MATCHING.md`
- `docs/RESEARCH_MATCHETE_AND_MASS_EIGEN.md`
- `docs/RESEARCH_EFT1_THRESHOLDS_AND_WEINBERG.md`
- `docs/RESEARCH_NUMERICAL_FITS_AND_BENCHMARKS.md`
- `docs/RESEARCH_GENERAL_SCALAR_AND_PILOTS.md`
- `docs/REPOSITORY_MAINTENANCE_HISTORY.md`
- `docs/RESEARCH_HISTORY_INDEX.md` (always generated if there is something eligible to move)

Each volume has a table of contents and original source filename and SHA-256 for every embedded research record, with the original text expandable per section. Clickable links between moved Markdown records are retargeted within the volumes; commands in historical fenced code blocks are retained unchanged.

If a retained source file or test references a candidate filename, that original file is **kept**. This check follows chains of references among retained docs. Plain narrative filename mentions in Markdown are logged without automatically blocking; executable/config references and active Markdown links are blockers. Dynamic paths, binary files and references recorded only inside `output/` remain outside the source-reference proof; original paths are still preserved by the external archive.

## Run on Windows PowerShell

Use the active project `.venv`:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

python -X utf8 -m pytest -q tests/python/test_historical_docs_consolidation.py --tb=line
if ($LASTEXITCODE -ne 0) { throw 'Historical docs tests failed' }

# Read-only plan. Writes reports under output/cleanup_audit only.
python -X utf8 -m tools.ConsolidateHistoricalDocs
if ($LASTEXITCODE -ne 0) { throw 'Historical docs audit failed' }
Get-Content -Encoding UTF8 .\output\cleanup_audit\docs_consolidation_v1\summary.md

# Inspect all classifications, including active-reference blockers.
Import-Csv .\output\cleanup_audit\docs_consolidation_v1\inventory.csv |
    Group-Object classification | Select-Object Name, Count
```

Do not apply unless `summary.md` lists sensible totals and the blocked/retained files make sense. *The intended end-state is fewer docs, but the tool will not force a target count by moving referenced files.*

```powershell
python -X utf8 -m tools.ConsolidateHistoricalDocs `
    --apply --confirm CONSOLIDATE_HISTORICAL_DOCS
if ($LASTEXITCODE -ne 0) { throw 'Docs consolidation failed; do not commit' }

python -X utf8 -m pytest -q --tb=line
if ($LASTEXITCODE -ne 0) { throw 'Regression tests failed; do not commit' }

git status --short
git diff --check
```

Do not rerun `--apply` after a successful consolidation: the created volume and index filenames are deliberately treated as existing destinations and will trigger a fail-closed stop.

## External archive and recovery

Application creates a dated directory **outside** the project, for example:

`C:\thesistemp\T3_historical_docs_archive_20261029T...Z\`

It preserves original project-relative paths such as `docs/GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md` and verifies the byte-for-byte SHA-256 before deleting any original. The archive contains `manifest.json` identifying all original files, corresponding volumes and their hashes. Generated new volumes are verified against the input contents before originals are removed. In-process failure triggers original-file restoration from the archive and removal of generated volumes.

For individual restoration, use `Copy-Item` from `<archive>\docs\<filename>` back to the original project path. Avoid restoring files after the changes have been committed without first inspecting references and Git status.

## Interpretation limits

A green Python test suite checks imports and assertions, not live Wolfram/Matchete/RGBeta reproducibility. This consolidation makes **no physical certification claim** and modifies no mathematical source implementation. The external archive and Git history are required to reconstruct historical filenames appearing in saved provenance records.

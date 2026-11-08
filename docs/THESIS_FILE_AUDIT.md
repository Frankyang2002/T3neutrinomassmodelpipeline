# File-by-file thesis relevance and dependency audit

## Scope

This review is intended for the **current local checkout**, not for rewriting Git history. The initial source inspection used GitHub `main` commit `806cdf08608dbe3f5e0dcd8af785cfad16f26af6` (2026-11-08). The latest audited checkout HEAD is recorded automatically in each generated plan; new changes are checked locally rather than assuming the old commit still applies.

The project has roughly 9,800 tracked files, including more than 550 files under `scripts/`, 18 historical `tools/` utilities, and thousands of saved numerical and symbolic outputs. A single filename or absence of static imports **cannot establish scientific irrelevance**. The local auditor therefore distinguishes inventory from archival approval.

## What the audit does

`tools/ThesisRepositoryAudit.py` runs without nonstandard Python dependencies. It does not import or execute research scripts. It inventories **every Git-tracked path** and visible untracked paths and assigns a directory/content role. For every research script and utility, it records:

- Original relative path, SHA-256, size, role and a brief source docstring/description.
- A research topic derived from the filename and its source location.
- Incoming textual references from executables, tests, configuration, saved results and documentation; the CSV records exact source files and lines.
- Whether Git reports the file unchanged, and whether it belongs to a protected workflow.
- A decision: `KEEP`, `KEEP_REFERENCED`, `KEEP_PROVENANCE`, `HOLD_LOCAL_CHANGE`, `HOLD_INCOMPLETE_SCAN`, `REVIEW_DOCUMENTED`, `REVIEW_STANDALONE`, or `SUGGEST_ARCHIVE`.

The scanner reads the text of tracked and locally untracked files, not only the `scripts/` folder. It ignores previously generated audit reports and auto-generated script/document indexes: an index listing every script is a navigation aid, not evidence that all scripts are executed.

Static scanning has limits. Dynamic imports, generated Wolfram calls, external notebooks and private local commands may not be discoverable. `REVIEW_STANDALONE` means **human review required**, not unused code. An oversized/unreadable text input disables automatic archival recommendations until coverage is resolved.

## Run the full file inventory

From `C:\thesistemp\T3neutrinomassmodelpipeline` in PowerShell:

```powershell
python -X utf8 tools/ThesisRepositoryAudit.py

$dir = Join-Path (Split-Path (Get-Location) -Parent) 'T3_external_archives/thesis_repository_review'
Import-Csv "$dir/scripts_and_tools.csv" |
    Group-Object status | Sort-Object Count -Descending |
    Format-Table Count,Name -AutoSize

Import-Csv "$dir/scripts_and_tools.csv" |
    Where-Object { $_.status -eq 'SUGGEST_ARCHIVE' } |
    Select-Object path, topic, reason, doc_refs |
    Format-Table -AutoSize

Import-Csv "$dir/scripts_and_tools.csv" |
    Where-Object { $_.status -eq 'REVIEW_STANDALONE' } |
    Select-Object path, description |
    Out-GridView
```

Output files are stored in the sibling `C:\thesistemp\T3_external_archives\thesis_repository_review` folder, **outside Git**. They are not automatically added to the repository or GitHub.

Output files:

| File | Description |
| --- | --- |
| `C:\thesistemp\T3_external_archives\thesis_repository_review\all_files.csv` | Every tracked file plus discovered untracked file, with roles and byte counts |
| `C:\thesistemp\T3_external_archives\thesis_repository_review\scripts_and_tools.csv` | Per-file classified inventory for scripts and tools |
| `C:\thesistemp\T3_external_archives\thesis_repository_review\reference_witnesses.csv` | Exact incoming reference evidence, separated by source type |
| `C:\thesistemp\T3_external_archives\thesis_repository_review\plan.json` | HEAD-bound plan, SHA-256 values and classification evidence |
| `C:\thesistemp\T3_external_archives\thesis_repository_review\summary.md` | Count summary and suggested candidates |
| `C:\thesistemp\T3_external_archives\thesis_repository_review\reviewed_selection.txt` | Explicit archive selection, initially limited to unreferenced one-time utilities and verified junk |

The explicitly inspected one-time *script* candidates are `scripts/InstallT3RGBetaDoubletMassFix.py` (the heavy-doublet epsilon fix is present in inspected `T3RGBetaModel.wl`), `scripts/RenameT3ModelFolders.py`, `scripts/RemoveOldGeneralScalarDisplay.ps1`, `scripts/RetireOldStudyDisplays.ps1`, `scripts/RetireSupersededT3HTMLGenerators.ps1` and `scripts/maintenance/Remove-LegacyReadmes.ps1`. They are only *candidates*: code/test or saved-output references hold them, and documented calls require explicit review. No general script-name rule produces archive approval.

**Inspect the CSV and evidence before removing anything.** Do not select numerical study scripts solely because they are named `Run...Batch`, `...Fix`, `...Old` or `.experimental.py`. Several preserve unique EFT matching, basis and covariance checks.

## Reversible external archive

After reviewing the selection file, archive exactly those paths:

```powershell
python -X utf8 tools/ThesisRepositoryAudit.py `
    --apply `
    --selection "C:\thesistemp\T3_external_archives\thesis_repository_review\reviewed_selection.txt" `
    --confirm ARCHIVE_REVIEWED_THESIS_FILES

.\scripts\RunCoreTests.ps1

git status --short
git diff --check
```

For a standalone research script whose source and numerical role have been *individually checked*, add its exact path to the selection file, and explicitly supply `--allow-research-scripts`. To archive a path still cited **only in documentation**, first inspect and update the affected documentation; then explicitly add `--allow-documented`.

The archival operation rejects Git-HEAD drift, modified files, incoming executable/test/config references, saved-result references, changes to reference evidence since preview and mismatching Git blobs. All selected bytes are copied outside the checkout and SHA-256 verified **before** the working copies are removed. The archive includes an original-path manifest. If a failure interrupts file removal, retain the archive for recovery.

To restore, substitute the exact verified path printed by the archival command:

```powershell
python -X utf8 tools/ThesisRepositoryAudit.py --restore `
    'C:\thesistemp\T3_external_archives\thesis_review_YYYYMMDDTHHMMSSffffffZ'
```

Restoration refuses to overwrite any existing path. The script does not stage, commit, push or rewrite Git history.

## What should remain in the thesis codebase

- Physical model and interaction construction (`common/`, `model/`, `Lagrangian/`).
- Representation algebra, UV RGEs and intermediate scalar EFT (`RGE/`).
- Numerical state, fermion and scalar threshold interfaces, SMEFT running, fitting and benchmarks (`Numerical/`, `configs/`, `physics/`).
- Tests of symmetry, flavour, sign, normalisation, heavy-fermion decoupling, scalar decoupling and numerical observables.
- Reproducible, source-backed figures and saved configurations used in the thesis. Preserve original inputs and hashes even for conditional results.

The calculation's required order remains

\[
\mathrm{UV}\;\xrightarrow{\text{integrate out all heavy fermions}}\;
\mathrm{EFT}_{\text{SM+scalars}}\;\xrightarrow{\text{integrate out all heavy scalars}}\;
\mathrm{SMEFT}.
\]

Do **not** label an inventory of general interactions as fully matched or a conditional one-scale Weinberg result as the completed fermion-first two-threshold EFT. The primary fixed-boundary reduced-scalar fits and the exploratory general-coupling studies have different physical certification scopes; consult `docs/RESULTS.md`, `docs/SCALAR_SECTOR.md` and `docs/GENERAL_SYMBOLIC_PHYSICS_AUDIT.md`.

## Follow-up after first archival pass

After reviewing `scripts_and_tools.csv`, use the report to select historical study drivers that no longer support the thesis. Do this in small reviewed sets, rerunning the 189-test core suite each time. **Do not mass-archive `output/`**: saved native matching inputs, fitted JSONs, provenance and figure source data can be runtime dependencies. A subsequent separate audit can trace each retained thesis figure/benchmark to the exact saved data it needs.

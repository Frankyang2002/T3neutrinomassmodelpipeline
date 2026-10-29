# Repository Maintenance History

> Historical research record. This document combines **original source
> Markdown** under original-path headings. Local Markdown links to moved
> files are updated; exact original bytes are in the external archive.
> This is not a new certification
> of physics, a rewritten conclusion, or a new computational workflow.
> Original bytes and SHA-256 values are stored in the external archive.

## Contents

- [Descriptive filenames for numbered T3 research stages](#source-descriptive-research-filenames-md) — `docs/DESCRIPTIVE_RESEARCH_FILENAMES.md`
- [Retain only the selected fitted candidates](#source-fit-candidate-retention-md) — `docs/FIT_CANDIDATE_RETENTION.md`
- [Retiring one-time root README notes](#source-legacy-notes-cleanup-md) — `docs/LEGACY_NOTES_CLEANUP.md`
- [Filename migration: historical Stage 23 report provenance](#source-renamed-quartic-provenance-compatibility-md) — `docs/RENAMED_QUARTIC_PROVENANCE_COMPATIBILITY.md`
- [Reports/output/full — retention audit](#source-reports-full-retention-audit-md) — `docs/REPORTS_FULL_RETENTION_AUDIT.md`
- [Historical script dependency audit (read-only)](#source-script-dependency-audit-md) — `docs/SCRIPT_DEPENDENCY_AUDIT.md`
- [One-time artifact cleanup (local, archive-first)](#source-one-time-artifact-cleanup-md) — `docs/ONE_TIME_ARTIFACT_CLEANUP.md`
- [T3 repository cleanup: historical-document reference audit](#source-repository-cleanup-audit-md) — `docs/REPOSITORY_CLEANUP_AUDIT.md`
- [Script dependency audit — refinement of the 77 broad-reference modules](#source-script-dependency-refinement-md) — `docs/SCRIPT_DEPENDENCY_REFINEMENT.md`

---

<a id="source-descriptive-research-filenames-md"></a>

## Original document: `docs/DESCRIPTIVE_RESEARCH_FILENAMES.md`

**Original SHA-256:** `24f85daaa2b71a7e0a1e0dbe81b93d71d0f2aa1ea22f12bd6a74d2258e5f9b5d`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Descriptive filenames for numbered T3 research stages</summary>

# Descriptive filenames for numbered T3 research stages

## Scope and source

Read-only inspected baseline: Git commit `fe0f511104cce8d7d914d114e3fd884eb47c6a95` (29 October 2026). This migration executes **only on the user's local checkout**. It does not access or modify GitHub. The supplied migration is a **full standalone Python file**, not a patch. It generates renamed copies of the full original files from the user's own checkout.

The existing numbered research names refer to genuinely different calculations. The rename changes **source filenames and their exact references**, not the calculations themselves.

| Original source | Descriptive filename |
|---|---|
| `scripts/Stage30CoupledA0UVRunning.py` | `scripts/CoupledA0UVRunning.py` |
| `scripts/Stage31CoupledUVRunning.py` | `scripts/CoupledUVRunning.py` |
| `scripts/Stage33FermionThresholdHandoff.py` | `scripts/FermionThresholdHandoff.py` |
| `scripts/Stage34IntermediateEFTRunning.py` | `scripts/IntermediateEFTRunning.py` |
| `scripts/Stage35ScalarThresholdHandoff.py` | `scripts/ScalarThresholdHandoff.py` |
| `scripts/Stage36FinalWeinbergRunning.py` | `scripts/FinalWeinbergRunning.py` |
| `scripts/RunGeneralT3Stage16.py` | `scripts/RunGeneralT3FermionFirstThresholds.py` |
| `scripts/RunGeneralT3Stage17.py` | `scripts/RunGeneralT3EFT1ScalarPairMixing.py` |
| `scripts/RunGeneralT3Stage18.py` | `scripts/RunGeneralT3EFT1WilsonSourceTransport.py` |
| `scripts/RunGeneralT3Stage19.py` | `scripts/RunGeneralT3WeinbergSourceTransport.py` |
| `scripts/Stage22TensorAlgebra.py` | `scripts/RGBetaTensorAlgebra.py` |
| `docs/GENERAL_T3_STAGE17_STATUS.md` | `docs/GENERAL_T3_EFT1_SCALAR_PAIR_MIXING.md` |
| `tests/python/test_general_t3_stage36_final_weinberg.py` | `tests/python/test_general_t3_final_weinberg.py` |

The full mapping of all eligible `StageNN` filenames (source scripts, PowerShell/Wolfram launchers, research notes, tests) is generated from actual **tracked files** in the local checkout. The tool refuses naming collisions or unexpected target paths. Supporting source references to the exact module, script and document filenames are rewritten locally.

### Intentionally excluded

- `output/**` and `Reports/output/**`, especially **saved scientific results**, historical stage-named RGE/matching exports and original numerical provenance.
- `docs/archive/**` and test fixtures, whose content is historical evidence.
- Batch-numbered filenames, since `Batch` labels are not part of this request.
- Narrative historical labels (`Stage 16`, `Stage 17`), physics variable names, and old *saved-data* paths (`output/general_t3_stage17/...`): the chronology and data contracts must remain reconstructible.
- JSON inputs and exported data, which must retain their original source fingerprints.

**Important:** Because some historical JSON files record script or documentation paths, renaming the active scripts cannot retroactively make every old provenance path valid. Original source files are preserved in the external backup and in Git history; do not rewrite saved numerical provenance to hide those historical references. This migration therefore does not certify old artifact replay. Regenerate reports from current paths only after running tests and checking source hashes.

## Local PowerShell procedure

Start from the repository root and use the existing Python virtual environment:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

# Optional: verify audit tool tests first.
python -X utf8 -m pytest -q tests/python/test_descriptive_research_file_names.py --tb=line

# Read-only plan, with a machine-readable mapping saved under output/.
python -X utf8 -m tools.RenameNumberedResearchFiles `
    --report output/cleanup_audit/descriptive_research_filenames/plan.json

# Review this JSON before applying.
Get-Content -Encoding UTF8 output/cleanup_audit/descriptive_research_filenames/plan.json |
    Select-Object -First 75

# Apply (modifies local checkout; NO GitHub writes).
python -X utf8 -m tools.RenameNumberedResearchFiles `
    --apply --confirm RENAME_SOURCE_FILES
if ($LASTEXITCODE -ne 0) { throw "Rename failed; inspect backup and error" }

# All existing tests; no symbolic recalculation or optimization.
python -X utf8 -m pytest -q --tb=line
if ($LASTEXITCODE -ne 0) { throw "Tests failed; do not commit rename" }

git status --short
```

The command rejects a changed HEAD, source files with uncommitted local modifications, existing rename destinations and a pre-existing backup location. It permits unrelated uncommitted deletions from the previous documentation cleanup. **Do not run `git reset --hard`** to make the rename work.

## Backup and restoration

Before editing anything, the tool copies every affected original file to the directory adjacent to the project:

`C:\thesistemp\T3_stage_filename_backup_fe0f511\originals\...`

The external backup's `manifest.json` records the SHA-256 of each original, new destination and rewritten bytes. The copy is verified before source renaming begins. On a failed application, the tool attempts to restore the original files automatically and leaves the backup intact.

To restore manually, first inspect `manifest.json`, remove only the renamed destinations listed under `renames`, and then copy the contents of `originals/` back into the project root while preserving paths. Do not delete unrelated user changes. If tests fail, retain the backup and do not commit until the failure is diagnosed.

## Testing scope

The package's seven tests cover naming maps, Python import/CLI path rewrites, preservation of saved-data paths, document titles, and full backup/application in a synthetic Git checkout. They are **not** a replacement for the Windows repository's full Python tests or external Wolfram/Matchete/RGBeta regression. Paths constructed dynamically from stage numbers may require manual follow-up if the full tests reveal them.

### Files added by this package

- `tools/RenameNumberedResearchFiles.py` — standalone complete local migration tool.
- `tests/python/test_descriptive_research_file_names.py` — seven offline tests.
- `docs/DESCRIPTIVE_RESEARCH_FILENAMES.md` — this guide.

No previously archived documents, physics scripts, saved data, numerical configurations or GitHub files are included in the ZIP.

</details>

---

<a id="source-fit-candidate-retention-md"></a>

## Original document: `docs/FIT_CANDIDATE_RETENTION.md`

**Original SHA-256:** `f16d391b1de7086e74f83895f58cc7a331b47679db3c97ba7f30e896813e85d1`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Retain only the selected fitted candidates</summary>

# Retain only the selected fitted candidates

Audited against public repository main commit `166ed4c57b8c32eb0bb55e5b0cc9e30b67e7bb64` (20 Oct 2026). **Local-only**. Does not update GitHub, recompute fits, or modify physics.

The active `best_candidate()` in `Numerical/orchestration/studies/fitted_rge/RunFittedSensitivityComparison.py` scans three fitting stages in order:

1. `fit_feasibility`
2. `fit_feasibility_retry`
3. `fit_feasibility_focused`

For each model it selects the candidate with the **lowest finite `best.chi2`** from its `feasibility_summary.json`, requiring an adjacent `fitted_candidate_config.json`. No full fitting-stage folder is redundant.

## Best candidates in the inspected saved data

| Selected stage | Models | Count |
|---|---|---:|
| Original (`fit_feasibility`) | T3A-m4, T3A-m2, T3A-p0; T3B-m3, T3B-p1; T3D-m2, T3D-p0, T3D-p2 | 8 |
| Retry (`fit_feasibility_retry`) | T3B-m1; T3C-m3, T3C-p1; T3E-m4, T3E-m2, T3E-p0, T3E-p2 | 7 |
| Focused (`fit_feasibility_focused`) | T3C-m1 | 1 |

**Nine superseded model-level folders** can be archived if your local numerical results still agree with the above: original fit folders for B-m1, C-m1, C-m3, C-p1, E-m4, E-m2, E-p0, E-p2; plus the retry folder for C-m1.

## Install and run in PowerShell

From the project root, after placing the ZIP next to `pipeline.py`:

```powershell
Expand-Archive .\T3_FIT_CANDIDATE_CLEANUP.zip -DestinationPath . -Force
Unblock-File .\scripts\ArchiveSupersededFits.ps1
.\scripts\ArchiveSupersededFits.ps1                # preview only
.\scripts\ArchiveSupersededFits.ps1 -Apply -WhatIf # PowerShell no-op simulation
.\scripts\ArchiveSupersededFits.ps1 -Apply         # move, never permanently delete
```

The archive is under `.cleanup_quarantine/superseded_fit_candidates/<UTC timestamp>/Reports/output/full/...` with an `archive_manifest.csv` identifying original and retained stage/chi2. The script recomputes winners from your actual **local** summary JSON, requires 16 unique valid models, leaves any unexpected directory contents alone, and checks for a conflicting saved `fitted_rge` stage before archiving that model's old candidates.

**Do not archive `fit_feasibility`, `fit_feasibility_retry`, or `fit_feasibility_focused` in their entirety.** The main rebuild and other study runners call `best_candidate()` and require the selected candidates at their original paths.

Do not modify `batch_summary.json` or `retry_summary.json`: these are historical summaries and are not re-generated by a presentation-only report rebuild. Old `previous_summary`/`previous_retry_summary` strings in retry/focused summary JSON can refer to previous attempts' original locations; the archive manifest preserves the new locations. The script does not rewrite numerical data or provenance fields.

The cleanup is **reversible** by moving the archived model directory back to its original path. The script has no deletion command. It does not remove unrelated empty directories; use the previously supplied full-report empty-folder cleanup for those.

After the move, verify the live fit selection using the existing plotting/rebuild workflows without running `--run` on fitting or RGE study generators. Review `git status --short` before any commit; no GitHub changes are made.

</details>

---

<a id="source-legacy-notes-cleanup-md"></a>

## Original document: `docs/LEGACY_NOTES_CLEANUP.md`

**Original SHA-256:** `1e2c9faf458326e9efea4f0147f456f89548efb4b3586703f4f75a89993ab89d`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Retiring one-time root README notes</summary>

# Retiring one-time root README notes

This is a conservative second-pass audit for `README_*.txt` in the repository root.
It is independent of the earlier console-log archival. Its only automatic
retirement candidates are **unreferenced, Git-tracked, unmodified** root README
notes. They are copied and SHA-256-checked in an **external archive** before
removing the repository originals. It does not rewrite any scientific result.

`README_THESIS_COMPARISONS.txt` is protected: it describes thesis comparison
assets. `README.txt`, `README.md`, `SHA256_MANIFEST.json`, all three root
scientific-output ZIPs, the four legacy PowerShell utilities, the two existing
`report_cleanup_archive/` folders, all source code, and all generated
numerical and report outputs are not retired by this tool.

From the repository root:

```powershell
python -X utf8 -m pytest -q tests/python/test_legacy_notes_cleanup.py --tb=line
python -X utf8 -m tools.AuditLegacyNotes
Get-Content -Encoding UTF8 .\output\cleanup_audit\legacy_notes_v2\summary.md
```

Inspect the inventory. To archive **only** the marked `ARCHIVE_CANDIDATE` notes:

```powershell
python -X utf8 -m tools.AuditLegacyNotes --apply --confirm ARCHIVE_UNREFERENCED_LEGACY_NOTES
if ($LASTEXITCODE -ne 0) { throw 'Legacy note archive failed.' }
python -X utf8 -m pytest -q --tb=line
```

A timestamped `T3_legacy_notes_archive_*` folder is created in the project
parent folder, outside Git, with a manifest recording original project paths,
SHA-256 values and project HEAD. The tool refuses an existing archive target.

The text-reference scanner is deliberately conservative: a literal mention in
source documentation blocks automatic archival even when it is only historical.
It does not scan inside binary scientific ZIP bundles or retained numerical
output folders. These are explicitly **protected**. Review the remaining
historical PowerShell utilities and output bundles independently; do not infer
that a lack of source imports makes them safe to delete.

</details>

---

<a id="source-renamed-quartic-provenance-compatibility-md"></a>

## Original document: `docs/RENAMED_QUARTIC_PROVENANCE_COMPATIBILITY.md`

**Original SHA-256:** `96b4be315322be3d1f86a5162068771b837934828cea11d41370a107ddb22ea6`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Filename migration: historical Stage 23 report provenance</summary>

# Filename migration: historical Stage 23 report provenance

The research filename migration changes the Matchete tensor export script from
`scripts/ExportGeneralT3Stage22Matchete.wl` to
`scripts/ExportGeneralT3ScalarTensorMatchete.wl`. The original saved report
`output/general_t3_stage22/stage23_quartic_prefactor_candidates.json` still
records the old pathname in `source_verified.matchete_source` as historical
provenance. The current `scripts/QuarticPrefactorAudit.py` correctly describes
its renamed source instead.

`NativeQuarticActionAudit.py` now permits precisely this **single string
metadata difference** while continuing to:

- recalculate the entire Stage 23 report from Stage 22 inputs;
- validate their byte-level SHA-256 hashes via `stage23_run`;
- require exact equality for every other saved Stage 23 field;
- reject missing or tampered physics arrays, tensor blocks and status flags;
- perform the original native RGBeta assembly verification.

No saved JSON or physical source is rewritten. A failed hash check or any other
change to the candidate report still causes a fail-closed error. This is a
provenance-path compatibility exception, **not** a new Matchete/RGBeta coupling
certification.

## Local verification

```powershell
python -X utf8 -m pytest -q tests/python/test_renamed_quartic_action_provenance.py --tb=line
python -X utf8 -m pytest -q tests/python/test_general_t3_native_quartic_action.py tests/python/test_general_t3_quartic_action_identity.py tests/python/test_general_t3_triplet_intertwiner.py tests/python/test_general_t3_triplet_native_structure.py --tb=line
python -X utf8 -m pytest -q --tb=line
```

</details>

---

<a id="source-reports-full-retention-audit-md"></a>

## Original document: `docs/REPORTS_FULL_RETENTION_AUDIT.md`

**Original SHA-256:** `3b5c794c04f5bcc7c9ee8ea4e0d505d97773bd52f3352d81f4b8b52246ad639b`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Reports/output/full — retention audit</summary>

# Reports/output/full — retention audit

**Repository basis:** read-only inspection of GitHub `main` at commit `166ed4c57b8c32eb0bb55e5b0cc9e30b67e7bb64` on 2026-10-20. This is **not** a scan of the user's Windows filesystem. Local untracked files/empty folders may differ.

## Supported physics input families

1. **Common-input:** UV ordinary couplings, heavy parameters and scales inherited from fitted `T3-A-m2` by `RunCommonInputRGE.py`, including the explicitly coded class-E half-`lambdaT3` convention.
2. **Individual fits:** independently optimised Yukawa matrices/parameters selected from feasibility, retry or focused fit candidates. Some studies share a scanned scalar coupling across these fits; that does **not** make them obsolete E-p0/common-Yukawa studies.

## Retain all current tracked numerical study families

| Under `Reports/output/full/` | Tracked files | Basis / why keep |
|---|---:|---|
| `common_input_rge/` | 145 | Common A-m2 reference: diagnostics and running figures |
| `common_input_sensitivity/` | 17 | Common-reference gamma scans used by current envelope HTML |
| `common_input_coupling_scan/` | 28 | Common-reference coupling variations, RGE/frozen results |
| `common_input_fermion_scale/` | 23 | Common-reference fermion mass/threshold scan |
| `fit_feasibility/` | 34 | Fitted configurations and fit selection inputs |
| `fit_feasibility_retry/` | 18 | Best-fit selection inputs (not obsolete solely because called retry) |
| `fit_feasibility_focused/` | 2 | Best-fit selection inputs (not obsolete solely because called focused) |
| `fitted_rge/` | 131 | Individually fitted running diagnostics and current figures |
| `fitted_sensitivity/` | 56 | Active wide-range fitted gamma scan and thesis figures |
| `fermion_scale/` | 7 | Individually fitted mass-ratio results; current renderer reads `scan.json` |
| `staged_rg_comparison/` | 2 | **Retain** saved 000,100,110,111 numerical endpoint comparisons; historical internal folder name |
| `weinberg_running_comparison/` | 2 | **Retain** saved Weinberg evolved/frozen endpoint comparison data |
| `common_scalar_comparison/` | 28 | **Individual fitted Yukawas**; a common lambda grid, not four historical benchmarks |
| `common_scalar_rg_controls/` | 12 | Running vs frozen control points corresponding to common scalar scans |
| `common_scalar_beta_dominance/` | 5 | Derived UV beta decomposition linked to current scalar study |
| `heavy_fermion_scale_comparison/` | 12 | **Individual fitted Yukawas** plus common-scalar points; retain numerical Study B data |
| `representation_comparison/` | 2 | Compact aggregation of fitted UV/EFT1/Weinberg endpoint contributions |
| `common_scalar_diagnostics/` | 2 | 80-point scalar study audit/summary; not used by current rebuild, but real data |
| `fitted_sensitivity_all/` | 18 | Earlier narrow `gamma=0.8...1.2`, 21-point fitted scan; not used by current rebuild, but distinct saved numerical samples; archive only if intentionally discarding that scan |

For the 19 tracked groups: **544 files total; zero 0-byte files and zero empty tracked directories**. These statistics do not say whether a Windows checkout has additional untracked/ignored files or empty directories.

## Retire (if present locally)

The four fixed benchmark and old shared-parameter/E-p0 report lineage:

- `four_benchmark_correlations/`
- `common_parameter_comparison/`, `common_parameter_comparison_v2/`
- `common_yukawa_optimization/`, `common_yukawa_least_squares/`, `common_yukawa_final_summary/`
- `common_yukawa_correlation_scans/`, `common_yukawa_three_panel_scans/`
- `scalar_coupling_correlations/`
- `benchmark_summary/`, `yukawa_grouping/`
- `comparison/`, `within_model_comparison/` (old umbrella/four-scenario lineage)

These 13 directory names are absent from the audited GitHub full tree; they may still exist in a local checkout. The accompanying script will **quarantine**, not irreversibly delete, those directories when `-Apply` is explicitly requested.

Also remove old `staged_rg_*.png` from `staged_rg_comparison/<scenario>/` and `weinberg_running_*.png` from `weinberg_running_comparison/<scenario>/` **only** when the corresponding five current display PNGs exist and are nonempty. Do not delete these studies' `scan.json` and `scan.csv`.

## Manual review (not automatic removal)

If found locally: `beta_dominance/`, `observable_sensitivity/`, `common_scalar_pilot/`, `hypercharge/`, `dimensions/`, `smoke/` and unknown study roots. The current main rebuild doesn't use them, but absence from the main workflow alone does not prove all contained data/artifacts are disposable. Check for unique saved numerical or symbolic content first.

Do **not** automatically remove `normalized_fitted_sensitivity/`: `tools/RebuildCurrentReports.py` currently generates it from `fitted_sensitivity/`, even though another historical cleanup script names it for retirement. Keeping or moving that derived output is a separate report-layout decision.

**Avoid running `scripts/CleanReportOutput.ps1 -Apply` for this selective job.** It intentionally deletes all PNGs and HTML under `Reports/output`, including actively retained figures, and its scope is much broader.

## Install and run (PowerShell at repository root)

```powershell
Expand-Archive -LiteralPath .\T3_REPORTS_FULL_AUDIT.zip -DestinationPath . -Force
.\scripts\AuditAndPruneReportsFull.ps1
.\scripts\AuditAndPruneReportsFull.ps1 -Apply -WhatIf
.\scripts\AuditAndPruneReportsFull.ps1 -Apply
```

The script is preview-only by default; deletion requires explicit `-Apply` (add `-Confirm` for individual confirmations). Retired trees are moved to `.cleanup_quarantine/reports_full/<timestamp>/` outside the report tree. Empty files are deleted and empty directories are removed bottom-up. The script does **not** delete unrecognized nonempty directories, current numerical studies, fitted benchmark JSON/CSV, or files outside `Reports/output/full` (aside from creating the quarantine for recoverability).

For an independent exact local inventory:

```powershell
Get-ChildItem .\Reports\output\full -Directory | Select-Object Name
Get-ChildItem .\Reports\output\full -Recurse -File -Force |
    Where-Object Length -eq 0 | Select-Object FullName
Get-ChildItem .\Reports\output\full -Recurse -Directory -Force |
    Where-Object { -not (Get-ChildItem -LiteralPath $_.FullName -Force | Select-Object -First 1) } |
    Select-Object FullName
```

</details>

---

<a id="source-script-dependency-audit-md"></a>

## Original document: `docs/SCRIPT_DEPENDENCY_AUDIT.md`

**Original SHA-256:** `f27a060e33585829294ae50c3fb91f057ab6f5bf410fc131ef768f4c21dbaedf`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Historical script dependency audit (read-only)</summary>

# Historical script dependency audit (read-only)

**Scope:** identification of direct and indirect dependencies before any possible archiving of old numerical or symbolic scripts. This tool does **not** delete, rename, archive, import, or execute project physics scripts. It only reads repository text and writes an audit report.

At the inspected GitHub commit `fe0f511104cce8d7d914d114e3fd884eb47c6a95` (29 October 2026), `scripts/` contains 99 Python files, including `__init__.py`. The tool inventories all non-`__init__` top-level Python script modules **from your own local checkout**, not from a hard-coded snapshot.

## Installation

The package adds these two files, preserving project-relative paths:

- `tools/AuditScriptDependencies.py`
- `docs/SCRIPT_DEPENDENCY_AUDIT.md`

From the local repository root in PowerShell:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Join-Path $HOME 'Downloads\T3_Script_Dependency_Audit.zip'
if (!(Test-Path $zip)) { throw "ZIP not found: $zip" }
if ((Test-Path .\tools\AuditScriptDependencies.py) -or (Test-Path .\docs\SCRIPT_DEPENDENCY_AUDIT.md)) { throw 'Audit file(s) already exist. Review before overwriting.' }
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
python -X utf8 -m tools.AuditScriptDependencies
Get-Content -Encoding UTF8 .\output\cleanup_audit\script_dependencies_v1\summary.md | Select-Object -First 55
```

The script does not require Wolfram, Matchete, RGBeta, a live internet connection, or numerical reruns. It does not run any `scripts/` module.

## Generated local files

- `output/cleanup_audit/script_dependencies_v1/summary.md`: human-readable classification and import-dependency list.
- `output/cleanup_audit/script_dependencies_v1/inventory.csv`: every local `scripts/*.py` module (except `__init__.py`), classification, source hash and reference counts.
- `output/cleanup_audit/script_dependencies_v1/references.csv`: evidence by source path, line number and snippet.
- `output/cleanup_audit/script_dependencies_v1/audit.json`: full machine-readable metadata including skipped files and parse/read errors.

## Classification meaning

- `RETAIN_DECLARED_WORKFLOW`: a documented active or currently required numerical/reporting entrypoint. This is a conservative preservation policy, not a guarantee that the script's numerical results are physically certified.
- `RETAIN_STATIC_CODE_OR_TEST_REFERENCE`: source imports the script or it is referenced in Python/tests. A test referencing an old stage runner counts as a dependency until separately reviewed.
- `REVIEW_DOCUMENTED_ENTRYPOINT`: referenced in documentation, shell scripts, or configuration; check whether the invocation remains relevant.
- `REVIEW_FROZEN_PROVENANCE`: appears in saved numerical output; retain the reproducibility reference unless an external archive and durable path are arranged.
- `REVIEW_NO_STATIC_REFERENCE`: no reference found in the scanned text. **This is not permission to delete**; external invocation or dynamically constructed names can be missed.

The AST parser identifies ordinary `import scripts.X`, `from scripts import X`, `from scripts.X import ...`, relative script imports, and literal `importlib.import_module('scripts.X')` or `runpy.run_module('scripts.X')` expressions. It also scans source, docs, shell commands, and saved outputs for literal script names. A text mention is less conclusive than an import. The scanner skips `.git`, virtual environments, cache directories, and its own generated report directory.

Input text files exceeding the default 8 MiB limit are reported as *skipped*, not assumed clean; increase the limit using `--max-text-mb N` when necessary.

## Do not retire yet

In particular, `BuildT3ReducedScalarFinalResults.py` and `BuildReducedScalarExistingInteractiveReports.py` import `RunT3GeneralScalar48Campaign.py`. The script's historical name therefore does **not** make it obsolete. The original figure/report builders, matching evaluators, threshold-running scripts and source-hash-verified fit records may have similar dependencies.

Review `references.csv` and the direct-import relationships before proposing any historical script archive. The frozen benchmark inputs and report generators must remain reproducible. Once a shortlist is reviewed against the local checkout, any archival should use an **external archive with original relative paths and SHA-256 manifests**, followed by `python -X utf8 -m pytest -q --tb=line`. This package performs no such archival.

</details>

---

<a id="source-one-time-artifact-cleanup-md"></a>

## Original document: `docs/ONE_TIME_ARTIFACT_CLEANUP.md`

**Original SHA-256:** `1f87ff36056a36c15a38ab6b514505777f1793a5dfb2109fe2290c88e9b77c3f`

**Status:** Historical record; verify current source and core documentation
before using its old commands or provisional conclusions.

<details>
<summary>Expand original: One-time artifact cleanup (local, archive-first)</summary>

# One-time artifact cleanup (local, archive-first)

This utility archives **only five explicit root-level run artifacts**, subject to local Git cleanliness, exact filename matching, simple content validation and an all-text reference scan:

- `overnight_refinement.log`
- `reduced_interactive_reports.log`
- `reduced_interactive_reports_retry.log`
- `t -q --tb=line`
- `t -q --tb=lineqqq`

The logs record historical console execution; archive them externally instead of discarding evidence. The two unusual `t -q ...` filenames were observed to contain terminal pager help at commit `fe0f511`; content changes block automatic archival.

**Not automatic deletion targets:** any `output/` or `Reports/output/` files, historical fitted configs, equations, scientific source code, tests, Wolfram exports, `README_*.txt` historical implementation notes, root `T3_*.zip` source/output bundles, existing `report_cleanup_archive/` backups, or old PowerShell launcher scripts. The audit lists several of these under `review_only` for a later explicit decision.

## Running in PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -X utf8 -m pytest -q tests/python/test_one_time_artifact_cleanup.py --tb=line
python -X utf8 -m tools.AuditOneTimeArtifacts
Get-Content -Encoding UTF8 .\output\cleanup_audit\one_time_files_v1\summary.md
```

Review the summary. Only files marked `ARCHIVE_CANDIDATE` will be handled by `--apply`. An item with a textual reference, local modification, no matching Git index entry, unusual content, or a scan error becomes `REVIEW_ONLY`. The scan is static: dynamic filename construction is not proven absent.

```powershell
python -X utf8 -m tools.AuditOneTimeArtifacts --apply --confirm ARCHIVE_ONE_TIME_ROOT_FILES
if ($LASTEXITCODE -ne 0) { throw "One-time archival failed" }
python -X utf8 -m pytest -q --tb=line
git status --short
```

On application, the script recomputes the local audit, creates a new dated external archive in the repository's parent directory, copies and SHA-256-verifies every eligible file, saves `manifest.json`, then removes the working-tree copies. Git displays the intended deletions. **No Git push occurs.** Keep the archive and manifest. For optional explicit archive location pass `--archive 'C:\thesistemp\T3_one_time_archive_custom'`; the folder must not exist already and must be outside the repo.

Do **not** automatically move or delete the README notes, root ZIPs, or source scripts merely because they are old or are only referenced in tests; they can contain historical validation and regeneration instructions. The scientific output folders are outside this tool's deletion scope.

</details>

---

<a id="source-repository-cleanup-audit-md"></a>

## Original document: `docs/REPOSITORY_CLEANUP_AUDIT.md`

**Original SHA-256:** `0526e977149314a1ca77a3d2878ccf6e2270765a8a21eb8575c895ba79cdbe60`

**Status:** Historical record; verify current source and core documentation
before using its old commands or provisional conclusions.

<details>
<summary>Expand original: T3 repository cleanup: historical-document reference audit</summary>

# T3 repository cleanup: historical-document reference audit

**Phase:** documentation audit and optional external retirement.  **Source baseline:**
`f21d15d2d38f4f313a0f2214bc52f90e1469dce2` (verified read-only on GitHub).
The installed local checkout may contain the uncommitted authoritative documents
`docs/SCALAR_SECTOR.md`, `docs/NUMERICAL.md`, and `docs/RESULTS.md`; this utility
uses the **actual local files**, not the remote snapshot, for its decisions.

## Scope and safety boundary

The audited group consists of 14 older explanatory documents from the scalar
relevance batches, zero-boundary probes, and reduced-scalar final reporting.
The intended successors are:

- `docs/SCALAR_SECTOR.md`: physical scalar potential, 50 additional quartic and
  eight quadratic operator slots, 46/4 fixed-order vertex selection, radiative
  generation, quadratic mixing and scope limitations.
- `docs/NUMERICAL.md`: numerical evolution, Yukawa fitting, report generation,
  interactive RGE and sensitivity diagnostics, convergence/provenance rules.
- `docs/RESULTS.md`: final frozen reduced-scalar benchmark outcomes, comparison
  with common A(-2) Yukawas, and historical UV-running diagnostic limitations.

**The archive is outside the repository**, by default at:

`C:\thesistemp\T3_cleanup_archive\documentation_v1\`

Nothing in `Lagrangian/`, `RGE/`, `Numerical/`, `physics/`, `scripts/`,
`tests/`, `configs/`, `output/`, `Reports/` or the fitted result directories is
retired by this tool. In particular, production and historical scripts are
not classified by filename. Some apparently historical scripts remain imports
of current reporting code; deleting them now would break reproducibility.

The audit searches current local text files, both tracked and untracked, for
historical Markdown filenames and stems. A candidate is blocked by references
from a **retained** file. References between documents archived together are
allowed. A candidate which becomes blocked can in turn block another candidate:
this transitive dependency is checked before any move.

The three authoritative master documents often list old filenames solely as
**historical provenance**. Plain mentions in those three are recorded as
non-blocking; active Markdown links or navigation instructions remain blockers.
This special rule does not apply to any other document, Python script, test,
configuration, or report. The utility excludes itself, this audit guide and
its generated report files from reference counting, and never follows symlinked
source trees. Dynamic path construction or references hidden inside binary
artifacts cannot be proved absent by text search. Review the report before
archiving even files marked safe.

## Dry-run commands (PowerShell; no unsigned .ps1 execution)

From the project root with `.venv` active:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -X utf8 -m tools.AuditDocumentationRetirement
Get-Content -Encoding UTF8 .\output\cleanup_audit\documentation_v1\summary.md
```

The audit creates only:

```text
output/cleanup_audit/documentation_v1/
    audit.json          # candidate hashes, source commit, all reference occurrences
    references.csv      # references, source filenames, line numbers, type
    summary.md          # safe/blocked/missing list and example blockers
```

The generated files are audit evidence, **not** physics inputs. A dry run moves
nothing. It does not affect the previously verified Python test baseline
(2,144 passed / 10 skipped / 1 warning in the user's checkout).

## After reviewing the report

If every safe file has been reviewed, the optional explicit command is:

```powershell
python -X utf8 -m tools.AuditDocumentationRetirement --apply --confirm ARCHIVE_SAFE_ONLY
```

This archives **only** the currently safe subset. It refuses to run if any
master document is absent or fails the SHA-256 fingerprint of the delivered
consolidated versions (with LF/CRLF normalised), if a source file could not be scanned,
if the archive destination is inside the checkout or already exists, or if a
candidate changed after the audit. For each archive member, its original SHA-256
is recorded; the copy's SHA-256 is checked before its original is removed.
The external archive writes `manifest.json` and keeps the original project
relative path (e.g. `docs/T3_SCALAR_VERTEX_BATCH2.md`).

Do **not** run `--apply` merely because the scanner reports some files safe:
first check `summary.md`, particularly whether a planned thesis source link
should remain accessible inside the repository. Files flagged blocked stay
in place. This tool does not automatically rewrite their references.

## Restoring an archived document

The archive does not modify a historical file's bytes. For a single retired
file, restore it manually from the archived path, for example:

```powershell
Copy-Item -LiteralPath 'C:\thesistemp\T3_cleanup_archive\documentation_v1\docs\T3_SCALAR_VERTEX_BATCH2.md' `
  -Destination '.\docs\T3_SCALAR_VERTEX_BATCH2.md'
```

Before restoring, ensure that destination is not already present and compare
the restored file's `Get-FileHash -Algorithm SHA256` against the manifest.
The original is retained until its verified archive copy has been written.

## Required post-archival validation

```powershell
git --no-pager status --short
python -X utf8 -m pytest -q --tb=line
```

The expected prior baseline is **2,144 passed, 10 skipped, one warning**.
Any new failure or changed skip count requires investigation. Following this
phase, use the reference audit results to update the remaining historical
links with **complete replacement files**, then audit again; do not delete
historical batch scripts or source numerical data as part of a documentation
cleanup.

**Installation prerequisite:** Both replacements from
`T3_Numerical_Results_Consolidation_f21d15d.zip` and the earlier
`docs/SCALAR_SECTOR.md` must already be installed. A dry run still reports
references if they are not, but `--apply` is disabled until all three master
documents match their delivered versions. This prevents retiring their
original sources while an old master document is still installed.

</details>

---

<a id="source-script-dependency-refinement-md"></a>

## Original document: `docs/SCRIPT_DEPENDENCY_REFINEMENT.md`

**Original SHA-256:** `29093ef82d8e70e7c52c592974ecdafe1b9100d4becbaf3bff5cc96fe5095613`

**Status:** Historical record; verify current source and core documentation
before using its old commands or provisional conclusions.

<details>
<summary>Expand original: Script dependency audit — refinement of the 77 broad-reference modules</summary>

# Script dependency audit — refinement of the 77 broad-reference modules

**Purpose:** reclassify the first-pass `RETAIN_STATIC_CODE_OR_TEST_REFERENCE` scripts without mistaking textual mentions for imports. This is a **read-only local audit**. No scripts, tests, results, PDFs, numerical fit records, or GitHub content are modified. No archival action is supplied.

## Provenance and scope

The first-pass `tools/AuditScriptDependencies.py` reported **98** script modules: 11 declared workflows, 77 with import/test/code mentions, and 10 documented-only entrypoints on the user's checkout at `fe0f511104cce8d7d914d114e3fd884eb47c6a95` (29 October 2026). The new audit rescans the **current local checkout** rather than relying on these counts or an old cached CSV. It uses the installed first-pass scanner to cover project text and independently parses all Python files to find concrete executable dependencies.

The refinement distinguishes:

1. **Syntactic Python imports:** `import scripts.X`, `from scripts import X`, `from scripts.X import ...`, and single-level relative imports inside `scripts/`.
2. **Literal dynamic imports:** literal `importlib.import_module('scripts.X')`, `__import__('scripts.X')`, or `runpy.run_module('scripts.X')`.
3. **Literal execution calls:** known `subprocess`, `os.system`, or `runpy.run_path` calls with an actual literal script/module target.
4. **Tests:** actual test imports/execution versus tests that merely mention a name in source strings or comments.
5. **Weak evidence:** Python comments/string mentions, shell wrappers, documentation, configuration, and frozen-output provenance. These are kept in separate columns; they are not automatically promoted to executable dependencies.

The closure of explicitly preserved workflows is computed by following confirmed script-to-script AST/import/execution edges, including transitive edges. Other executable dependencies and test imports remain protected by their own categories.

**Limitations:** constructed module names, reflection, subprocess commands assembled from variables, externally invoked entrypoints, and stale documentation may not be resolvable statically. `REVIEW` means **manual investigation**, *never* safe for deletion. Source metadata in historical JSON must not be edited to make a script appear unused. The earlier two full-study convenience wrappers remain potential review candidates, not automatic archive targets.

## Install (PowerShell)

Extract the ZIP at the project root, preserving exactly these new paths:

- `tools/RefineScriptDependencyAudit.py`
- `tests/python/test_t3_refined_script_dependencies.py`
- `docs/SCRIPT_DEPENDENCY_REFINEMENT.md`

The original `tools/AuditScriptDependencies.py` must remain installed. The new script imports its read-only `examine()` function; it does **not** replace or modify it.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Join-Path $HOME 'Downloads\T3_Script_Dependency_Refinement_v2.zip'
if (!(Test-Path $zip)) { throw "ZIP not found: $zip" }
$targets = @('tools/RefineScriptDependencyAudit.py', 'tests/python/test_t3_refined_script_dependencies.py', 'docs/SCRIPT_DEPENDENCY_REFINEMENT.md')
foreach ($target in $targets) { if (Test-Path $target) { throw "Already exists: $target" } }
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_t3_refined_script_dependencies.py --tb=line
python -X utf8 -m tools.RefineScriptDependencyAudit
Get-Content -Encoding UTF8 .\output\cleanup_audit\script_dependencies_refined_v2\summary.md | Select-Object -First 100
```

## Output files

Under `output/cleanup_audit/script_dependencies_refined_v2/`:

- `summary.md`: counts and named review candidates from the original 77, plus the active workflow transitive dependency closure.
- `inventory.csv`: all 98 script modules with prior and refined classes, exact hard versus weak evidence counts, and SHA-256 fingerprints.
- `strong_evidence.csv`: each resolvable import or literal execution, with referring source file, line and source excerpt.
- `audit.json`: all above information and first-pass weak references for additional review.

Run the full Python suite **only if production files are subsequently changed or scripts are archived**; no physics computation is performed by the audit.

Do not archive the 77-script group from classification totals alone. The next decision should be made from actual names in `summary.md` and the references in `strong_evidence.csv`, while preserving tested entrypoints, live imports, archived derivations, fitting provenance, and reproducibility.

</details>

---

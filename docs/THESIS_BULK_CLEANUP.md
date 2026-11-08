# T3 thesis: reversible bulk cleanup of generated LaTeX build files

This is **local-only**. The script does not access or modify GitHub, and does not
remove any Python/Wolfram/PowerShell code, scientific data, fitting output,
matching result, input benchmark, `.tex`, or `.pdf` file.

## Evidence (verified on GitHub main at 328decb8886ad3527f7cbb66667ae292f3095f4e)

- 9,794 tracked files after the eight prior local script retirements.
- 97 tracked LaTeX compiler byproducts in `Reports/output/` with a tracked
  same-stem `.tex` input: 25 `.aux`, 19 `.out`, 3 `.toc`, 25 `.log`, and
  25 `.compile.log`; approximately **665,322 bytes**. These are **generated
  build artifacts**, not the source TeX or final PDFs.
- 1,060 groups of byte-identical Git blobs, 2,250 additional file paths,
  approximately 40.46 MB of duplicate worktree bytes. These are **only
  reported**, never automatically removed: matching calculations and saved
  provenance may depend on each path.

`git` already deduplicates identical file contents in its object storage, so
removing paths is not equivalent to reclaiming their apparent worktree bytes
from `.git`.

## Install these three complete files into your LOCAL repository

```
tools/ThesisBulkCleanup.py
tests/python/test_thesis_bulk_cleanup.py
docs/THESIS_BULK_CLEANUP.md
```

## Install, preview and apply (PowerShell)

Download `T3_Thesis_Bulk_Cleanup_Round4.zip` to your Downloads folder, then:

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip = "$HOME\Downloads\T3_Thesis_Bulk_Cleanup_Round4.zip"
Expand-Archive -LiteralPath $zip -DestinationPath $repo -Force
Set-Location $repo

# Verify new tool before committing the full replacement files locally.
python -X utf8 -m pytest -q tests/python/test_thesis_bulk_cleanup.py
git add -- tools/ThesisBulkCleanup.py tests/python/test_thesis_bulk_cleanup.py docs/THESIS_BULK_CLEANUP.md
git diff --cached --check
git commit -m 'Add reversible TeX build-artifact cleanup and tests'

git status --short
python -X utf8 tools/ThesisBulkCleanup.py preview
python -X utf8 tools/ThesisBulkCleanup.py duplicates

$external = 'C:\thesistemp\T3_external_archives\thesis_bulk_cleanup'
Import-Csv "$external\preview\candidates.csv" |
    Group-Object { [IO.Path]::GetExtension($_.path) } |
    Select-Object Name, Count

# Archive only the compiler files in the exact, freshly generated plan.
python -X utf8 tools/ThesisBulkCleanup.py apply `
    --plan "$external\preview\plan.json" `
    --confirm ARCHIVE_TEX_BUILD_ARTIFACTS

.\scripts\RunCoreTests.ps1
python -X utf8 -m pytest -q tests/python/test_thesis_bulk_cleanup.py
git status --short
git diff --check
```

`preview` and `apply` require a clean **tracked** working tree. The initial
local commit ensures this. No push to GitHub is performed.

After the archival run, check `git diff --stat`, core tests and audit. If
satisfactory, locally commit the generated-artifact removals:

```powershell
git add -u -- Reports/output
git diff --cached --stat
git diff --cached --check
git commit -m 'Archive generated LaTeX compiler byproducts'
git status --short
python -X utf8 tools/ThesisRepositoryAudit.py
```

## Restore

Each archive is under `C:\thesistemp\T3_external_archives\thesis_bulk_cleanup\tex_build_<timestamp>`.
It contains `manifest.json` with the git SHA, 97+/- exact paths and SHA-256
checksums, and `files/` retaining the project-relative directory structure.

```powershell
python -X utf8 tools/ThesisBulkCleanup.py restore `
    --archive 'C:\thesistemp\T3_external_archives\thesis_bulk_cleanup\tex_build_<timestamp>'
```

The restore validates all archive copies and refuses overwriting any existing
file. A standard LaTeX rebuild also regenerates its compiler byproducts.

## Scientific preservation

The important EFT chain remains unchanged:

**UV → integrate out every heavy fermion → scalar-only EFT → integrate out
all heavy scalars → SMEFT.**

The 177 documentation-only/standalone scientific scripts and 2,250 identical-
content path candidates are **not deletion-certified**. In particular,
identical `*.wxf`, `*.json` and run-specific outputs may have distinct
provenance or be loaded by dynamically constructed paths. An explicit
reference/manifest audit is required before moving any of them.

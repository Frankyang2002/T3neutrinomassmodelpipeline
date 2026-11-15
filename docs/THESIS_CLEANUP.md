# Local thesis-first repository curation

This workflow reorganises the *working tree*, not the GitHub repository. A copy of each removed tracked file is retained outside the checkout with its SHA-256 and original relative path. It does not remove scientific code or numerical outputs automatically.

**Source audited:** GitHub `main` at `91f8d7b0f81417df751632bacb436156a04872de` and the uploaded 109-page November 2026 thesis. The local checkout may have extra files or newer changes, so every operation is recomputed from the installed checkout.

## 1. Install files, confirm checkout, run tests

Extract the accompanying ZIP into the **local** `T3neutrinomassmodelpipeline` folder, preserving paths. Do not unpack it in GitHub or use GitHub write operations.

```powershell
Set-Location C:\path\to\T3neutrinomassmodelpipeline
git status --short
git rev-parse HEAD
python -X utf8 -m pytest -q tests/python/test_thesis_repository_curation.py
```

Commit or back up existing tracked changes before applying archival. A dirty tracked checkout is rejected; untracked newly installed scripts are not deleted by the curation tool.

## 2. Read-only inventory

```powershell
python -X utf8 scripts/maintenance/CurateThesisRepository.py plan
```

The program prints an absolute `PLAN:` path and writes `plan.json` and `inventory.csv` outside the checkout, by default in its sibling directory `T3_thesis_external_archive/plans/<HEAD>/`.

The classification has four meanings:

| Decision | Meaning |
| --- | --- |
| `KEEP` | Protected thesis code, figures, benchmark paths or core scientific documentation |
| `READY_ARCHIVE` | Historical notes or TeX build auxiliaries with **no incoming literal references from retained tracked files** |
| `BLOCKED_REFERENCE` | Historical-looking material that is referenced by source, data or retained documentation |
| `REVIEW` | Scientific sources or outputs whose thesis relevance cannot safely be inferred from filenames |

All 16 models' source files, scalar EFT machinery, accepted fits, A(-4)/D(-2) comparison data, running-stage and sensitivity figures and provenance are preserved.

## 3. Apply only automatically safe historical-document cleanup

Insert the actual `PLAN:` location printed above:

```powershell
python -X utf8 scripts/maintenance/CurateThesisRepository.py apply `
  --plan "C:\path\to\T3_thesis_external_archive\plans\<HEAD>\plan.json" `
  --confirm ARCHIVE_THESIS_NONESSENTIAL
```

This operation copies and verifies all selected files *before* removing them from the checkout. It does not `git add`, `git commit` or contact GitHub. It prints a `manifest.json` path for restoration.

## 4. Remove additional material *only after source-level review*

The thesis is broader than individual file imports. A standalone plotting script or symbolic generator may have produced a saved thesis figure without being imported by the main pipeline. Therefore `REVIEW` items are never removed by the default command.

To archive a *confirmed non-thesis* item, create a CSV outside the checkout containing **specific audited paths**, one per line:

```csv
path,decision,justification
output/EXAMPLE_UNUSED_DATA.json,ARCHIVE,Verified not referenced by thesis figures or downstream saved calculations
```

The example path is illustrative; do not use it without reviewing a real path from the inventory. Run with a **different** confirmation token:

```powershell
python -X utf8 scripts/maintenance/CurateThesisRepository.py apply `
  --plan "C:\path\to\T3_thesis_external_archive\plans\<HEAD>\plan.json" `
  --review-file "C:\path\to\reviewed_non_thesis.csv" `
  --confirm ARCHIVE_MANUALLY_REVIEWED_THESIS_FILES
```

The utility checks that every manually selected item is an audited `REVIEW` item and refuses protected thesis paths, unknown paths or items referred to literally by retained tracked files. **This is a dependency guard, not proof that a scientific artifact is unused.** Manual provenance review remains mandatory.

## 5. Validate, then commit locally

```powershell
git status --short
python -X utf8 pipeline.py --help
python -X utf8 full16_symbolic.py --help
.\scripts\RunCoreTests.ps1
```

If the checks pass and the deletion inventory matches the intended thesis scope, stage/commit using your usual **local** Git workflow. Do not push as part of this cleanup task.

## 6. Restore from the external archive

Use the full `manifest.json` path printed by the operation:

```powershell
python -X utf8 scripts/maintenance/CurateThesisRepository.py restore `
  --manifest "C:\path\to\T3_thesis_external_archive\archived\<run>\manifest.json"
```

Restoration verifies the stored hashes and refuses to overwrite an existing file. Never delete the external archive before the thesis results and figures have been independently reproduced from the retained inputs.

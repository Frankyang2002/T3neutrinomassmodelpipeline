# T3 repository consolidation — local only

The project has already passed a one-point ordinary T3-B end-to-end smoke test. This cleanup is deliberately conservative: it must not break production input paths.

## Contents and installation

- `scripts/maintenance/CleanRepository.py` — full standalone cleanup script (new file).
- `CLEANUP_README.md` — this installation guide.

Extract to `C:\thesistemp\T3neutrinomassmodelpipeline\`. Nothing is pushed to GitHub.

## Run

From the repository root in PowerShell:

```powershell
$env:PYTHONUTF8 = "1"
python scripts/maintenance/CleanRepository.py
```

Review `cleanup_preview.json` before applying. It lists *every* planned move and copied file; it also lists probes held because other files reference them.

```powershell
python scripts/maintenance/CleanRepository.py --apply
```

The archive is created in a **sibling directory**, for example `C:\thesistemp\T3_repository_archives\cleanup_YYYYMMDD_HHMMSS_UTC\`, keeping your repository root clean. A SHA-256 manifest records every archived file.

## What changes

- Copies `output/` and `Reports/output/` to the external archive, **keeps the working copies**, because existing configs and staged reports can require them at runtime. No report is destroyed.
- Moves root-level historical README/stage notes into the archive. Generates `docs/PROJECT_GUIDE.md` that indexes their sections and describes the implemented path and limits.
- Moves timestamped `.backup-*` files, temporary top-level logs, and unreferenced probe scripts into the external archive. Referenced probes are retained and listed in the preview.
- Leaves production physics, configs, tests and generated RGE data intact. Does not remove the whole `Numerical/experimental/` directory, because real code and tests import it.
- Does not rewrite all source references to words like "iteration" or "stage"; they are historical identifiers in executable APIs, filenames and configs and should not be mass-replaced.

## Smoke-test again

```powershell
python scripts/maintenance/SmokeTestFullPipeline.py
```

Check `output/audit/end_to_end_smoke/smoke_result.json` for `PASS_END_TO_END_OPERATIONAL`. If it stops passing, restore and inspect dependency paths before any further cleanup.

## Restore moved files

```powershell
python scripts/maintenance/CleanRepository.py --restore "C:\thesistemp\T3_repository_archives\cleanup_YYYYMMDD_HHMMSS_UTC"
```

Restore refuses to overwrite existing files, and verifies source archive checksums. Files copied (not moved) are already present. The previous `docs/PROJECT_GUIDE.md`, if any, is backed up as `preexisting_PROJECT_GUIDE.md` in the archive.

**Limitation:** The cleaner has not run on your local checkout. This is a safe, verifiable first cleanup, not a claim that every legacy dependency or redundant research script is removed.

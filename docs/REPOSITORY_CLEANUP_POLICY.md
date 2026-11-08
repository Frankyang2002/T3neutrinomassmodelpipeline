# Repository maintenance policy

## What stays authoritative

- `README.md` and `docs/MODEL_AND_PIPELINE.md` document the **general** symmetry-allowed model and the intended fermion-first EFT ordering.
- `Lagrangian/`, `RGE/`, `Numerical/`, `physics/`, `common/`, `model/`, `configs/`, `tests/`, and runnable generation scripts remain until a source and dependency review demonstrates they are safe to retire.
- A generated JSON, CSV, Mathematica/Wolfram source, Matchete/RGBeta export, saved fitted input, manifest or threshold coefficient may be a **live numerical dependency** despite living in `output/` or `Reports/output/`. Do not delete such files solely because of their directory.

## What can be archived

- Historical root README/INSTALL files, batch notes, outdated project guides and old technical Markdown records **if no retained source/config/test directly names them**.
- Temporary backup copies, redundant ZIP exports and stray root-level scratch files after a reference check.
- Presentation-only copies (`png`, `pdf`, `html`, `svg`) under generated output directories, via an **additional opt-in flag**, after a source/config/test reference check. This is a heuristic, not a proof of runtime non-dependence.

The archive lives in a dated sibling directory outside the project. Every moved file is preserved with its original relative path and SHA-256; restoration refuses to overwrite conflicting files. Git history is secondary provenance, not a substitute for an external byte-preserving archive during the operation.

## Use

PowerShell from the **local** repository root, after extracting the package files:

```powershell
$env:PYTHONUTF8 = '1'
python -X utf8 -m pytest -q tests/python/test_extreme_repository_cleanup.py
python -X utf8 scripts/maintenance/ExtremeRepositoryCleanup.py
Get-Content -Encoding UTF8 output/cleanup_audit/extreme_repo_plan.json
```

The default is **inventory-only**. Examine `held_for_references` and `general_source_flags` in the JSON before applying. For the reviewable documentation and backup-file archive:

```powershell
python -X utf8 scripts/maintenance/ExtremeRepositoryCleanup.py --apply --confirm ARCHIVE_REVIEWED_FILES
```

For the much larger optional presentation-only artifact archive, review a fresh plan separately before invoking:

```powershell
python -X utf8 scripts/maintenance/ExtremeRepositoryCleanup.py --include-presentations
python -X utf8 scripts/maintenance/ExtremeRepositoryCleanup.py --include-presentations --apply --confirm ARCHIVE_REVIEWED_FILES
```

Then check:

```powershell
python -X utf8 -m pytest -q --tb=line
git status --short
git diff --check
python -X utf8 -m scripts.RunT3GeneralSymbolicStudy --self-test
```

A green Python test suite **does not certify** a regenerated full Matchete/RGBeta/general-numerical calculation. External generator reproducibility must be tested in the configured Wolfram environment.

To restore a specific archive created by this script:

```powershell
python -X utf8 scripts/maintenance/ExtremeRepositoryCleanup.py --restore 'C:\path\to\T3_external_archives\cleanup_YYYYMMDDTHHMMSSZ'
```

## Code retention and general-theory migration

This tool intentionally **does not delete executable `.py`/`.wl` code** based on names such as `legacy`, `probe`, `experimental` or `stage`. Python imports, Mathematica `Get[]`, saved scripts, CLI usage and dynamic filename construction can depend on it. The `general_source_flags` report identifies references to reduced or opt-in paths for an independent source-code migration.

The required migration is to make the complete operator basis the model definition, then verify downstream coupling tensors and EFT matching against that basis. Re-labeling a restricted run as general or flipping all environment flags without a mass-basis and RGE check is explicitly excluded.

No tool in this package writes to GitHub. It also does not amend scientific result files.

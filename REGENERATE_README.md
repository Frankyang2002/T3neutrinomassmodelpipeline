# T3 generated-output archive and rebuild

This operation is **local only**; it does not touch GitHub.

## Install
Extract the ZIP at `C:\thesistemp\T3neutrinomassmodelpipeline`.

New project file: `scripts/maintenance/ArchiveGeneratedResults.py`.
New project file: `REGENERATE_README.md`.

## Preview

```powershell
python scripts/maintenance/ArchiveGeneratedResults.py
Get-Content .\generated_archive_preview.json
```

## Archive and empty generated-output directories

```powershell
python scripts/maintenance/ArchiveGeneratedResults.py --apply
```

Old `output/` and `Reports/output/` are moved, not deleted, into a dated
`C:\thesistemp\T3_repository_archives\generated_outputs_*` directory.
The archive includes `manifest.json`. Source code, fitted configurations,
Matchete/RGBeta generated source outside these output folders, and tests remain.

**Warning:** The previously successful smoke test used existing generated RGE
and matching JSON. Moving outputs makes those files unavailable. The smoke test
cannot be expected to pass until you regenerate all required artifacts or restore
the archived files. This is a deliberate fresh-build test, not a regression.

## Rebuild

Inspect `python pipeline.py --help` and the repository's production README.
Regenerate the symbolic Matchete/RGBeta outputs with the installed Mathematica
environment and the appropriate production study/model selection, then execute
the production numerical pipeline and smoke test. Do not re-use old matching
data silently when asserting a fresh regeneration.

Regenerating **everything** is not a single verified command: some projects
require external Wolfram/Matchete/RGBeta installations and model-specific
settings. This package does *not* claim those steps are automated.

## Restore

If a rebuild stalls, restore to empty output directories:

```powershell
python scripts/maintenance/ArchiveGeneratedResults.py --restore "C:\thesistemp\T3_repository_archives\generated_outputs_YYYYMMDD_HHMMSS_UTC"
```

Replace the example with the actual archive path printed by `--apply`.
Restore refuses to overwrite nonempty directories, including new output files.
Move newly generated files elsewhere before restoring.

**Important:** If any files under `output/` or `Reports/output/` are tracked
by Git, moving them will show deletions in `git status`. Do not commit those
deletions until you explicitly decide how to handle generated artifacts.

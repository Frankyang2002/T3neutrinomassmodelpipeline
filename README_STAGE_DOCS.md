# Legacy notes retirement

Place `tools/RetireLegacyNotes.py` under the repository's `tools/` directory.

The utility targets 24 explicitly listed legacy documentation files, **not** `README.md`, `Reports/thesis_draft`, source code or numerical outputs. The earlier consolidation is a first pass and is **not proven lossless**. Review the new four docs before `--apply`.

From the repository root, use:

```powershell
python -m tools.RetireLegacyNotes
Import-Csv .\legacy_notes_retirement_preview.csv | Group-Object status | Select-Object Name,Count
python -m tools.RetireLegacyNotes --apply
# To undo:
python -m tools.RetireLegacyNotes --restore
```

On `--apply`, the script checks for all four consolidated documents, records SHA-256 checksums, and moves notes into `.cleanup_quarantine/legacy_notes/` with `manifest.json`. The script does not modify Git's index, commit, or push.

WARNING: Existing markdown links referring to these files may break. After relocation, review source-code references and incoming links before committing. Do not run `git add -A` blindly on an existing set of staged cleanup changes.

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

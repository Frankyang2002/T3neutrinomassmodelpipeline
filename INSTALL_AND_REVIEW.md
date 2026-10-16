# Documentation consolidation review

This ZIP contains complete replacement files at these exact project paths:

- `README.md`
- `docs/PHYSICS.md`
- `docs/NUMERICAL.md`
- `docs/RESULTS.md`

It is a **first-pass concise consolidation** of the public repository documentation, not an independently validated physics update or a line-by-line lossless replacement of every existing note.

## Install locally (PowerShell)

```powershell
# From C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_Documentation_Consolidation.zip" -DestinationPath . -Force
Get-Content .\README.md -TotalCount 15
Get-ChildItem .\docs\*.md
```

Move the downloaded ZIP to another location or adjust the path if necessary. Existing `INFO_*.md`, historical README files, and archived notes are **not** deleted.

## What remains to review before removing old notes

1. Ensure specialist file-by-file notes (`INFO_IMPLEMENTATION.md`, `INFO_RGE.md`, `INFO_TESTS.md`) contain no still-important details missing from the concise new documents.
2. Confirm recent local amendments to `Reports/notes/T3_Numerical_Findings.md` if local files differ from the GitHub `main` version.
3. Recheck paths and CLI flags against the actual checked-out code; some old README files describe replacement packages rather than production commands.
4. Do not untrack or delete historical notes until the technical content and citation references have been reviewed.
5. Keep documentation work separate from the large already-staged Git deletion cleanup.

A documentation-only check:

```powershell
git diff -- README.md
git status --short
```

No GitHub write operations were performed to produce this package.

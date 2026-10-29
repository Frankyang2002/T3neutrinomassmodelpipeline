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

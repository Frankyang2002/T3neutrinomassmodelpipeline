# Stage 12 — post-quarantine audit and narrow Git hygiene

Project paths:
- `.gitignore` (full replacement)
- `tools/CheckReportCleanup.py` (new, read-only tool)

Commands from repository root:

```powershell
python -m py_compile tools/CheckReportCleanup.py
python -m tools.CheckReportCleanup
Import-Csv .\report_cleanup_check.csv | Group-Object kind,detail | Select-Object Name,Count
# Inspect Git's index; tracked generated files remain tracked even if ignored:
git status --short
git ls-files Reports/output | Measure-Object -Line
```

The audit checks which Reports/output files are tracked, whether they still exist, and some literal references to missing TeX. It cannot prove that dynamically constructed paths are unused.

`.gitignore` excludes only generated displays, figure-review binary copies, and local cleanup material. It intentionally does NOT ignore Reports/output/full because some scripts still read/write numerical scan data there. Ignore rules do not untrack files already committed; do not run a mass `git rm --cached` until the numerical-data migration and reference audit are complete.

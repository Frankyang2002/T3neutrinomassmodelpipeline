# Stage 10 — Exact duplicate TeX audit

Install `tools/AuditDuplicateTex.py` in the project root, preserving directories.

PowerShell (run in repository root):

```powershell
python -m tools.AuditDuplicateTex
Import-Csv .\duplicate_tex_audit.csv | Group-Object classification | Select-Object Name, Count
Import-Csv .\duplicate_tex_audit.csv | Where-Object classification -eq 'EXACT_DUPLICATE' | Select-Object -First 20 path, canonical_path
```

This is **read-only**. SHA-256 exact matches are candidate redundant files, not automatic deletion approvals. Keep `.tex` sources until generation/reference dependencies are assessed.

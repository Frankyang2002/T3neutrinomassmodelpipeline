# Stage 8 — TeX/PDF report audit

Project path: `tools/AuditTexPdfReports.py`

Run from the repository root:

```powershell
python -m tools.AuditTexPdfReports
Import-Csv .\tex_pdf_report_audit.csv | Group-Object classification | Select-Object Name, Count
```

Scope only `comparison` if desired:

```powershell
python -m tools.AuditTexPdfReports --root Reports/output/full/comparison --manifest comparison_tex_pdf_audit.csv
```

No deletion occurs. `.tex` sources are explicitly retained until their purpose is verified. PDF/TeX same-stem pairs are flagged, not assumed to be duplicates.

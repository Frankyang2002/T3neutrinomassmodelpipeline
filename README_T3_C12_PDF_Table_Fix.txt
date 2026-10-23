Purpose
-------
Reformat the existing analytical GroupFactors PDF so the C12 information is
shown as comparative run-block tables, matching the style already used for C5.

Files
-----
- Reports/C12GroupFactorPDF.py
- scripts/BuildCombinedAnalyticalReports.py
- tests/test_c12_group_factor_pdf_tables.py

Install
-------
Copy these files into the repository root, preserving paths.

PowerShell commands
-------------------
python -m pytest -q tests/test_c12_group_factor_pdf_tables.py

# Rebuild the existing analytical PDFs from saved runs
python scripts/BuildCombinedAnalyticalReports.py --check
python scripts/BuildCombinedAnalyticalReports.py

# If you only want to rewrite/recompile the existing C12 PDF table without
# rebuilding all analytical reports:
python -m Reports.C12GroupFactorPDF `
    --tex "Reports\output\analytical\GroupFactors\GF_EFT_1_after_F.tex" `
    --compile

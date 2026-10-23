C12 GROUP FACTORS PDF — REPAIR OF THE MISSING-SECTION ERROR
==========================================================

Bug fixed
---------
Previous Reports/C12GroupFactorPDF.py expected a pre-existing per-model C12
section in the TeX file. That section is not always present, so the command
raised ReportRewriteError. The replacement reads the saved Matchete/CG and
one-loop tensor-RGE data directly and builds C5-style model comparison tables
without parsing a previous C12 section.

This changes report presentation only; fixed-order and coupled numerical
running, operator matching, and the existing C5 group-factor tables are not
changed.

Install
-------
Copy every file in this ZIP into the repository root, preserving paths:

- Reports/C12GroupFactorPDF.py                  FULL REPLACEMENT
- Reports/PipelineReports.py                    FULL REPLACEMENT
- scripts/BuildCombinedAnalyticalReports.py     FULL REPLACEMENT
- RGE/group_factors/validation/C12RepresentationGroupFactors.py FULL REPLACEMENT
- tests/test_c12_group_factor_pdf.py            FULL REPLACEMENT
- tests/test_c12_representation_group_factors.py FULL REPLACEMENT
- tests/test_c12_group_factor_block_tables.py   NEW

From C:\thesistemp\T3neutrinomassmodelpipeline in PowerShell:

    python -m pytest -q tests/test_c12_group_factor_pdf.py tests/test_c12_representation_group_factors.py tests/test_c12_group_factor_block_tables.py

QUICK UPDATE: existing analytical TeX/PDF using saved models (no physics rerun):

    python -m Reports.C12GroupFactorPDF `
        --tex "Reports\output\analytical\GroupFactors\GF_EFT_1_after_F.tex" `
        --source "output\full\comparison\smallY_smallL" `
        --compile

    Start-Process "Reports\output\analytical\GroupFactors\GF_EFT_1_after_F.pdf"

FULL ANALYTICAL REPORT REBUILD: existing 16-model study from saved data:

    python scripts/BuildCombinedAnalyticalReports.py --check
    python scripts/BuildCombinedAnalyticalReports.py

If the saved run data are under another study, use --source to set its folder.
For a single representation, specify --data-dir and --label:

    python -m Reports.C12GroupFactorPDF `
        --tex "Reports\output\analytical\GroupFactors\GF_EFT_1_after_F.tex" `
        --data-dir "output\T3_B_alpha_m1\data" `
        --label "T3-B-m1" --compile

Group factors
-------------
Each structure row is extracted from the saved *total* C12 beta after
projecting the y1* y2* matching channel, with separate model columns.
If the total beta is not collinear on the matched C12 seed support, only
separately aligned source contributions are printed and an explicit
beta_perp != 0 row and caveat are added. Missing tensor inputs are never
silently replaced with guessed coefficients.

The C5 / direct-Weinberg group-factor tables already in the existing
GroupFactors PDF remain unchanged.

Local validation
----------------
21 focused Python tests passed. A smoke-generated TeX report containing
all 16 presentation model columns compiled successfully with pdflatex,
without overfull lines. This used synthetic sample coefficients for layout
validation, not the user's saved 16-model numerical output.

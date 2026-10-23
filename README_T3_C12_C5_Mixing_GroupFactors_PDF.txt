T3 GroupFactors PDF: C12 non-closure and explicit C12 <-> C5 mixing
================================================================

This is a reporting-only replacement. No source physics, RGEs, or numerical
pipeline is changed.

FULL REPLACEMENT FILE
---------------------
Reports/C12GroupFactorPDF.py

NEW TEST FILE
-------------
tests/test_c12_c5_mixing_pdf.py

Behaviour
---------
Updates the EXISTING TeX/PDF report:
Reports/output/analytical/GroupFactors/GF_EFT_1_after_F.tex
Reports/output/analytical/GroupFactors/GF_EFT_1_after_F.pdf

Each C12 run block now includes separate rows for:
  - lambda5 C12 -> beta_C5
  - conjugate(lambda5) C5 -> beta_C12

The forward factor is reused from the existing analytic R_W calculation
(RGE/group_factors/validation/ValidateDirectWeinbergGroupFactors.py).
Its independent Wigner/CG phase-normalisation validation remains open.

The reverse factor has been checked for T3-B, alpha=-1 only:
  8*sqrt(3)/3 (ordered-tensor two-operator convention).
For the remaining representations, the PDF prints 'not verified'.

The previous beta_perp != 0 row has been removed from the group-factor table.
Any non-collinearity is instead described as an operator-closure warning
BELOW the table. This quantity is not C12 <-> C5 mixing.

The previous SM+Weinberg C5 report tables remain unchanged.

Install
-------
Extract files into C:\thesistemp\T3neutrinomassmodelpipeline and overwrite
Reports/C12GroupFactorPDF.py when asked.

PowerShell from the project root:

python -m pytest -q tests/test_c12_c5_mixing_pdf.py tests/test_c12_group_factor_pdf.py tests/test_c12_group_factor_block_tables.py

python -m Reports.C12GroupFactorPDF `
    --tex "Reports\output\analytical\GroupFactors\GF_EFT_1_after_F.tex" `
    --source "output\full\comparison\smallY_smallL" `
    --compile

Start-Process "Reports\output\analytical\GroupFactors\GF_EFT_1_after_F.pdf"

Alternatively, your normal combined analytical report builder invokes this
same PDF-updating module after rebuilding the existing report.

Verification
------------
25 local unit/regression tests passed (including old projection/reporter
regressions). The TeX compiled in a representative smoke test. The local
environment does not contain your complete fitted benchmark repository, so
rerun the command on your Windows checkout to produce the actual PDF.

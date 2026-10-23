T3 C12 REPRESENTATION-SPECIFIC GROUP FACTORS IN EXISTING PDF REPORTS
====================================================================

ROOT: C:\thesistemp\T3neutrinomassmodelpipeline

Full replacement / new files (extract the ZIP over the root):
  RGE/group_factors/validation/C12RepresentationGroupFactors.py   [new; reused from previous HTML bundle]
  Reports/C12GroupFactorPDF.py                                    [new]
  Reports/PipelineReports.py                                      [complete replacement]
  scripts/BuildCombinedAnalyticalReports.py                       [complete replacement]
  tests/test_c12_group_factor_pdf.py                              [new]

No changes to Wilson tensor RGEs, numerical solvers, or physical matching.
No HTML files are created or edited by this package.

WHAT IT DOES
------------
- Projects the saved one-loop Wilson-tensor beta onto the matched C12 tensor,
  using the same real-scalar/SU(2) Clebsch-Gordan seed as the existing RGE.
- Isolates the mixed conjugate(y1)*conjugate(y2) channel and checks exact
  seed-support proportionality separately for each group-factor source.
- Inserts one model-labelled C12 section into the SAME existing GroupFactors
  TeX, then recompiles the SAME named PDF (not a side report).
- Reports non-collinear terms explicitly, and distinguishes seed-support
  collinearity from complete operator closure.
- Supports multiple different representations in the across-model report;
  each model uses its own saved tensors and RGBeta representation metadata.
- T3-B alpha=-1 also displays the previously verified paper-basis coefficient
  if the exact source-derived regression succeeds.
- Hooks both ordinary Reports/PipelineReports.py and the standalone combined
  analytical report generator, scripts/BuildCombinedAnalyticalReports.py.

REQUIRED SAVED INPUT FILES IN EACH MODEL'S data/ DIRECTORY
-----------------------------------------------------------
  eft1_after_F_wilson_seed.json
  eft1_rgbeta_rge.json
  eft1_wilson_rge.json

Missing files are explicitly labelled unavailable in the report, not replaced
by inferred coefficients. A tensor proportionality failure is also reported.
The saved tensor RGE uses one SM lepton generation; full three-generation
operator-space closure has not been established.

TESTS (PowerShell, from project root)
-------------------------------------
python -m pytest -q tests/test_c12_group_factor_pdf.py

UPDATE THE EXISTING T3-B REPORT NOW; NO FULL PHYSICS RERUN
----------------------------------------------------------
$Tex = "Reports\output\analytical\GroupFactors\GF_EFT_1_after_F.tex"
$Data = "output\single\T3_B_alpha_m1\data"
if (!(Test-Path (Join-Path $Data "eft1_wilson_rge.json"))) {
    $Data = "output\T3_B_alpha_m1\data"
}
python -m Reports.C12GroupFactorPDF `
    --tex $Tex `
    --data-dir $Data `
    --label "T3-B alpha=-1"
Start-Process ($Tex -replace '\.tex$', '.pdf')

This adds a single-model C12 section to the selected existing report.
If GF_EFT_1_after_F.tex is a 16-model comparison, update ALL 16 together
via the combined report generator below to avoid showing only T3-B.

REBUILD THE COMBINED ANALYTICAL REPORT WITH ALL MODELS
-------------------------------------------------------
python scripts/BuildCombinedAnalyticalReports.py --check
python scripts/BuildCombinedAnalyticalReports.py
Start-Process "Reports\output\analytical\GroupFactors\GF_EFT_1_after_F.pdf"

The combined script requires the same 16 saved analytical model runs as before.
It regenerates the existing analytical report set and then replaces only the
C12 section of the F-first GroupFactors TeX/PDF. No Matchete or RGBeta reruns.

NORMAL FULL PIPELINE
--------------------
Reports/PipelineReports.py now also injects this section when the normal
stage-aware GroupFactors PDF report generator runs. Existing PDF names remain
unchanged. Standalone numerical opt-in benchmarks do not themselves regenerate
the full analytical PDFs.

NOTE ON LATEX
-------------
Requires pdflatex or latexmk, as for all existing GroupFactors PDFs.
If no PDF compiler is installed, the TeX will be updated, but the CLI raises
an error rather than claiming the PDF has been rebuilt.

T3-B coupled C12 <-> C5 RGE in GroupFactors HTML reports
=========================================================

Scope: T3-B alpha=-1 (dS1=dS2=2, dF=1), all F generations matched
simultaneously at mu_F. Reporting only: no RGE or physics files changed.

Install this ZIP at:
    C:\thesistemp\T3neutrinomassmodelpipeline

FULL FILES (not patches):
    Reports/GroupFactorCoupledHTML.py                   NEW
    Reports/PipelineReports.py                          REPLACE
    Numerical/orchestration/T3BIntegratedCoupledPipeline.py  REPLACE
    tests/test_t3b_group_factor_html.py                NEW
    README_T3B_GroupFactors_HTML.txt                   NEW

Test (PowerShell, project root):
    python -m pytest -q tests/test_t3b_group_factor_html.py

FASTEST: rebuild HTML companions using existing GroupFactors TeX files,
without rerunning any numerical calculation, Matchete or RGBeta:
    python -m Reports.GroupFactorCoupledHTML --report-root "Reports\output\analytical" --t3b

Open the F-first EFT HTML:
    Start-Process "Reports\output\analytical\GroupFactors\GF_EFT_1_after_F.html"

The new HTML appears beside its existing GF_EFT_1_after_F.tex/.pdf and
includes beta_C12, beta_C5, both mixing directions and group-factor table.
It is offline readable, with copyable LaTeX and links to the full PDF.
UV and final-EFT stages get their appropriate (not copied) RGE summaries.

FUTURE integrated coupled benchmark runs automatically refresh HTML
companions IF saved GroupFactors TeX files already exist:
    python -m Numerical.orchestration.T3BIntegratedCoupledPipeline `
      --config "Numerical\benchmarks\fitted_reference\T3_dS1_2_dS2_2_dF_1_alpha_m1.json" `
      --output-dir "output\single\T3_B_alpha_m1" `
      --enable-coupled --max-order 4

The standard physics pipeline's Reports/PipelineReports.py also writes HTML
companions when it generates stage-aware GroupFactors TeX/PDF reports.

IMPORTANT: The HTML expression describes the retained TWO-OPERATOR
leading-log-resummed equation, not complete higher-order EFT matching.
The displayed T3-B coefficients must not be applied to other model classes.
The order-one production fixed-order result remains authoritative.

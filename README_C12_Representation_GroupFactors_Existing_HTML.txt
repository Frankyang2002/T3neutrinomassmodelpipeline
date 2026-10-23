T3-B and representation-specific C12 group factors in EXISTING GroupFactors HTML
================================================================================

Purpose
-------
Update the report HTML already in your project (not a new HTML report).
Use saved Matchete EFT1 Wilson CG seed, RGBeta representation metadata and
saved one-loop scalar-only master-tensor beta to extract the group factors
for the mixed C12 ~ conjugate(y1)*conjugate(y2) matching channel.

Full replacement/new project paths in the ZIP
---------------------------------------------
NEW:         RGE/group_factors/validation/C12RepresentationGroupFactors.py
NEW:         Reports/C12GroupFactorHTML.py
REPLACEMENT: Reports/GroupFactorCoupledHTML.py
REPLACEMENT: Numerical/orchestration/T3BIntegratedCoupledPipeline.py
NEW:         tests/test_c12_representation_group_factors.py
This README belongs at the repository root.

The existing Reports/PipelineReports.py from the previous report HTML
integration ZIP remains unchanged. It calls GroupFactorCoupledHTML, which now
uses the source-derived factors for any F-first ordinary T3 models whose
saved tensor inputs are present. The integrated optional coupled runner
passes its own numerical data directory to the same report HTML writer.

PowerShell - run from C:\thesistemp\T3neutrinomassmodelpipeline
---------------------------------------------------------------

python -m pytest -q tests/test_c12_representation_group_factors.py tests/test_t3b_group_factor_html.py

$Html = 'Reports\output\analytical\GroupFactors\GF_EFT_1_after_F.html'
$DataDir = 'output\single\T3_B_alpha_m1\data'
if (!(Test-Path (Join-Path $DataDir 'eft1_wilson_rge.json'))) {
    $DataDir = 'output\T3_B_alpha_m1\data'
}
Test-Path $Html
Test-Path (Join-Path $DataDir 'eft1_wilson_rge.json')

python -m Reports.C12GroupFactorHTML `
    --html $Html `
    --data-dir $DataDir `
    --label 'T3-B alpha=-1'

Start-Process $Html

If your EXISTING HTML is in a different report directory, set $Html to
that exact existing file. The program intentionally refuses to create a
new report if --html does not already exist.

Or regenerate all saved GroupFactors HTML companions, preserving any existing
HTML content (exact source-derived factors if --data-dir is supplied):

python -m Reports.GroupFactorCoupledHTML `
    --report-root Reports\output\analytical `
    --t3b `
    --data-dir $DataDir `
    --preserve-existing-html

The user's normal integrated opt-in command continues to work:

python -m Numerical.orchestration.T3BIntegratedCoupledPipeline `
    --config 'Numerical\benchmarks\fitted_reference\T3_dS1_2_dS2_2_dF_1_alpha_m1.json' `
    --output-dir 'output\single\T3_B_alpha_m1' `
    --enable-coupled --max-order 4

What is source-derived, and what remains unverified
---------------------------------------------------
- The Matchete tree-level CG contraction fixes each C12 tensor component.
- The saved master Wilson-tensor beta provides the one-loop component betas.
- Mixed y1* y2* terms are isolated exactly and divided by their matched
  tensor components, with a collinearity check on ALL C12 seed-support components.
- Individual beta sectors are shown separately; a term is marked
  non-collinear if no single coefficient can be extracted.
- For T3-B alpha=-1 ONLY, the verified lambda12Inv1/Inv2 to
  lambda12/lambda12Adj mapping is applied and compared to
  -3*g2^2 + 2*lambda12 - lambda12Adj/2 + ye*conjugate(ye).
- For other F-first T3 representations, the report shows actual
  source-derived factors in their saved RGBeta tensor basis, without copying
  T3-B coefficients or inventing a quartic-basis conversion.
- Output records the number of components outside the C12 matching support.
  This calculation is NOT a complete three-generation operator-closure proof.
- Neither physical Wilson evolution nor standard production matching is
  changed by this presentation-only extension.

Validation in the available container: 21 local tests passed (including
existing HTML writer regression tests). The real saved T3-B tensor seed and
full-benchmark reporting are to be checked by running the commands locally.

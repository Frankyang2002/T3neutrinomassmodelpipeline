THESIS BENCHMARK + COMPUTED F-HARD TRIANGLE (LOCAL ONLY)
=========================================================

Intended project path of the new full file:
  scripts/RunT3ThesisHardTriangle.py

This package is a NEW wrapper. No existing source file is modified.

Verified against read-only public GitHub source structure (13 Nov 2026):
- scripts/BuildT3GeneralisedThesisTables.py generates the actual 32 thesis
  fixed-Yukawa runs (D_-2 and A_-4 donors across 16 model representations).
- Numerical/parallel_t3/complete/all16_rough_yukawa_fit.py currently sets
  hard_c5_zero_benchmark=True during Stream1 and immediately runs EFT1.
- Numerical/parallel_t3/wilson/fermion_hard_triangle_batch1.py implements
  install_hard_boundary(...), which injects the computed hard C5 into
  the fermion-scale input .npz and preserves the LLSS component.

The new wrapper calls the original thesis generator, intercepts its
fresh Stream1 boundary after it is built, injects the existing hard-triangle
computation using that trial's B background, and continues with exactly the
same subsequent RGE, scalar matching, SMEFT and PMNS steps as the original.

RUN COMMANDS (PowerShell; run from your local repository root)
-------------------------------------------------------------
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

# Extract this ZIP with its project-relative scripts/ path:
Expand-Archive -LiteralPath "$HOME\Downloads\T3_thesis_hard_triangle_local.zip" -DestinationPath . -Force

# New output directory MUST not exist. The donor fit reproduces the exact
# source Yukawas recorded by the old thesis tables' provenance.json.
.\.venv\Scripts\python.exe scripts\RunT3ThesisHardTriangle.py `
  --reference-results output\parallel_t3\A4_Dm2_targeted_fit01 `
  --output output\parallel_t3\thesis_generalised_hard_triangle_v1 `
  --wilson-solver RK2 --wilson-steps 2 --mass-bilinears zero `
  --mu-low 100 --gamma-factor 1

# After success inspect the evidence and compare with the original baseline:
$old = Import-Csv output\parallel_t3\thesis_generalised_tables_v1\all32_observables.csv
$new = Import-Csv output\parallel_t3\thesis_generalised_hard_triangle_v1\all32_observables.csv
$new | Format-Table reference,model,status,m3_eV,dm21_1e-5_eV2,dm3l_1e-3_eV2 -AutoSize
Get-Content output\parallel_t3\thesis_generalised_hard_triangle_v1\hard_triangle_injection_ledger.json
Get-Content output\parallel_t3\thesis_generalised_hard_triangle_v1\provenance.json

OUTPUTS
-------
output/parallel_t3/thesis_generalised_hard_triangle_v1/
  - results.json, provenance.json, summary.json
  - Table_B1.tex through Table_B4.tex
  - all32_observables.csv
  - hard_triangle_injection_ledger.json
  - hard_fermion/D_-2/<model>/fermion_hard_triangle.{json,npz}
  - hard_fermion/A_-4/<model>/fermion_hard_triangle.{json,npz}

The original thesis_generalised_tables_v1/ remains unchanged.

CRITICAL SCOPE
--------------
1. This injects ONE hard topology only: heavy-fermion scalar-quartic triangle.
   It does NOT implement complete renormalized one-loop fermion matching.
2. The original conditional GammaCC-to-Weyl mapping, selected/reduced scalar
   potential, and candidate anomalous dimensions are not independently
   certified by this wrapper.
3. DO NOT claim the submitted thesis includes this term until the local run
   completes 32/32 and the resulting numerical tables are substituted into
   the LaTeX thesis. Figures in sections 5.1-5.3 and Appendix B remain based
   on the old hard-zero benchmark until regenerated with matching included.
4. The fixed Yukawa donors are reused without refitting. Chi2 values can
   change; do not silently retain statements claiming all fits remain <0.1.
5. Every failure is blocking. Check the ledger and the per-run tracebacks;
   never replace an exception by a zero matching coefficient.
6. The script has been syntax checked in this environment, but the full
   local repository and its runtime dependencies/saved B inputs are not
   available here. Its 32/32 numerical completion is NOT verified.

ACCURATE THESIS CLAIM AFTER THE ABOVE CONDITIONS ARE MET
---------------------------------------------------------
"The numerical pipeline includes the finite one-loop scalar-quartic
heavy-fermion triangle contribution to the Weinberg coefficient at the
fermion threshold, in addition to the tree-level LLSS boundary, intermediate
EFT running, and simultaneous scalar-threshold matching. The direct
fermion-hard contribution remains a selected topology; complete one-loop
fermion matching, including other hard structures and counterterms, is
outside the scope of this benchmark."

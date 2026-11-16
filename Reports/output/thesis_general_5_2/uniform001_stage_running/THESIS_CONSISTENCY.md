# Comparison with `T3_Honours_Thesis (3).pdf` (106-page uploaded draft)

## What still agrees with the new 0.01 transfer tables

- Section 5.1's qualitative observation that most saved low-scale C5 flavour matrices are almost proportional across representations under fixed D(-2) or A(-4) Yukawas survives. Of the 32 transfers, only the four B(-1)/C(-1) cases have Weinberg residuals of order 0.1–0.3. The other 28 are below ~5.5e-4.
- B(-1) and C(-1) still show substantial non-proportional contributions and changed PMNS mixing angles with fixed cross-model Yukawas.
- The D(-2) reference is more transferable under the chosen approximate chi-squared cutoff: 5/16 versus 2/16 for A(-4) (including respective self-reference rows).
- The new C5 coefficient for fixed D(-2) Yukawas transferred to A(-2) is `k=0.99814924`, `epsilon=9.37308e-5`, rather than the thesis Section 5.1's approximate `k≈0.9972`. The group-factor-only expression is a schematic approximation, not an exact equality.

## Corrections before treating the PDF as an up-to-date thesis

1. **Section 4.1, PDF page 64**: The paragraph saying unlisted representation-dependent scalar couplings are zero is wrong for the new benchmark. Every *allowed BSM scalar quartic* other than lambdaT3 is +0.01 at muUV. Conjugate Yukawas are zero **at the UV boundary**, but their operators remain symmetry-allowed. Mixed scalar bilinears are dimensionful and retain their existing mass-squared prescription; SM Higgs quartic lambdaH is not 0.01.
2. **Section 4.2, PDF pages 64–65**: The new accepted reference fits keep lambdaT3 fixed (A–D 0.1; E 0.05) while refitting the 18 real entries of Y1 and Y2. Existing text describes the older 19-parameter Sobol/local hybrid. Acknowledge the original source fits as *seeds*, not accepted new-benchmark results.
3. **Section 5.1, PDF pages 68–70**: Replace the outdated single-reference numbers, notably k≈0.9972, and change the 'Appendix C / T3-A-m2 benchmark' wording to the new D(-2) and A(-4) reference transfers in Appendix B. Use 32/32 current tables. The statement '9 ... direct proportional; 4 ... almost proportional' is not defined by the new residuals; do not keep these counts without defining a threshold.
4. **Section 5.2, PDF pages 70–74**: Both numerical heatmaps (Figures 5.2, 5.3) are from the old general-fit/frozen-stage data, *not* the new uniform-0.01 reference campaign. Do not retain their exact percentages (e.g., 'at most 2%') without rerunning. After generating the new figures, interpret the new CSV instead. Distinguish the conditional one-scale calculation from independently evolved EFT1 Wilson transport.
5. **Section 5.3, PDF pages 73–75**: The new T1/T2/T3 trajectory report is on GitHub at `Reports/output/thesis_general_5_3/zee_uniform001_ovals_v1/full_sensitivity/`; ensure the thesis figure includes the new PNG and its provenance, not the old general-fit image. The main B(-1)/C(-1) mechanism remains plausible but must be checked against its new trajectory slopes when making quantitative claims.
6. **Chapter 3, PDF pages 56–60**: Its text claims sequential Matchete matching and intermediate-EFT Wilson running *are implemented*. The 0.01 numerical programme documents conditional **one-scale** full-general C5. State that the strict UV -> all fermions -> all scalars -> SMEFT sequence is the intended EFT construction, not a certified numerical prediction of the current 0.01 plots.
7. **Conclusion, PDF pages 77–78**: Qualitative conclusions about representation dependence and B(-1)/C(-1) remain consistent. Avoid saying the new running study reproduces older percentage sizes until its 16 model rows pass baseline checks.

## Graph data provenance verified on GitHub

- Older heatmaps from `Reports/output/thesis_qualitative/representation_running/csv/frozen_stage_observable_changes.csv`, validation pinned to commit `6ebe9937...`, old fitted run `thesis16_merged_v1`.
- New reference fits: `output/studies/thesis_qualitative/uniform_scalar_001/campaigns/uniform001_all16_reference_v1/reference/` (accepted 16/16).
- New D(-2)/A(-4) comparison: `Reports/output/thesis_general_5_1/uniform001_fixed_yukawa_tables/uniform001_Dm2_Am4_fixed_v1/` (32/32).
- New coupling-sensitivity ovals: `Reports/output/thesis_general_5_3/zee_uniform001_ovals_v1/full_sensitivity/` (separate from running stage heatmaps).

There are **no verified 16-model running-stage ablation outputs under the new scalar-0.01 reference** on GitHub at the time of this check. The provided script computes them locally and refuses to render thesis graphs from incomplete data.

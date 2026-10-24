# Existing GitHub comparison studies — verified thesis figure placement

**Repository inspected:** `Frankyang2002/T3neutrinomassmodelpipeline`, `main`, commit `648a3cb564bd6b0d1dc91ca355e5f542f20491f9` (24 October 2026). This analysis uses the repository's **already saved data**; no UV RGE, matching, optimisation, or EFT integration was rerun here.

## Main conclusion of the source audit

**The UV-on versus UV-frozen comparison already exists; so does the full 16-model high-scalar, UV-frozen study.** Do not rerun those expensive numerical workflows for the thesis. The results need careful reporting, since the saved UV-sensitivity studies use different fitted baselines and the exact same-seed ablation is not a successful neutrino fit.

The current source files are:

- `Numerical/analysis/CompareUVFermionSensitivity.py` and `output/studies/fermion_start_lambda_comparison/uv_sensitivity_comparison/{summary.json,summary.csv,local_slopes_by_model.csv}`. Includes all 16 models × 3 rescaling directions × 5 observables × 2 fit modes = **480 slope records**. Paired in the sense of common model and direction; **not** paired at the same numerical Yukawa fit or same gamma grid.
- `Numerical/orchestration/DiagnoseLambda5e10UV.py` and `output/studies/lambda_5e-10_refit/diagnostics/a_m2_exact_seed_uv_ablation/ablation.json`. One **same-seed** study of `000`, `100`, `110`, `111`. Baseline is T3-A(-2) at `lambdaT3=5e-10`; full-stage `chi2≈15715` is a poor fit.
- `Numerical/analysis/AnalyseUVScalingBreaking.py`, `DecomposeUVScalingBreaking.py`, and `CorrelateUVBreakingWithSensitivity.py` plus saved outputs under `output/studies/lambda_5e-9_refit/`. Instantaneous UV beta-function scaling-breaking diagnostic and correlations with sensitivity; correlations **do not establish causation**.
- `output/studies/fermion_start_lambda_comparison/{5e-9,0p5}/comparison/{optimised,shared_a_m2}/summary.json`. Both cases contain all 16 model results. The `0p5` study has `lambdaT3=0.5` for A–D, **0.25 for E**; the `5e-9` study similarly uses `2.5e-9` for E.
- `Reports/output/display_fermion_start/0p5/running/optimised/` and `/shared_a_m2/`, both with 16 saved EFT1 and Weinberg trajectories.
- `output/studies/fermion_start_lambda_comparison/0p5/eft1_perturbativity_screen.json`. All **32** sampled quartic trajectories passed the simple `max|quartic|<4pi` magnitude screen. It is *not* a full stability/unitarity/perturbativity proof.

## Which new figures to use

All files go to `Reports/output/thesis_assets/figures/` and are available in **PDF and PNG**. The numbering continues from the previously delivered `fig01`–`fig10` thesis package. Figure numbers are provisional; update `\label` references to match the LaTeX thesis.

| Thesis section | Figure | Use | What is proven / crucial caveat |
|---|---|---|---|
| **5.3 UV evolution and breaking of matching degeneracies** | `fig11_saved_uv_vs_fermion_sensitivity` | **Main text**, alongside a concise slope table | The saved UV-refit sensitivity slopes are much larger than the nearly zero fermion-start slopes. It is a **descriptive** comparison of independently fitted baseline points with different gamma spacings; do not label as a same-input UV on/off experiment. |
| **5.3 (or Appendix D)** | `fig12_same_seed_staged_ablation` | **Supporting** same-seed control | Switching on `UV`, then `EFT1`, then `Weinberg` clearly separates effects. The benchmark has **large chi-square** and is not a viable fit; avoid using it as your principal phenomenological result. |
| **5.4 High versus low scalar coupling** | `fig13_low_high_scalar_fit_quality` | **Main text** | There are 16 saved fits per case/mode. With the study's `chi2≤5.86` screen: **15/16** low-scalar optimised and **13/16** high-scalar optimised benchmarks pass. Fitting success varies; don't interpret `status=Success` as chi-square acceptance. Fitted observables are fit targets, not independent predictions. |
| **5.5 Representation dependence at shared Yukawas** | `fig14_shared_yukawa_representation_effect` | **Main text** | Within each coupling case the A(-2) Yukawa matrices are imposed for all 16 models. The predicted `Δm²_3ℓ` varies widely by representation, especially for C and E. The A(-2) donor matrices differ **between** low- and high-coupling cases, so overlap of the two curves is a consistency observation, not a same-input scalar-coupling test. |

Tables are in `Reports/output/thesis_assets/tables/`:

- `table_uv_sensitivity_maxima.csv` (supports Fig. 11),
- `table_staged_ablation.csv` (supports Fig. 12),
- `table_existing_scalar_comparison.csv` (supports Figs. 13 and 14, all 16 model values).

## How this changes the thesis plan

- **Section 5.2 — matching rescaling:** keep analytic `T1/T2/T3` derivation and the saved sensitivity scan. Prove the *leading* bilinear scaling algebraically, then show the saved UV sensitivity; distinguish initial boundary invariance from invariance of integrated observable results.
- **Section 5.3 — UV running:** use **Fig. 11** as the full 16-model descriptive sensitivity comparison; use **Fig. 12** to show the stage ablation with the low-quality-fit caveat. Introduce the saved UV beta-product breaking diagnostics only if useful to explain the observed UV difference, not as an independent causal proof.
- **Section 5.4 — scalar strength:** use **Fig. 13** and discuss the fitted Yukawa compensation `lambda×Y1×Y2≈constant` and fit limitations. The current fits establish viability trends, **not** a pure high-vs-low scalar perturbation holding Yukawas fixed.
- **Section 5.5 — representations:** use **Fig. 14** to show a controlled **within-case** representation change with shared Yukawas. This is stronger than comparing independent best fits, which are intentionally calibrated to the same targets.
- **Section 5.6 — flavour:** retain previously delivered Figs. 04, 06, 07, 10. The last two show hard/direct matching and the *final* SM+Weinberg neutrino running; physical intermediate PMNS predictions inside the scalar-only EFT are not claimed.

## Minimal remaining work

1. **No full 16-model UV or high-scalar recomputation required.** The records exist.
2. Verify/quote the fit-selection rule and full numerical input table in Chapter 4.
3. If you want one additional stronger UV result: run the **same-seed staged ablation** at a successful `5e-9` optimised benchmark (such as A(-2)), using the existing `RunFittedStagedRGComparison.evaluate_staged` interface and saved artifacts. This is a suggested small validation, **not a prerequisite** for a coherent results chapter.
4. Write the interpretation with stated limitations rather than starting a broader phenomenology study.

## Reproduction (PowerShell, from repository root)

```powershell
Expand-Archive -LiteralPath .\T3_VERIFIED_EXISTING_COMPARISONS_FIGURES.zip -DestinationPath . -Force

# Regenerate four figures and three CSV tables from CURRENT LOCAL saved results:
python -m Numerical.analysis.BuildExistingRGComparisonThesisFigures --source local

# Or reproduce this included verified GitHub snapshot with no other saved data required:
python -m Numerical.analysis.BuildExistingRGComparisonThesisFigures --source snapshot

Start-Process 'Reports\output\thesis_assets\figures\fig13_low_high_scalar_fit_quality.pdf'
```

**No tracked source from the repository is overwritten**, except that extracting into your local repository adds the new `Numerical/analysis/BuildExistingRGComparisonThesisFigures.py` script, one guide under `docs/`, and figures/tables/data under `Reports/output/thesis_assets/`. Nothing changes on GitHub.

# Thesis figure placement — recommended final order

Based on `docs/T3_THESIS_PLAN.md` and the existing saved, optimised **UV-frozen** `lambdaT3 = 5e-9` results.
The contents are **actual saved outputs** except Figure 01 (schematic), Figure 02 (analytical scaling identity), and Figure 10 (vector composition of the already-generated Figures 08 and 09).

**Main narrative:** T3 field choices and BSM couplings determine matching; leading matching has continuous coupling-rescaling degeneracies; different parts of RG evolution can break these degeneracies; neutrino mass splittings can run substantially while PMNS angles barely change, consistent with approximate alignment of the fitted mass-weighted symmetric Yukawa tensors.

## Main thesis: exact placement (one primary location per figure)

| Section in your existing plan | Figure PDF basename | Priority | Place immediately after discussing | Correct interpretation |
|---|---|---|---|---|
| **3.3 Scalar-only EFT and Weinberg generation** | `fig01_t3_effective_theory_flow` | CORE | The three threshold scales and EFT stages | **Schematic:** common fermion decoupling, scalar EFT, Weinberg EFT. Cross-reference in 4.2 rather than repeating. |
| **3.5 Scaling directions and predictions** | `fig02_analytic_coupling_scaling` | OPTIONAL | The transformations `lambda -> gamma lambda`, `Y1,Y2 -> gamma^-1/2 Y1,Y2` | **Analytical identity**, not a fitted or RG-evolved sensitivity plot. An equation alone may suffice. |
| **5.1 Fit viability** | `fig03_optimised_fit_quality_16_models` | OPTIONAL | Defining fit target, `chi2`, and which saved models meet it | Shows fit objective for 16 UV-frozen fitted configurations. Avoid using tiny chi2 differences to rank viable theories. Put the full 16-model table in **Appendix C**. |
| **5.3 Where quantum corrections matter** | `fig07_t3b_c5_threshold_components` | CORE | Explaining `C5 = C5_hard + C5_direct` | Quantifies relative one-loop contributions **at the scalar threshold**, only T3-B(-1). **Not** a comparison of UV on/off. |
| **5.6.1 Neutrino mass-basis alignment** | `fig04_five_class_flavour_alignment` | CORE | Defining the ratio `norm(offdiag(U.T @ C12 @ U))/norm(C12)` | Shows mass-basis near-diagonality for five fitted class representatives. Not an automatic result of Majorana symmetry. |
| **5.6.2 Why symmetrisation matters** | `fig05_t3b_yukawa_product_matrices` | OPTIONAL | Contrasting the raw product `A` with `(A + A.T)/2` | Matrix-by-matrix T3-B(-1) explanation. Keep if you have space; otherwise use it in Appendix D. |
| **5.6.3 Heavy-state cancellations** | `fig06_t3b_heavy_state_cancellations` | CORE | Writing `C12 = sum_r F_r` | Explains small off-diagonals as fitted cancellations across heavy states, not one small Yukawa. |
| **5.6.4 Consequence for observable running** | `fig10_mass_vs_mixing_running_combined` | CORE | Comparing the evolution of mass splittings and PMNS angles | **Combined two-panel plot**; use this **instead of** including figs 08 and 09 separately. Only final SM+Weinberg interval, 1 TeV -> 100 GeV. |
| **5.6.4 alternative** | `fig08_t3b_weinberg_mass_splitting_running` and `fig09_t3b_weinberg_angle_running` | ALTERNATIVE | As above, in two separate floats | Same data as fig10. **Do not insert all three**. |

### Chapter-level reading order

- **Chapter 1 (Introduction):** no figure is essential. State questions and goals; optionally reference, but do not repeat, fig01.
- **Chapter 2 (Background):** no uniquely justified original plot in the saved-data package. A T3 loop topology and representation summary would be useful, but need sourced theoretical diagrams; not fabricated here.
- **Chapter 3 (Analytical):** fig01 and optionally fig02. These show the EFT structure and analytical scaling argument, not numerical validation.
- **Chapter 4 (Numerical methodology):** cross-reference fig01 from §4.2. No need for a second pipeline diagram, unless the chapter requires one illustrating precise numerical solver/fit sequence.
- **Chapter 5 (Results):** fig03 (optional) in §5.1; fig07 in §5.3; figs04, 05 (optional), 06, 10 in §5.6. Figures 08+09 are **alternatives** to 10.
- **Chapter 6 (Conclusions):** no new figure; summarise results and limitations.

## Results sections without sufficient plots **yet**

| Section | Figure to generate once matching evidence is available | Evidence required |
|---|---|---|
| **5.2 Coupling rescaling** | `chi2`, mass splittings and angles **versus** `gamma` for the actual `T1/T2/T3` parameter trajectories | Sensitivity outputs evaluated with/without RG evolution; must distinguish UV-frozen and UV-running scenarios. `fig02` is only an analytical illustration. |
| **5.3 UV evolution effect** | Controlled **UV-on versus UV-frozen** endpoint differences, with identical, matched inputs | Both numerical evaluations. `fig07` does not establish the UV effect. |
| **5.4 Low/high scalar coupling** | `lambdaT3=5e-9` versus `lambdaT3=0.5`: before and after Yukawa refitting | Successful optimised high-scalar output. Plot fixed-Yukawa and refitted comparisons separately. No results are invented here. |
| **5.5 Representation dominance** | Representation-by-representation relative changes with **shared** versus optimised Yukawas | Controlled common-input comparisons, not just optimised chi2 values. |

**Recommended minimum main-thesis set if short on pages:** figs 01, 07, 04, 06, 10 (5 figures). Add fig02 (analytic motivation), fig03 (fit context), fig05 (matrix explanation) if space permits. For a numerical thesis, the future §5.2 and §5.4 figures should become higher priority than optional figs02/03/05.

## Tables: recommended locations

| Table file (LaTeX/CSV) | Intended location | Use |
|---|---|---|
| `table01_uv_frozen_fits_16_models` | **Appendix C**, cited in §5.1 | Sixteen optimized results and target observables. |
| `table02_five_class_flavour_alignment` | **§5.6.1** or Appendix D | Exact off-diagonal fractions for five representative classes; fig04 provides visual summary. |
| `table03_t3b_c5_decomposition` | **§5.3**, alongside fig07 | The three Frobenius norms and hard/direct ratio. |
| `table04_t3b_weinberg_running_endpoints` | **§5.6.4**, alongside fig10 | Actual endpoint values and relative changes of observables. |
| `table05_t3b_heavy_cancellation` | **Appendix D**, cited in §5.6.3 | Signed heavy contribution values underlying fig06. |

## Scientific caveats to retain in the captions/text

- `C5` is **symmetric**, not generally diagonal. The near-diagonality of fitted `C12` is in the eigenbasis of the *full matched* `C5`.
- The tensor `C12` is a mass-weighted, **symmetrised** combination of *two different* Yukawas. The individual Yukawas need not be small, aligned, or diagonal.
- Hard and direct contributions can have different heavy-mass weighting, so exact `C5`-`C12` alignment is not guaranteed.
- The final Weinberg-EFT mass-splitting and angle curves (figs08/09/10) do not show physical PMNS mixing-angle evolution above the scalar threshold.
- Every saved numerical plot is from the **UV-frozen, low-scalar, fitted** sample. Never label it as a verified UV-on/off or `lambdaT3=0.5` comparison.
- Five-class figures use one fitted real-Yukawa benchmark *per class*, not all 16; the other numerical output fig03 covers all 16 fit-quality values.

## Reproducibility and project paths

Extract this ZIP from the project root. **All figures** go in `Reports/output/thesis_assets/figures/`; **tables** in `Reports/output/thesis_assets/tables/`; this guide and the thesis plan in `docs/`. The existing reproducible generator is at `Numerical/analysis/BuildT3ThesisAssets.py`. The added panel-composition script is `Numerical/analysis/MakeT3ThesisMassMixingPanel.py`. Both use saved data; neither modifies production running or fits.

Using PowerShell in the project root:

```powershell
Expand-Archive -LiteralPath .\T3_THESIS_FIGURES_BY_CHAPTER.zip -DestinationPath . -Force
# Figures are ready as PDF+PNG. Optional: reconstruct plots using bundled saved snapshot
python -m Numerical.analysis.BuildT3ThesisAssets --source snapshot
# Optional: compose the side-by-side fig10 from fig08 and fig09 PDF pages
python -m Numerical.analysis.MakeT3ThesisMassMixingPanel
```

For `\includegraphics`, paths below assume the **LaTeX compiler root is the repository root**. If your thesis is a separate LaTeX project, copy the selected figures into its `figures/` folder and adjust those paths. Tables use `\usepackage{booktabs}`, figures use `\usepackage{graphicx}`.

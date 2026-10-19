# Chapter 4 thesis exports

Generated using `Numerical/orchestration/ExportChapter4Figures.py` from saved numerical data only. PNG and PDF versions are exported.

- **4.1:** `figures/fig_4_1_representation_comparison.*`, `tables/table_4_1_baseline_16.csv`. Sixteen common-input predictions, grouped A–E.
- **4.2:** `figures/fig_4_2_rg_corrections.*`, `tables/table_4_2_rg_corrections_16.csv`, `latex/figure_4_2_rg_corrections.tex`. Signed running-versus-frozen effects at factor 1; LaTeX includes the definition equation and figure caption.
- **4.3:** `figures/fig_4_3_fermion_mass_dependence.*`; six-model scan of **both** heavy mass `MF` and matching scale `muF`.
- **4.4:** `figures/fig_4_4_fermion_rg_dependence.*`; fractional effects using the same scan.
- **4.5:** `figures/fig_4_5_coupling_dependence.*`; baseline-normalised **running** delta-m21-squared for Yukawa and lambda scans against k^4 and k^2. This plot deliberately avoids redundant atmospheric panels.
- **4.6:** `figures/fig_4_6_observable_correlations.*`; A-m2, B-m3, D-p0, common-input UV sensitivity directions T1/T2/T3. **Blue regions are approximate Gaussian NuFIT ellipses, not official correlated likelihood contours. Model dash-dot ovals are empirical 2-sigma descriptive envelopes, NOT model confidence regions.** No RGE rerun.
- **4.7:** `figures/fig_4_7_c5_contributions.*`, `tables/table_4_3_c5_contributions_6.csv`. Frobenius norms and matrix-level interference. These are not additive contributions to mass-squared differences.

Labels show $\Delta m_{31}^{2}$ for the stored normal-ordering `delta_m3l_sq_ev2`; original field names are untouched. Class E uses its distinct lambdaT3 normalisation. Figures and tables cover *common-input*, not individually refitted, points.

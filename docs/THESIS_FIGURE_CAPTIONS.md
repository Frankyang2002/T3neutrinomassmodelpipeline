# Suggested thesis figure captions

The paths below are relative to `Reports/output/thesis_assets/figures/`. In LaTeX:

```latex
\usepackage{graphicx}
\usepackage{booktabs}
% In the figure environment:
\centering
\includegraphics[width=0.92\linewidth]{Reports/output/thesis_assets/figures/fig04_five_class_flavour_alignment.pdf}
\caption{...}
\label{fig:t3-flavour-alignment}
```

**Figure 01 — Effective theory flow.** "Hierarchical T3 calculation with simultaneous decoupling of all heavy fermions at $\mu_F$, generation of the mixed Yukawa Wilson coefficient $C_{12}$, direct operator mixing into the Weinberg coefficient in the scalar-only EFT, scalar matching at $\mu_S$, and final SM+Weinberg running. Schematic only."

**Figure 02 — Analytic scaling.** "Leading bilinear matching rescaling $\lambda_{T3}\rightarrow\gamma\lambda_{T3}$ and $Y_{1,2}\rightarrow\gamma^{-1/2}Y_{1,2}$, leaving $\lambda_{T3}Y_1Y_2$ invariant algebraically. Curves are analytic and do not include RGE evolution or constraints on admissible parameters."

**Figure 03 — Saved fits.** "Optimised fits of the 16 T3 variants at the UV-frozen fermion threshold for $\lambda_{T3}=5\times10^{-9}$. The plotted $\chi^2$ is the fitting objective, not a likelihood-based model comparison. The T3-C(-1) and T3-C(-3) points have larger residual objectives than the rest."

**Figure 04 — Flavour alignment.** "Off-diagonal Frobenius norm of the unsymmetrised Yukawa product $A=Y_1^*M_F^{-1}Y_2^\dagger$ and its Majorana-symmetrised part $C_{12}=(A+A^T)/2$, evaluated in the orthogonal eigenbasis of each benchmark's matched $C_5$. Five optimised UV-frozen real-Yukawa benchmarks are shown. Symmetry alone does not enforce diagonality."

**Figure 05 — Yukawa matrices.** "Unsymmetrised Yukawa product versus the symmetrised $C_{12}$ for the UV-frozen T3-B(-1) benchmark. Entries are in units of $10^{-6}\,\mathrm{GeV}^{-1}$, in the same eigenbasis of the complete matched $C_5$. The symmetrisation strongly suppresses off-diagonal entries for this fitted point."

**Figure 06 — Heavy-state cancellation.** "The three heavy-fermion contributions to the off-diagonal $C_{12}$ entries in the eigenbasis of the matched T3-B(-1) $C_5$. Significant cancellations in their sum explain why the symmetrised mass-weighted Yukawa product is nearly diagonal in this basis."

**Figure 07 — Scalar-threshold components.** "Frobenius norms of the hard-matching and direct $C_{12}\rightarrow C_5$ contributions to the matched Weinberg coefficient for the UV-frozen T3-B(-1) benchmark. The norm of the sum differs from the sum of norms because of interference. These are one-loop fixed-order contributions, not a coupled Wilson-evolution result."

**Figure 08 — Mass-squared running.** "Relative changes in neutrino mass-squared splittings over the final SM+Weinberg EFT interval, $\mu=1\,\mathrm{TeV}$ to $100\,\mathrm{GeV}$, for the optimised UV-frozen T3-B(-1) benchmark. Both splittings decrease by approximately $8.08\%$ in the saved numerical trajectory."

**Figure 09 — Angle running.** "Changes in mixing angles over the same SM+Weinberg EFT interval and benchmark as Figure 08. The largest absolute drift is below $10^{-4}$ degrees. These are Weinberg-stage physical predictions; the scalar-only interval is not represented here by a physical PMNS trajectory."

## Model-comparison caveats

Avoid saying: "We proved the running is flavour-invariant", "Majorana symmetrisation forces a diagonal matrix", or "strong scalar coupling always makes $C_{12}$ running larger." The plots show an approximately aligned **fitted** flavour tensor for selected benchmarks and extremely small PMNS drift in the final Weinberg stage. High-coupling and UV-running comparisons need separate numerical evidence.

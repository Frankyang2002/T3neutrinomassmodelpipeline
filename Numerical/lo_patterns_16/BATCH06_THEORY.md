# Chat A Batch 6 — heavy-resolved CG × scalar-loop factorization

## Objective and source scope

Test whether the all-16 conditional one-loop amplitude is proportional, at the **coefficient-tensor** level, rather than merely after selecting special fitted Yukawa textures. This uses the existing candidate left-handed Weyl vertex convention from `scripts/BuildT3DirectOneLoopInputs.py` and the exact restricted 2×2 scalar blocks from Batch 3. It does **not** establish a physically normalized Matchete/RGBeta amplitude, modify production code, or perform EFT.

All results are conditional on the input map, and the 0/16 physical-normalization certificate is unchanged. The complete 16-model source hashes and Batch-1/3/4/5 lineage are checked before writing output.

## 1. Remove Yukawas before comparing representations

The stored source Weyl vertices are of the form

\[
Z_{iIs}=\sum_{r=1}^3Y_{1,ir}^*C^{(1)}_{rIs},\qquad
Y_{jIt}=\sum_{r=1}^3Y_{2,jr}^*C^{(2)}_{rIt}.
\]

The component tensors \(C^{(1,2)}\) are extracted from each frozen model's existing *full-rank* fitted Yukawas by linear solving, exactly as in Batch 5. Their reconstruction is checked against original vertices; no group factors are invented. \(I,J\) are candidate Weyl fermion components, \(s,t\) original scalar components, \(r\) a heavy-generation index, and \(i\) an external lepton-flavour index.

For a heavy singular-mass eigenspace \(g\), let \(P_g\) be the invariant spectral projector of \(MM^\dagger\). The independently reconstructed heavy chirality-flip numerator is \(M^\dagger P_g\). For each equal-charge mixing pair \(q=(s_q,t_q)\), write

\[
\mathcal{M}^2_{S,q}=\begin{pmatrix}a_q&\delta_q^*\\\delta_q&b_q\end{pmatrix},\quad
\lambda_{q\pm}=\frac{a_q+b_q}{2}\pm\sqrt{(a_q-b_q)^2/4+|\delta_q|^2}.
\]

The exact loop matrix-function cross entry is

\[
K_q(M_g^2)=\delta_q\,\frac{B_0(0;M_g^2,\lambda_{q+})-B_0(0;M_g^2,\lambda_{q-})}{\lambda_{q+}-\lambda_{q-}}.
\]

Then define the **Yukawa-independent**, potentially complex \(3\times3\) source-generation coefficient

\[
B^{(R)}_{g,rs}=\sum_q K_q(M_g^2)
  C^{(1)}_{r,I,s_q}(M^\dagger P_g)_{IJ}C^{(2)}_{s,J,t_q}.
\]

The complete candidate mass matrix is *exactly* (for the restricted source assumptions)

\[
(M_\nu^{(R)})_{ij}=\frac{1}{16\pi^2}\sum_g\bigg[
  \bigl(Y_1^* B_g^{(R)}Y_2^{*T}\bigr)_{ij}+
  \bigl(Y_1^* B_g^{(R)}Y_2^{*T}\bigr)_{ji}\bigg].
\]

The code reconstructs all 16 full mass matrices and **each** of their three heavy-eigenspace contributions from these coefficients. This is a distinct contraction pathway from Batch 3's direct evaluation with fixed lepton vertices.

## 2. What constitutes a factorization theorem at fixed spectra?

Choose reference model \(R_0\), here `A_-2`, and evaluate the complex proportionality coefficient across **all 27 complex entries** (3 heavy groups × 3 × 3 source-generation indices):

\[
\kappa_R^{\text{best}}=\frac{\sum_{g,r,s} B_{g,rs}^{(0)*}B_{g,rs}^{(R)}}{\sum_{g,r,s}|B_{g,rs}^{(0)}|^2},\quad
\varepsilon_R=\frac{\|\{B_g^{(R)}-\kappa_R B_g^{(0)}\}_g\|_F}{\|\{B_g^{(R)}\}_g\|_F}.
\]

Additionally fit independent \(\kappa_{R,g}\) in each group. Equality of group scales matters. In the limit of vanishing \(\varepsilon_R\) and group-scale spread, the relation

\[
B_g^{(R)}=\kappa_R B_g^{(0)}\quad\forall g
\]

**algebraically implies** proportional mass matrices \(M_\nu^{(R)}=\kappa_R M_\nu^{(0)}\) **for any pair of Yukawa matrices** substituted into the same candidate vertex convention, at that exact set of scalar and heavy masses. It does *not* prove the relation for arbitrary masses/potentials, much less prove a physical amplitude equality across different gauge representations.

By construction these source-generation coefficients can distinguish two mass matrices that happen to share their singular values or that happen to be proportional only at one fitted Yukawa texture.

## 3. Why factorization is plausible, and when it can fail

For small mixing,

\[
K_q(M_g^2)=\delta_q\,
\frac{B_0(0;M_g^2,a_q)-B_0(0;M_g^2,b_q)}{a_q-b_q}
+O(\delta_q|\delta_q|^2).
\]

If \(a_q=a\) and \(b_q=b\) for all contributing charges, the divided difference is a common function of \(M_g\). The candidate CG/yukawa tensor contractions may therefore be absorbed into representation coefficients. If two charge sectors have identical \(|\delta_q|\), the *exact* scalar functions also coincide; two aligned terms can add as in the E frozen benchmarks.

Failure modes include **non-identical scalar mass splittings between charge components**, inequivalent scalar mixing magnitudes so that the nonlinear loop function does not factor, and heavy-state coefficients that depend differently on \(M_g\) in Dirac versus Majorana constructions. The numerical report evaluates rather than presupposes these effects.

## 4. Five source-conditional mass scenarios

| Scenario | In-memory change | Interpretation |
|---|---|---|
| `baseline` | Common magnitude \(|\lambda_{T3}|=0.10\); E mixing ×2 | Compare existing representation coefficients at controlled mixing; does not change saved frozen E point |
| `scalar_common_shift` | Multiply all S1 diagonal squared masses by 1.28 and S2 by 0.78 | Universal mass-parameter change preserving multiplet degeneracy |
| `heavy_generation_tilt` | Multiply sorted distinct Takagi masses by 0.71, 1.19, 1.64 | Vary heavy mass **eigenspaces**, preserving Weyl mass symmetry and their projectors |
| `large_scalar_mixing` | Multiply S1–S2 mixing by 20 | Stress the \(O(\delta|\delta|^2)\) corrections, provided no tachyons occur |
| `charge_split_diagnostic` | Split scalar component diagonal masses by index, holding original vertices and mixing | Deliberately outside the restricted SU(2)-multiplet-degenerate mass ansatz. A charge-conserving mass-matrix algebra exercise, **not** a demonstrated gauge-invariant UV deformation |

The exact scalar kernel and invariant heavy projectors are recalculated independently for every scenario. All original JSONs are left intact. The report checks the exact tensor formula against Batch 3's independently constructed full result at **each of 80 model/scenario points**.

## 5. Mirror relations and model-by-model conclusions

Alongside comparison of all 16 models to `A_-2`, the report compares seven pairs from Batch 5: the A/D exchanges, B and C reflection pairs, and both E reflected pairs. These are now measured at the heavy-coefficient level instead of only final neutrino-matrix level. If a pair remains proportional after nondegenerate mass changes, the restricted source evidence for a representation-level relation is stronger. If it fails for `charge_split_diagnostic`, that identifies precisely why the multiplet-degenerate approximation matters.

The report contains **actual** heavy-group complex multipliers, one coefficient matrix per heavy group (all numerical entries in `batch06_factorization.json`), scalar loop kernels by charge, linear-vs-exact errors and specific E-channel norms. Conclude only what those computed numbers support. Cross-package B1/B2 phase and normalization certification is separate and remains unresolved in Chat A.

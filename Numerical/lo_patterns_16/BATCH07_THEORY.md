# Chat A — Batch 7: factorization theorem, CG weights and corrections

## Scope

This batch consolidates Batches 1–6. It does **not** introduce EFT matching, change the UV input builder or certify the Matchete/RGBeta-to-Weyl dictionary. All physics statements are **conditional on the candidate Weyl vertices**, restricted scalar quadratic, fixed renormalization convention and selected benchmark/deformation scenarios. The physical-normalization certificate remains **0/16**.

The read-only runner uses **the saved Batch-6 coefficient tensors**, the saved Batch-1 frozen source inputs and the source-audited CG table. It does not refit Yukawas and does not call the original direct UV evaluator. It does recompute closed-form scalar kernels and their derivatives independently in memory to test the approximation.

## 1. The tensor behind the flavour pattern

For heavy Weyl mass matrix \(M_F=M_F^T\), let \(P_g\) denote a projector onto the complete eigenspace of the same eigenvalue of \(M_FM_F^\dagger\). This projector is unambiguous under rotations of degenerate states. Define the component contraction

\[
T_{qg,rs}=C^{(1)}_{r,I,s_q}(M_F^\dagger P_g)_{IJ}C^{(2)}_{s,J,t_q},
\qquad
B^{(R)}_{g,rs}=\sum_q K_q(M_g^2)T_{qg,rs}.
\]

For each isolated equal-charge scalar mixing block,

\[
\mathcal M_{S,q}^2=\begin{pmatrix}a_q&\delta_q^*\\\delta_q&b_q\end{pmatrix},
\quad \lambda_{q\pm}={a_q+b_q\over2}\pm\sqrt{(a_q-b_q)^2/4+|\delta_q|^2},
\]

\[
K_q(M_g^2)=\delta_q\frac{B_0(0;M_g^2,\lambda_{q+})-B_0(0;M_g^2,\lambda_{q-})}{\lambda_{q+}-\lambda_{q-}}.
\]

The **candidate** source-generation tensor contributes to the full neutrino matrix through

\[
M_\nu^{(R)}=\frac{1}{16\pi^2}\sum_g \Big[Y_1^* B_g^{(R)}Y_2^{*T}+(Y_1^* B_g^{(R)}Y_2^{*T})^T\Big].
\]

The complex bilinear expression and **the absence of any separate physical four-spinor normalization** are intentional; do not infer the latter from this formula.

## 2. Fixed-spectrum factorization theorem

**Statement.** For two models at a *matched mass-spectrum experiment*, suppose their extracted Yukawa-independent coefficient tensors obey

\[
B_g^{(R)}=\kappa_R B_g^{(0)}\qquad\forall\;g
\]

with **one** complex constant \(\kappa_R\) shared by all distinct heavy-mass eigenspaces. Then for any common Yukawas \(Y_1,Y_2\),

\[
M_\nu^{(R)}(Y_1,Y_2)=\kappa_RM_\nu^{(0)}(Y_1,Y_2),
\qquad \frac{m_i^{(R)}}{m_j^{(R)}}=\frac{m_i^{(0)}}{m_j^{(0)}}
\]

when the denominator is nonzero. This follows by linearity of the symmetrized bilinear contraction. **This is an algebraic theorem about the candidate coefficient tensors**, not an assertion that different gauge representations have identical physical amplitudes. The converse is not needed, and matching singular-value ratios alone does not prove this hypothesis.

For an approximate equality, put \(E_g=B_g^{(R)}-\kappa_RB_g^{(0)}\). The deterministic Frobenius bound is

\[
\|\Delta M_\nu\|_F\le\frac{2}{16\pi^2}\|Y_1\|_2\|Y_2\|_2\sum_g\|E_g\|_F.
\]

Batch 6 measured the concatenated \(B_g\) coefficient residual across **all 3 heavy eigenspaces × 3 × 3 entries**, not merely a final flavour-matrix match.

## 3. All 16 source-audited group weights

These expressions are from `output/general_observable/lo_direct_one_loop/source_audits/charge_resolved_loop_weights.json`, whose convention is `Y1 * Y2 * conjugate(T3MixCG)` on the corresponding source weight-basis indices. Their **signs are not physically certified four-spinor amplitudes**.

| Model | Active electric charges | Per-channel group weight | Total source group weight |
|---|---|---|---|
| A_-4 | −2 | −2√3/3 | −2√3/3 |
| A_-2 | −1 | −2√3/3 | −2√3/3 |
| A_+0 | 0 | −2√3/3 | −2√3/3 |
| B_-3 | −1 | −2√3/3 | −2√3/3 |
| B_-1 | 0 | −2√3/3 | −2√3/3 |
| B_+1 | +1 | −2√3/3 | −2√3/3 |
| C_-3 | −1 | +2/3 | +2/3 |
| C_-1 | 0 | +2/3 | +2/3 |
| C_+1 | +1 | +2/3 | +2/3 |
| D_-2 | 0 | −2√3/3 | −2√3/3 |
| D_+0 | +1 | −2√3/3 | −2√3/3 |
| D_+2 | +2 | −2√3/3 | −2√3/3 |
| E_-4 | −2, −1 | +√6/3 each | +2√6/3 |
| E_-2 | −1, 0 | +√6/3 each | +2√6/3 |
| E_+0 | 0, +1 | +√6/3 each | +2√6/3 |
| E_+2 | +1, +2 | +√6/3 each | +2√6/3 |

Consequently, **magnitudes** relative to A_-2 are (A/B/D):1, C:\(1/\sqrt3\), E:\(\sqrt2\). The measured Batch-6 coefficient scales are close to these ratios when the scalar mixing is controlled at \(|\lambda_{T3}|=0.10\); the small deviations are numerical/source-conditional, not new group-theory constants. The audit phase and the candidate Weyl fitted complex scale are distinct objects, pending B1/B2.

## 4. Why multiplet degeneracy supports factorization

If the unmixed scalar squared masses satisfy \(a_q=a\), \(b_q=b\) for all active channels, the first divided-difference loop function is common across charges:

\[
K_{qg}=\delta_qD_g(a,b)+O(\delta_q|\delta_q|^2),\qquad
D_g(a,b)=\frac{f_g(a)-f_g(b)}{a-b},\quad f_g(x)=B_0(0;M_g^2,x).
\]

Then the **charge sum** acts on the CG/fermion tensor and may reduce to one representation coefficient. Crucially, scalar degeneracy by itself does not mathematically enforce \(B_g^{(R)}\propto B_g^{(0)}\): the fermion contractions, Weyl orientation and mixing coefficients must also have the required structure. This is precisely what Batch 6 measured rather than assuming.

The baseline coefficient residuals across 16 models reported by the local Batch-6 run were \(0\)–\(4.53\times 10^{-10}\); this is numerical evidence, not an exact symbolic identity.

## 5. Derivation of the first corrections

Set \(a_q=a+\epsilon_{a,q}\) and \(b_q=b+\epsilon_{b,q}\). For nondegenerate \(a\ne b\), differentiation of the scalar divided difference gives

\[
\partial_aD_g=\frac{f'_g(a)-D_g}{a-b},\qquad
\partial_bD_g=\frac{D_g-f'_g(b)}{a-b}.
\]

The mixing expansion of the exact two-scalar matrix function is odd in \(\delta\):

\[
K_q=\delta_q\left[D_g+|\delta_q|^2E_g+O(|\delta_q|^4)\right],\qquad
E_g=\frac{f'_g(a)+f'_g(b)-2D_g}{(a-b)^2}=f_g[a,a,b,b].
\]

Combining the two perturbations,

\[
\boxed{\begin{aligned}
B_{g,rs}&=\sum_q\delta_qD_g T_{qg,rs} \\
&\quad+\sum_q\delta_q(\epsilon_{a,q}\partial_aD_g+\epsilon_{b,q}\partial_bD_g)T_{qg,rs}\\
&\quad+\sum_q\delta_q|\delta_q|^2E_gT_{qg,rs}+\cdots.
\end{aligned}}
\]

The omitted terms include \(O(\delta\epsilon^2)\), \(O(\delta|\delta|^2\epsilon)\) and \(O(\delta|\delta|^4)\). Nondegenerate scalar masses are assumed for the displayed derivative expressions; the exact loop function has a smooth degenerate limit but is not used in this derivative truncation.

A charge-dependent splitting makes \(\partial_aD_g\) and \(\partial_bD_g\) reweight different charge channels. Their \(M_g\) dependence generally prevents absorbing the correction into one common \(\kappa_R\). Nonlinear mixing similarly reweights \(|\delta_q|^3\) rather than \(\delta_q\).

## 6. Test design and interpretation

The runner compares directly to saved Batch-6 **exact** coefficient tensors at both `charge_split_diagnostic` and `large_scalar_mixing`, with all three invariant heavy groups retained. It tests the first-order scalar shift around the exact baseline and the leading-versus-cubic mixing expansion.

The five Batch-6 stress maxima provided locally are:

| Scenario | Largest coefficient residual |
|---|---:|
| baseline | 4.52958914e−10 |
| scalar_common_shift | 4.35904865e−10 |
| heavy_generation_tilt | 2.36521908e−9 |
| large_scalar_mixing | 1.81783236e−7 |
| charge_split_diagnostic | 2.23487435e−4 |

The charge-splitting scenario is an **artificial algebraic scalar mass deformation**, not a demonstrated SU(2)-invariant source coupling. Its role is to find the limit of the restricted-mass factorization statement. The reported relative approximation errors are not physical uncertainties.

## 7. Reproducibility and open physics

Input lineage is checked by SHA-256 against saved Batch-1/5/6 reports and the original frozen candidate Weyl inputs. A changed audit convention, missing model, source-drifted JSON, false certification flag, or populated output directory is rejected. Nothing is written until all calculations complete successfully.

Physical source B1/B2 certification, UV coupling normalization, allowed general scalar potential, and eventual all-fermion-then-all-scalar EFT matching are **outside this batch**.

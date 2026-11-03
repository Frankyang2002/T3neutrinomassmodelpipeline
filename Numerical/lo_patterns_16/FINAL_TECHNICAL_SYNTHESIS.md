# Chat A — All-16 T3 radiative neutrino masses: final analytical synthesis

**Scope and status.** This document consolidates the eight completed Chat A batches for the **restricted direct-UV one-loop** T3 setup at the frozen reference scale \(\mu_0=10^7\,\mathrm{GeV}\). All numeric results below were **reported by the user's Windows runs**. This packaged narrative has not independently rerun the frozen 16-model calculation; the companion local verification program reads the actual saved outputs and generates a source-checked version with a machine-readable manifest. **Physical cross-package B1/B2 normalisation certification: 0/16.** Neither the absolute mass predictions nor inter-representation amplitude phases are claimed as source-certified physical predictions. No EFT or matching calculation is included.

## Executive scientific result

In the candidate Weyl convention and restricted scalar sector, the sixteen one-loop neutrino mass matrices admit an accurately factorized description at matched Yukawas and the controlled multiplet-degenerate spectrum:

\[
\boxed{M_\nu^{(R)}(Y_1,Y_2)\simeq \kappa_R\,\mathcal F(Y_1,Y_2;\{M_g\},a,b)}.
\]

Here \(\kappa_R\) is a model-dependent coefficient (whose **magnitude ratios** follow the source-audited CG contractions), while \(\mathcal F\) contains the dependence on the Yukawa flavour matrices and common heavy-state loop functions. The algebraic implication is exact **if** the underlying heavy-eigenspace coefficient tensors are proportional with a *single* complex multiplier for every heavy eigenspace. The condition is measured, not postulated, for the restricted benchmarks. Thus the frozen models' different normalised spectra primarily reflect their **different fitted Yukawas**. Deliberate charge-dependent scalar mass splitting breaks proportionality at first order; nonlinear scalar mixing gives higher-order corrections.

The key measured numbers from the previous local runs are:

| Test | User-reported result | Interpretation |
|---|---:|---|
| Provisional UV evaluations | 16/16 | Candidate Weyl one-loop amplitudes available |
| Exact scalar-kernel matrix reconstruction | 16/16; maximum relative error \(2.10\times10^{-11}\) | Confirms finite two-scalar kernel implementation against saved UV matrices |
| Batch-5 matched-texture full-matrix pairs | 28/28; largest residual \(1.48\times10^{-14}\) | Matrix proportionality for seven pairs and four texture choices |
| Batch-6 all-model heavy-coefficient fit | Maximum baseline residual \(4.529589\times10^{-10}\) | Yukawa-independent factorization at matched controlled spectrum |
| Batch-7 CG magnitude vs coefficient scale | Largest absolute discrepancy \(3.094433\times10^{-8}\) | Candidate group-weight magnitude pattern agrees numerically |
| Batch-8 scalar-splitting remainder order | 2.0047–2.0149 | Correct quadratic Taylor remainder along the tested path |
| Batch-8 transverse onset order | 1.0016–1.0084 (B–E) | Linear onset of nonproportionality |
| Physically source-certified amplitudes | **0/16** | B1/B2 cross-package dictionary not settled in this workstream |

All error measures here are **algebraic/numerical**, not theory uncertainties. The scalar-splitting path is a controlled, charge-conserving **diagnostic deformation**, not an established gauge-invariant extension of the Lagrangian.

## 1. Fields, representations, and restricted interaction

The five representation families and sixteen allowed \(\alpha\) values are:

| Family | \((d_{S_1},d_{S_2},d_F)\) | \(\alpha\) values | Heavy-fermion assignment |
|---|---|---|---|
| A | (1,3,2) | −4, −2, 0 | Dirac |
| B | (2,2,1) | −3, −1, +1 | Majorana only at −1; otherwise Dirac |
| C | (2,2,3) | −3, −1, +1 | Majorana only at −1; otherwise Dirac |
| D | (3,1,2) | −2, 0, +2 | Dirac |
| E | (3,3,2) | −4, −2, 0, +2 | Dirac |

The original interaction topology is represented schematically by \(y_1\bar L F^c S_1+y_2\bar L F S_2^\dagger+\lambda_{T3}HH S_1S_2^\dagger+\mathrm{h.c.}\). Index ordering, complex conjugation and signs used in this study follow the **candidate** all-left-handed Weyl builder. The restricted benchmark sets additional, unsupplied symmetry-allowed scalar couplings to zero. Its scalar bilinear \(B\) is zero, and the unmixed squared masses are constant within each scalar multiplet. The frozen family E uses \(|\lambda_{T3}|=0.05\), while A–D use \(0.10\); *controlled comparisons* rescale E mixing to 0.10. Frozen Yukawa matrices differ between benchmarks. Consequently, raw frozen masses must not be interpreted as pure representation effects.

## 2. The exact restricted one-loop scalar kernel

For each allowed equal-electric-charge scalar pair, the candidate scalar squared-mass block is

\[
\mathsf A_q=\begin{pmatrix}a_q&\delta_q^*\\\delta_q&b_q\end{pmatrix},\qquad
\lambda_{q\pm}=\frac{a_q+b_q}{2}\pm \frac12\sqrt{(a_q-b_q)^2+4|\delta_q|^2}.
\]

Here \(\delta_q\propto\lambda_{T3}v^2\), with the proportionality determined by the source CG tensor and the chosen Weyl convention. Define \(f_g(x)=B_0(0;M_g^2,x)_{\overline{\mathrm{MS}},\mathrm{finite}}\). The exact cross-propagator entry is

\[
\boxed{K_q(M_g^2)=\delta_q\frac{f_g(\lambda_{q+})-f_g(\lambda_{q-})}{\lambda_{q+}-\lambda_{q-}}.}
\]

It is smooth in the equal-mass limit. In this convention the finite zero-momentum function is

\[
f_g(x)=1-\ln\frac{M_g^2}{\mu^2}-\frac{x\ln(x/M_g^2)}{x-M_g^2},
\]

with its coincident-mass value obtained continuously. The common additive reference-scale dependence cancels in the divided difference. For \(|\delta_q|\) small compared with the squared-mass separation,

\[
D_g(a,b)=\frac{f_g(a)-f_g(b)}{a-b},\quad
K_q=\delta_qD_g(a,b)+\delta_q|\delta_q|^2E_g(a,b)+O(\delta_q|\delta_q|^4),
\]

\[
E_g(a,b)=\frac{f'_g(a)+f'_g(b)-2D_g(a,b)}{(a-b)^2}.
\]

These are **mass-squared** divided differences, so \(D_g\) has units \(\mathrm{GeV}^{-2}\). The cubic term improves the finite-mixing stress test relative error from roughly \(5.66\!\times\!10^{-5}\)–\(7.55\!\times\!10^{-5}\) to \(4.78\!\times\!10^{-8}\)–\(8.50\!\times\!10^{-8}\) (the previous run's large-mixing deformation); this does not assign a physical error bar.

## 3. Source-audited charge and CG table

The audit convention is the source contraction `Y1 * Y2 * conjugate(T3MixCG)` on the matched weight indices. Its signs **are not yet certified four-spinor amplitude signs**. It is valid here as a source-convention group-theory table.

| Model | Active internal electric charges \(Q\) | Candidate weight per channel | Candidate total |
|---|---|---|---|
| A_-4 | −2 | \(-2\sqrt3/3\) | \(-2\sqrt3/3\) |
| A_-2 | −1 | \(-2\sqrt3/3\) | \(-2\sqrt3/3\) |
| A_+0 | 0 | \(-2\sqrt3/3\) | \(-2\sqrt3/3\) |
| B_-3 | −1 | \(-2\sqrt3/3\) | \(-2\sqrt3/3\) |
| B_-1 | 0 | \(-2\sqrt3/3\) | \(-2\sqrt3/3\) |
| B_+1 | +1 | \(-2\sqrt3/3\) | \(-2\sqrt3/3\) |
| C_-3 | −1 | \(+2/3\) | \(+2/3\) |
| C_-1 | 0 | \(+2/3\) | \(+2/3\) |
| C_+1 | +1 | \(+2/3\) | \(+2/3\) |
| D_-2 | 0 | \(-2\sqrt3/3\) | \(-2\sqrt3/3\) |
| D_+0 | +1 | \(-2\sqrt3/3\) | \(-2\sqrt3/3\) |
| D_+2 | +2 | \(-2\sqrt3/3\) | \(-2\sqrt3/3\) |
| E_-4 | −2, −1 | \(+\sqrt6/3\) each | \(+2\sqrt6/3\) |
| E_-2 | −1, 0 | \(+\sqrt6/3\) each | \(+2\sqrt6/3\) |
| E_+0 | 0, +1 | \(+\sqrt6/3\) each | \(+2\sqrt6/3\) |
| E_+2 | +1, +2 | \(+\sqrt6/3\) each | \(+2\sqrt6/3\) |

Relative to A_-2, the **magnitudes** are

\[
\boxed{|\kappa_A|:|\kappa_B|:|\kappa_C|:|\kappa_D|:|\kappa_E|
=1:1:1/\sqrt3:1:\sqrt2.}
\]

Only E has two active topology charge channels. Batch 3 found that under its restricted multiplet-degenerate inputs, the two E channel **mass matrices** are aligned with Frobenius cosine 1 and each contributes approximately half the total. This is constructive *charge-channel* interference; it does not preclude destructive *heavy-eigenspace* interference.

## 4. Heavy-eigenspace coefficient theorem

Let \(M=M^T\) be the candidate Weyl heavy mass matrix, \(P_g\) the spectral projector onto a complete eigenspace of \(MM^\dagger\), and \(C^{(1,2)}\) the component tensors extracted from the source-derived Weyl vertices (not independently physically certified). Let \(s_q,t_q\) index a mixing pair in the original scalar basis. Define

\[
T^{(R)}_{qg,rs}=C^{(1)}_{rI s_q}(M^\dagger P_g)_{IJ} C^{(2)}_{sJ t_q},
\quad B^{(R)}_{g,rs}=\sum_qK_q(M_g^2)T^{(R)}_{qg,rs}.
\]

In the restricted input convention,

\[
\boxed{M_\nu^{(R)}=\frac1{16\pi^2}\sum_g\left[Y_1^*B_g^{(R)}Y_2^{*T}
+\left(Y_1^*B_g^{(R)}Y_2^{*T}\right)^T\right].}
\]

**Conditional factorization theorem:** if, for models \(R,R_0\) at matched spectra,

\[
B_g^{(R)}=\kappa_RB_g^{(R_0)}\quad\text{for every }g,
\]

with one complex \(\kappa_R\) independent of \(g\), then for **every common Yukawa pair**

\[
M_\nu^{(R)}=\kappa_R M_\nu^{(R_0)},\quad
m_i^{(R)}=|\kappa_R|m_i^{(R_0)},\quad
\frac{m_i^{(R)}}{m_j^{(R)}}=\frac{m_i^{(R_0)}}{m_j^{(R_0)}}.
\]

The proof is simply linearity of the preceding symmetrized bilinear map. For deviations \(E_g=B_g^{(R)}-\kappa_RB_g^{(R_0)}\), the rigorous Frobenius bound is

\[
\|\Delta M_\nu\|_F\le\frac{2}{16\pi^2}\|Y_1\|_2\|Y_2\|_2
\sum_g\|E_g\|_F.
\]

Batch 6 fits \(\kappa_R\) to the **27 complex entries** of the three 3×3 \(B_g\) tensors using \(\langle B_0,B_R\rangle/\langle B_0,B_0\rangle\). Maximum residual \(4.53\times10^{-10}\), and maximum spread of independently fitted heavy-group scales \(5.95\times10^{-10}\), at controlled baseline. This is an *empirical approximate satisfaction* of the theorem's hypothesis, **not a symbolic proof** that the CG identity holds at all masses or in all general potential terms. Changing an overall coefficient cannot change mass ratios; near-equal spectral ratios alone would have been insufficient evidence, which is why Batch 6 reconstructs \(B_g\) directly.

## 5. Frozen neutrino spectra and flavour patterns

The following eV masses are user-reported sorted Takagi singular values for the **sixteen different fitted frozen points** (not an ordering assignment and not physical-source-certified):

| Model | \(m_1\) [eV] | \(m_2\) [eV] | \(m_3\) [eV] | \(m_1/m_3\) | Heavy-state \(C_1\) |
|---|---:|---:|---:|---:|---:|
| A_-4 | 0.000474483 | 0.00758273 | 0.0436006 | 0.010882 | 0.057116 |
| A_-2 | 0.00721475 | 0.0110085 | 0.0485602 | 0.14857 | 0.480645 |
| A_+0 | 0.00721427 | 0.0110000 | 0.0484990 | 0.14875 | 0.480074 |
| B_-3 | 0.00721823 | 0.0109814 | 0.0483324 | 0.14935 | 0.480121 |
| B_-1 | 0.00452249 | 0.00980338 | 0.0504413 | 0.089659 | 0.208537 |
| B_+1 | 0.00721583 | 0.0109461 | 0.0480743 | 0.15010 | 0.505219 |
| C_-3 | 0.0269469 | 0.0281923 | 0.0550407 | 0.48958 | 1.000000 |
| C_-1 | 0.0374973 | 0.0385889 | 0.0622850 | 0.60203 | 0.285778 |
| C_+1 | 0.000365488 | 0.00826871 | 0.0477861 | 0.0076484 | 0.0131221 |
| D_-2 | 0.00721417 | 0.0109978 | 0.0484824 | 0.14880 | 0.478505 |
| D_+0 | 0.00721467 | 0.0110070 | 0.0485499 | 0.14860 | 0.480577 |
| D_+2 | 0.000474733 | 0.00758174 | 0.0435948 | 0.010890 | 0.0571179 |
| E_-4 | 0.00331940 | 0.00741694 | 0.0384885 | 0.086244 | 0.2139999 |
| E_-2 | 0.00471419 | 0.00867883 | 0.0424695 | 0.11100 | 1.000000 |
| E_+0 | 0.00464358 | 0.00864859 | 0.0425154 | 0.10922 | 1.000000 |
| E_+2 | 0.00332597 | 0.00744248 | 0.0386336 | 0.086090 | 0.213222 |

The frozen \(m_3\) span is approximately 0.03849–0.06229 eV. This range is **benchmark-dependent**, not a universal representation ranking. The close A/D reflected pairs arise with very similar fitted Yukawa matrices; their full-matrix frozen proportionality residuals are about \(2.81\times10^{-5}\)–\(5.40\times10^{-5}\), and fall to \(\sim10^{-16}\)–\(10^{-15}\) at matched textures. By contrast the frozen C_-3/C_+1 pair residual is \(0.98844\), dropping to \(\lesssim6\times10^{-16}\) under matched textures. This differentiates a **fitted-texture effect** from an overall source-conditional representation coefficient.

## 6. Why C_+1 is highly hierarchical: invariant heavy-mode projections

Write \(M_\nu=\sum_gH_g\) with \(g\) a **complete distinct heavy-mass eigenspace**, never an arbitrarily chosen eigenvector inside a degenerate space. For a singular pair \(u_i,v_i\), define

\[
a_{gi}=u_i^\dagger H_gv_i,\quad
G=\sum_g\|H_g\|_F,\quad
A_i=\frac{\sum_g|a_{gi}|}{G},\quad
C_i=\frac{|\sum_g a_{gi}|}{\sum_g|a_{gi}|}.
\]

Then

\[
\boxed{m_i=G A_i C_i.}
\]

\(A_i\) quantifies the projection of contributions onto the chosen singular mode; \(C_i\) quantifies their complex-phase coherence. Neither is a new observable: they are *basis-consistent diagnostics of the candidate decomposition*, and individual projections can be sensitive to nearly degenerate light singular values. The source-conditional C-family data were:

| C model | Fermion | \(A_1\) | \(C_1\) | \(m_1/m_3\) | Main diagnostic |
|---|---|---:|---:|---:|---|
| C_-3 | Dirac | 0.26878 | 1.00000 | 0.48958 | Heavy terms aligned; no cancellation of lightest-mode projections |
| C_-1 | Majorana | 0.56910 | 0.285778 | 0.60203 | Partial cancellation despite larger overall amplitude |
| C_+1 | Dirac | 0.31920 | 0.0131221 | 0.0076484 | Strong destructive heavy-eigenspace interference |

Thus \(C_+1\) retains only about 1.312% of its sum of projected amplitude magnitudes. The observed lightest mass ratio is

\[
\frac{m_1(C_{-1})}{m_1(C_{+1})}\approx102.6
\approx\underbrace{2.64}_{G\text{ ratio}}\times
\underbrace{1.78}_{A_1\text{ ratio}}\times
\underbrace{21.78}_{C_1\text{ ratio}}.
\]

The large **coherence** factor supplies the strongest numerical distinction. However C_-1 is Majorana while C_+1 is Dirac **and** the frozen Yukawas differ. This decomposition does not by itself attribute the hierarchy to one physical cause. Batch 5's *counterfactual* matched-texture swaps found \(m_1/m_3=0.0076484\) under C_+1 Yukawas and \(0.60203\) under C_-1 Yukawas in **all sixteen** controlled models. Therefore the chosen texture is sufficient to determine the hierarchy in this restricted calculation even when the underlying fermion construction is held fixed within each recipient model. These swaps are diagnostic calculations, not new fitted or physical models.

## 7. How factorization fails: first scalar corrections

Let \(a_q=a+t\epsilon_{a,q}\), \(b_q=b+t\epsilon_{b,q}\), keeping \(\delta_q\) fixed. Expanding the **leading-in-mixing** kernel gives

\[
\begin{split}
B_g^{(R)}(t)={}&\sum_q\delta_qD_gT^{(R)}_{qg}
+t\sum_q\delta_q(\epsilon_{a,q}\partial_aD_g+\epsilon_{b,q}\partial_bD_g)T^{(R)}_{qg}\\
&+\sum_q\delta_q|\delta_q|^2E_gT^{(R)}_{qg}+\cdots,\\
\partial_aD_g={}&\frac{f'_g(a)-D_g}{a-b},\qquad
\partial_bD_g=\frac{D_g-f'_g(b)}{a-b}.
\end{split}
\]

The remainder includes \(O(t^2\delta)\), \(O(t\delta|\delta|^2)\), and \(O(\delta|\delta|^4)\). The derivatives depend on heavy mass \(M_g\) differently from the common \(D_g\), so charge-dependent weights generally cease to be representable by one \(\kappa_R\).

**Important accuracy distinction.** The Batch-7 first-split approximation error (0.000309–0.000953 **relative to the entire exact coefficient tensor**) is not the fractional error in the much smaller factorization-breaking component (order 0.00002–0.0002). Therefore Batch 7 alone does **not** validate the split correction's transverse accuracy.

Batch 8 takes the derivative of the **exact finite-mixing kernel** rather than its leading term. Define \(R(t)=\{B_g^{(A_{-2})}(t)\}_g\), \(B(t)=\{B_g^{(R)}(t)\}_g\),

\[
\kappa(t)=\frac{\langle R(t),B(t)\rangle}{\langle R(t),R(t)\rangle},\qquad
D_R(t)=B(t)-\kappa(t)R(t).
\]

Accounting for the derivative of \(\kappa(t)\) is essential; a fixed multiplier would confuse overall rescaling with genuine factorization breaking. The exact derivatives obey

\[
B(t)=B(0)+t B'(0)+O(t^2),\qquad
D_R(t)=D_R(0)+t D_R'(0)+O(t^2).
\]

User-reported stress-test results (maximum residual in each scenario):

| Scenario | \(\max_R\epsilon_R\) | Interpretation |
|---|---:|---|
| Controlled baseline | \(4.529589\times10^{-10}\) | Very strong coefficient factorization |
| Common scalar shift | \(4.359049\times10^{-10}\) | Shared loop function changes; relation preserved |
| Nonuniform heavy masses | \(2.365219\times10^{-9}\) | Relation approximately preserved |
| Large scalar mixing (×20) | \(1.817832\times10^{-7}\) | Odd nonlinear mixing corrections resolvable |
| Charge-dependent scalar split | \(2.234874\times10^{-4}\) | Dominant factorization breaking on this diagnostic path |

The measured full-tensor Taylor remainder order at \(t=1/8\to1/4\) was 2.0047–2.0149. Where the transverse signal exceeds the numerical floor, its growth order was 1.0016–1.0084. A models remained at the nonproportional noise floor for this **particular** splitting path; this is not a universal symmetry protection.

| Family | Measured full remainder order | Measured transverse order | Endpoint nonproportionality |
|---|---:|---:|---:|
| A | 2.00592 | Noise floor | ≤\(7.0\times10^{-16}\) |
| B | 2.01225 | 1.00193 | \(1.084\times10^{-4}\) |
| C | 2.01225 | 1.00193 | \(1.084\times10^{-4}\) |
| D | 2.01488 | 1.00841 | \(2.235\times10^{-4}\) |
| E | 2.00466 | 1.00165 | \(2.241\times10^{-5}\) |

The independent exact-derivative-versus-central-difference test had maximum relative discrepancy \(1.5868\times10^{-5}\), below the test threshold. These measured powers establish local convergence along the specified diagnostic path; they do **not** make the splitting deformation gauge-invariant or certify the UV normalization.

## 8. Precise scope, open questions and handoff

**Established algebraically:** the exact restricted two-scalar matrix-function divided difference; the implication \(B_g^{(R)}=\kappa_RB_g^{(0)}\;\forall g\implies M_\nu^{(R)}=\kappa_RM_\nu^{(0)}\) at matched Yukawas; the Frobenius error bound; the lightest-mode amplitude identity; Taylor expansion of the exact finite-mixing kernel.

**Numerically verified in user-reported frozen runs:** 16/16 provisional source-Weyl amplitudes, exact kernel reconstructions, E channel alignment, matched-texture matrix proportionality, candidate CG magnitude pattern, heavy-eigenspace coefficient proportionality at the restricted baseline, and the selected scalar-splitting convergence orders.

**Not established:** an all-parameter SU(2)-invariant factorization theorem, gauge-invariance of the ad hoc charge splitting, a physically normalized Matchete/RGBeta-to-Weyl amplitude, certified source phases or possible different convention-dependent symmetry factors, or phenomenological viability of the *fitted* spectra. The 0/16 physical certification status is not changed by these mathematical/numerical checks. This report is **not** an EFT calculation.

**Provenance / reproducibility.** Local source files are `output/general_observable/lo_pattern_study16/provisional_01` and the seven `analysis_batch02_01` through `analysis_batch08_01` directories. Candidate group weights originate in `output/general_observable/lo_direct_one_loop/source_audits/charge_resolved_loop_weights.json`. Numerical audit outputs are not replicated in this report; the companion `BuildLOFinalSynthesis16.py` checks all eight results, source audit/CG rows and SHA links, then creates a fresh generated evidence table, source manifest and handoff. Individual full-precision complex mass matrices and source hashes remain in the original evidence directories.

**Future use:** preserve this as a provisional interpretation. Only after independent B1/B2 source-convention certification and an authorised regeneration should absolute normalisations and physically meaningful relative phases be updated. Do not conflate this source-conditional notebook theorem with full UV or SMEFT matching.

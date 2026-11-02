# Restricted T3: 16 conditional UV patterns and analytic checks

**Scope:** fixed `mu0 = 10^7 GeV`, restricted A–E/16 frozen benchmarks only. No EFT, running or thresholds. This is a reproducibility and conditional numerical study, **not** a physical Matchete/RGBeta/Weyl normalisation certificate. GitHub `main` inspected at `0b210f18f6a2cfb10fff03059e2ec4a73d2195a4`; the installed runner records the *actual local* Git SHA and source hashes.

## Source-derived numerical conventions

`Numerical/direct_one_loop/core.py` accepts the complete candidate left-Weyl action

\[
\mathcal L_Y=-\nu_a\psi_i(Y_{ais}\phi_s+Z_{ais}\phi_s^*)+\mathrm{h.c.},\qquad
\mathcal L_m=-\tfrac12\psi_iM_{ij}\psi_j+\mathrm{h.c.},
\]
\[
V_2=\phi^\dagger A\phi+\tfrac12(\phi^TB\phi+\mathrm{h.c.}),
\quad A=A^\dagger,\quad B=B^T.
\]

With \(\psi_i=U_{ik}\chi_k\), \(U^TMU=\operatorname{diag}(M_k)\), and
\(\phi_s=C_{sr}R_r\), \(R_r\) canonically real, the actual code evaluates

\[
g_{akr}=\sum_{is}U_{ik}(Y_{ais}C_{sr}+Z_{ais}C^*_{sr}),
\qquad
(m_\nu)_{ab}=\frac1{16\pi^2}\sum_{kr}M_k\,g_{akr}g_{bkr}
 B_0^{\overline{\rm MS}}(0;M_k^2,m_r^2).
\]

The last sum is over all declared vertices; `core.evaluate` makes unordered **vertex** pairs and explicitly adds the transposed mixed-vertex contribution, not an unexplained extra factor two. \(B_0=1-\log(M_k^2/\mu^2)-x\log x/(x-1)\), \(x=m_r^2/M_k^2\). The UV pole is proportional to \(\sum_{kr}M_k g_{akr}g_{bkr}\) and must cancel to produce a finite radiative mass. `native_uv_loop_reference.independent_reference` independently reconstructs the real Hessian and applies \(M^\dagger f(MM^\dagger)\) to eliminate reliance on a second Takagi implementation. A synthetic Majorana/complex-scalar split reference checks the \(1/(16\pi^2)\) loop convention and the **additional 1/2 from a real scalar component decomposition**, without claiming that the Matchete-to-Weyl dictionary has been verified.

`BuildT3DirectOneLoopInputs.py` provisionally maps `y1_HC` into **`Z`** (conjugate scalar) and `y2_HC` into **`Y`** (non-conjugate scalar), and adds scalar mixing as

\[
A_{S_2(t),S_1(s)}=-\frac{\lambda_{T3}v^2}{2}
\,[T3MixCG]_{H^0H^0,s,t},\qquad A_{S_1(s),S_2(t)}=A_{S_2(t),S_1(s)}^*.
\]

This is the **candidate input-builder orientation**, not independently verified physical Feynman-rule normalisation. Its potential-to-Hessian sign and the extra `Matchete`/`RGBeta` phase/normalisation dictionary remain pending. The original neutral Higgs component in the builder is indexed `[1,1]` with zero-based indexing. The restricted frozen models have `B=0`; extra symmetry-allowed operators absent from the restricted benchmark are set to zero.

## First-order scalar-mixing formula

Write \(A=A_0+\Delta\) with block-diagonal \(A_0=\operatorname{diag}(A_{S_1},A_{S_2})\), \(\Delta_{S_2,S_1}\propto\lambda_{T3}\), and \(B=0\). In a basis where the *unmixed* scalar blocks are diagonal, define \(\widehat Z_{aks}=\sum_i Z_{ais}U_{ik}\), \(\widehat Y_{akt}=\sum_i Y_{ait}U_{ik}\). At first order,

\[
(m_\nu)_{ab}^{(1)}=\frac1{16\pi^2}
\sum_{k,\,s\in S_1,\,t\in S_2}M_k
(\widehat Z_{aks}\widehat Y_{bkt}+\widehat Z_{bks}\widehat Y_{akt})
\,\Delta_{t s}\,\mathcal D_k(m_s^2,m_t^2),
\]
\[
\mathcal D_k(a,b)=\frac{B_0(0;M_k^2,a)-B_0(0;M_k^2,b)}{a-b},
\quad\mathcal D_k(a,a)=\partial_a B_0(0;M_k^2,a).
\]

This expression has the correct **input-basis** scalar orientation \(Z(S_1^*)\,\Delta_{S_2,S_1}\,Y(S_2)\). If a nontrivial scalar diagonalisation is needed within the unmixed blocks, transform the source-component vertices and \(\Delta\) into that basis before applying the formula. Nonzero common charge, nonvanishing Y1/Y2 CG and a nonvanishing mixing entry are all necessary for a component contribution. The full expression includes flavour-dependent heavy eigenmasses, finite split-scalar masses, and coherent interference.

For *this* restricted topology (with no same-vertex loop terms), scaling just one Yukawa is exactly linear, \(m_\nu(ty_1,y_2)=t m_\nu(y_1,y_2)\), likewise for \(y_2\), including complex phases and a negative \(t\). The scalar mixing has an odd-power expansion \(m_\nu(\lambda)=c_1\lambda+c_3\lambda^3+\cdots\) for this `B=0` block structure; in particular `t=-1` reverses the full complex mass while leaving Takagi masses unchanged. Linear \(\lambda\) scaling is only approximate at finite mixing. Unlike the Yukawas, changing the heavy mass moves both the explicit mass insertion and the loop-kernel arguments and has no universal exact power law.

## Native candidate CG contraction and all 16 charge channels

The existing `charge_resolved_loop_weights.json` audits **`T3Y1CG * T3Y2CG * conjugate(T3MixCG)`** in the same source-weight indices. These are the **audit's component tensor ordering and conjugation**, not a source-certified physical four-spinor loop. Its `indices_weight_basis_one_based` are preserved in `channel_weights.json`; neither the loop prefactor nor a cross-package field rephasing has been folded into these weights. For each row below, all other potential source charges have **zero topology weight** under the restricted channel selection, even when an individual Yukawa vertex alone supports them. The output includes the equal-charge S1/S2 component pairs with zeros, and the exact original-component charge partitions from `channels.partition_exact`.

| Model | Fermion | Allowed topology scalar charges \(Q\) | Native candidate weight per channel | Group-only total | Conditional interpretation |
|---|---|---|---|---|---|
| A_-4 | Dirac doublet | -2 | \(-2/\sqrt3\) | \(-2/\sqrt3\) | One doubly charged scalar-mixing channel; other source-component charges cannot complete topology |
| A_-2 | Dirac doublet | -1 | \(-2/\sqrt3\) | \(-2/\sqrt3\) | One singly charged channel; same CG magnitude as A_-4 |
| A_+0 | Dirac doublet | 0 | \(-2/\sqrt3\) | \(-2/\sqrt3\) | One neutral scalar channel; same candidate group factor |
| B_-3 | Dirac singlet | -1 | \(-2/\sqrt3\) | \(-2/\sqrt3\) | A single channel with B's doublet scalar pair |
| B_-1 | **Majorana singlet** | 0 | \(-2/\sqrt3\) | \(-2/\sqrt3\) | Neutral Majorana internal fermion; source spinor normalisation still pending |
| B_+1 | Dirac singlet | +1 | \(-2/\sqrt3\) | \(-2/\sqrt3\) | Positive charged channel, not a new group multiplicity |
| C_-3 | Dirac triplet | -1 | \(+2/3\) | \(+2/3\) | Singly charged scalar channel with triplet fermion representation metric |
| C_-1 | **Majorana triplet** | 0 | \(+2/3\) | \(+2/3\) | Neutral Majorana triplet; triplet metric and relative phase remain important |
| C_+1 | Dirac triplet | +1 | \(+2/3\) | \(+2/3\) | Single positive charged channel; candidate smaller CG magnitude than B |
| D_-2 | Dirac doublet | 0 | \(-2/\sqrt3\) | \(-2/\sqrt3\) | Triplet S1 and singlet S2 mix only in a neutral channel |
| D_+0 | Dirac doublet | +1 | \(-2/\sqrt3\) | \(-2/\sqrt3\) | One positively charged component channel |
| D_+2 | Dirac doublet | +2 | \(-2/\sqrt3\) | \(-2/\sqrt3\) | One doubly positive charged component channel |
| E_-4 | Dirac doublet | -2, -1 | \(+\sqrt6/3\) each | \(+2\sqrt6/3\) | Two allowed negative charge channels; coherent relative phase and loop kernels matter |
| E_-2 | Dirac doublet | -1, 0 | \(+\sqrt6/3\) each | \(+2\sqrt6/3\) | One charged and one neutral channel; weights alone cannot predict interference |
| E_+0 | Dirac doublet | 0, +1 | \(+\sqrt6/3\) each | \(+2\sqrt6/3\) | One neutral and one charged channel |
| E_+2 | Dirac doublet | +1, +2 | \(+\sqrt6/3\) each | \(+2\sqrt6/3\) | Two positive charge channels; topology multiplicity 2 in **CG basis only** |

All entries are **native source-audit candidate weights**, not physically verified CG loop coefficients. In particular, `E`'s two equal weights cannot be replaced by an unconditional factor of two in its physical neutrino mass. Model-dependent masses, flavour textures, accidental cancellations and scalar mixing can reverse any hierarchy suggested by group-only totals.

The distinct model family representations are

- A: \((d_{S1},d_{S2},d_F)=(1,3,2)\), \(\alpha=-4,-2,0\).
- B: \((2,2,1)\), \(\alpha=-3,-1,+1\).
- C: \((2,2,3)\), \(\alpha=-3,-1,+1\).
- D: \((3,1,2)\), \(\alpha=-2,0,+2\).
- E: \((3,3,2)\), \(\alpha=-4,-2,0,+2\).

With the paper's integer hypercharge convention as in the source catalog, the component electric charge for weight index \(w=0,\ldots,d-1\) is \(Q=(d-1)/2-w+Y\), where `S1` has \(Y=\alpha/2\), `S2` has \(Y=(\alpha+2)/2\), and `F` has \(Y=(\alpha+1)/2\). This is not the usual half-normalized SU2 expression with hypercharge added as \(Y/2\) unless the symbol `Y` is translated accordingly. All matching is performed before heavy Takagi rotations, since doubled-Dirac Takagi eigenstates need not carry a unique original-component charge.

## Two comparisons, deliberately separate

1. **Frozen benchmark view** (`all16_summary.csv`, `models/*.json`): actual fixed-scale 3x3 complex matrices and sorted Takagi eigenmasses from each supplied frozen benchmark, a coupled numerical observation *conditional* on its unverified Weyl mapping. Mismatch between two reference algorithms is reported independently. A model is blocked rather than fabricated when its input is absent or inconsistent.
2. **Formal common-kernel group view** (`common_kernel_group_only.csv`): set all loop functions, flavour Yukawa factors and off-diagonal mixing factors equal *by hypothesis* and compare only \(\sum_q W_q\). Then \(|W|/|W_A|=1\) for A/B/D, \(1/\sqrt3\) for C and \(\sqrt2\) for E. **This is not a physically legal common normalisation comparison until B1/B2 independently map the representations and their couplings.** Use the controlled within-model scans to test the kernel dependence instead of treating these formal ratios as neutrino-mass predictions.

## Reproducibility and negative controls

The runner reads 16 frozen benchmark JSONs and the corresponding 16 `provisional_inputs/{model}.json`. It validates the per-input recorded benchmark SHA256, representation, scale, source executable guards, source CG audit model coverage and component charge conservation. The per-model report records matrix in GeV/eV, masses, scalar/fermion eigenmasses, both pole diagnostics, reference/production differences, cancellation conditioning, and the exact original-charge partition. `individual_terms` from production use mass eigenstate indices; their sum is checked but individual eigenstates are not incorrectly assigned unique electric charges. Source-frozen charge partitions are the correct basis for channel interpretation.

`scans.csv` contains `t=-1,0,0.5,1,2` for each of `y1`, `y2`, `lambdaT3`, and `t=0.5,1,2` for full heavy Weyl-mass scaling, with invalid/tachyonic points blocked. `matrix_scaling_residual=||m(t)-t m(1)||/||m(1)||` distinguishes exact Yukawa scaling from only leading-order scalar-mixing scaling. Changing the sign of a complete complex matrix leaves Takagi masses invariant, so the check uses complex matrices rather than eigenmass ratios alone. The reference's pole subtraction and the production's absolute term sum expose numerically ill-conditioned cancellations.

**Physical evidence required before final numbers:** independent B1 witness of Matchete four-spinor to normalized left-Weyl Feynman rules (including possible Majorana symmetry factors and Dirac SU2 epsilon phases), B2 witness of SU2 component and RGBeta/Matchete coupling/dagger/portal mapping, and a source-and-benchmark hash-anchored coordinator certificate. This study refuses to promote a Boolean in an arbitrary JSON file into physical certification. The CLI `--physical-dictionary` intentionally fails closed until the coordinator certificate's verifiable schema exists.

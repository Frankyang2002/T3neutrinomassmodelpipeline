# Chat A Batch 7 — Conditional factorization theorem and first corrections

**Status:** derived and numerically checked within the frozen candidate Weyl dictionary. Physical Matchete/RGBeta normalization and phases are **not certified (0/16)**. No EFT matching or production-code changes.

## 1. Source-convention input and exact matrix identity

The candidate source vertices are \(Z=Y_1^*C_1\) and \(Y=Y_2^*C_2\), and the restricted scalar quadratic has \(B_{SS}=0\) with isolated 2×2 equal-charge blocks. For a full heavy-degenerate eigenspace projector \(P_g\),

\[\boxed{B^{(R)}_{g,rs}=\sum_q K_q(M_g^2)\,C^{(1)}_{r,I,s_q}(M_F^\dagger P_g)_{IJ}C^{(2)}_{s,J,t_q}}\]

\[K_q=\delta_q\,{f_g(\lambda_{q+})-f_g(\lambda_{q-})\over\lambda_{q+}-\lambda_{q-}},\quad f_g(x)=B_0(0;M_g^2,x).\]

\[\boxed{M_\nu^{(R)}={1\over16\pi^2}\sum_g\operatorname{Sym}\left(Y_1^*B_g^{(R)}Y_2^{*T}\right)},\quad\operatorname{Sym}(X)=X+X^T.\]

The coefficients B_g are 3×3 **source-generation** matrices, not neutrino-mass matrices. Their independent heavy-group projections avoid arbitrary phase choices inside degenerate fermion eigenspaces.

## 2. Fixed-spectrum factorization theorem

**Theorem.** If for all heavy groups \(B_g^{(R)}=\kappa_R B_g^{(0)}\) with one common complex \(\kappa_R\), then for *any* common 3×3 Yukawa pair

\[M_\nu^{(R)}(Y_1,Y_2)=\kappa_R M_\nu^{(0)}(Y_1,Y_2),\]

as follows immediately from linearity of the symmetrized bilinear contraction. The statement is exact for the specified coefficient tensors; a small measured residual only establishes approximate factorization at the checked spectra.

For errors \(E_g=B_g^{(R)}-\kappa_RB_g^{(0)}\), the Frobenius bound is

\[\|\Delta M_\nu\|_F\le {2\over16\pi^2}\|Y_1\|_2\|Y_2\|_2\sum_g\|E_g\|_F.\]

This is a deterministic bound for arbitrary common Yukawa matrices, not a fitted-texture observation. The theorem does not assert a physical equality between differently normalized theories.

## 3. Source CG-weight table (16 models)

The listed weights use the **source audit** contraction `Y1 * Y2 * conjugate(T3MixCG)`, not a certified four-spinor vertex rule. The source-audit signed phase must **not** automatically be identified with the fitted Weyl coefficient phase κ.

| Model | Q of active S1/S2 pairs | Source weight(s) | Total weight | |weight|/|A_-2| | measured |κ| | baseline tensor ε |
|---|---|---|---|---:|---:|---:|
| A_-4 | -2 | -2*sqrt(3)/3 | -2*sqrt(3)/3 | 1.000000000 | 1.000000000 | 6.453e-16 |
| A_-2 | -1 | -2*sqrt(3)/3 | -2*sqrt(3)/3 | 1.000000000 | 1.000000000 | 0.000e+00 |
| A_+0 | 0 | -2*sqrt(3)/3 | -2*sqrt(3)/3 | 1.000000000 | 1.000000000 | 6.358e-16 |
| B_-3 | -1 | -2*sqrt(3)/3 | -2*sqrt(3)/3 | 1.000000000 | 1.000000025 | 4.530e-10 |
| B_-1 | 0 | -2*sqrt(3)/3 | -2*sqrt(3)/3 | 1.000000000 | 1.000000025 | 4.530e-10 |
| B_+1 | 1 | -2*sqrt(3)/3 | -2*sqrt(3)/3 | 1.000000000 | 1.000000025 | 4.530e-10 |
| C_-3 | -1 | 2/3 | 2/3 | 0.577350269 | 0.577350284 | 4.530e-10 |
| C_-1 | 0 | 2/3 | 2/3 | 0.577350269 | 0.577350284 | 4.530e-10 |
| C_+1 | 1 | 2/3 | 2/3 | 0.577350269 | 0.577350284 | 4.530e-10 |
| D_-2 | 0 | -2*sqrt(3)/3 | -2*sqrt(3)/3 | 1.000000000 | 1.000000000 | 5.688e-16 |
| D_+0 | 1 | -2*sqrt(3)/3 | -2*sqrt(3)/3 | 1.000000000 | 1.000000000 | 5.130e-16 |
| D_+2 | 2 | -2*sqrt(3)/3 | -2*sqrt(3)/3 | 1.000000000 | 1.000000000 | 1.492e-15 |
| E_-4 | -1, -2 | sqrt(6)/3, sqrt(6)/3 | 2*sqrt(6)/3 | 1.414213562 | 1.414213531 | 3.928e-10 |
| E_-2 | 0, -1 | sqrt(6)/3, sqrt(6)/3 | 2*sqrt(6)/3 | 1.414213562 | 1.414213531 | 3.928e-10 |
| E_+0 | 1, 0 | sqrt(6)/3, sqrt(6)/3 | 2*sqrt(6)/3 | 1.414213562 | 1.414213531 | 3.928e-10 |
| E_+2 | 2, 1 | sqrt(6)/3, sqrt(6)/3 | 2*sqrt(6)/3 | 1.414213562 | 1.414213531 | 3.928e-10 |

**CG magnitude versus measured tensor scale:** largest absolute ratio mismatch 3.094e-08. The finite mismatch reflects conditional scalar-kernel effects/numerical precision at this benchmark; neither source audit nor fitted κ certifies a physical phase.

## 4. First factorization-breaking terms

Write \(a_q=a+\epsilon_{a,q}\), \(b_q=b+\epsilon_{b,q}\) and \(\delta_q\) for the cross-scalar mixing. For nondegenerate unmixed scalar squared masses define

\[D_g(a,b)={f_g(a)-f_g(b)\over a-b},\quad D_{g,a}={f_g\!\prime(a)-D_g\over a-b},\quad D_{g,b}={D_g-f_g\!\prime(b)\over a-b}.\]

The exact off-diagonal two-scalar loop kernel has the expansion

\[\boxed{K_{qg}=\delta_q\{D_g+\epsilon_{a,q}D_{g,a}+\epsilon_{b,q}D_{g,b}+|\delta_q|^2E_g+\cdots\}},\]

\[E_g={f_g\!\prime(a)+f_g\!\prime(b)-2D_g\over(a-b)^2}=f_g[a,a,b,b].\]

The terms retained in this displayed expression are first order in the charge-dependent diagonal shifts and cubic in mixing; mixed terms such as \(\delta|\delta|^2\epsilon\) and second-order splittings are omitted.

For each heavy group and charge, let \(T_{qg,rs}=C^{(1)}_{r,I,s_q}(M_F^\dagger P_g)_{IJ}C^{(2)}_{s,J,t_q}\). Then the leading splitting correction is

\[\boxed{\Delta B_{g,rs}^{\rm split} =\sum_q\delta_q(\epsilon_{a,q}D_{g,a}+\epsilon_{b,q}D_{g,b})T_{qg,rs}},\]

and the first nonlinear-mixing term is \(\Delta B_g^{\rm cubic}=\sum_q\delta_q|\delta_q|^2E_gT_{qg}\). These corrections need not be proportional to the reference tensor if they weight different charge sectors or heavy groups differently.

## 5. Numerical evidence and perturbative checks

| Scenario | Largest B_g factorization residual (all 16) |
|---|---:|
| baseline (worst B_-1) | 4.529589e-10 |
| scalar_common_shift (worst B_-1) | 4.359049e-10 |
| heavy_generation_tilt (worst C_-3) | 2.365219e-09 |
| large_scalar_mixing (worst C_-3) | 1.817832e-07 |
| charge_split_diagnostic (worst D_+2) | 2.234874e-04 |

The following compares the analytic **truncated** first-correction formulae with exact saved Batch-6 heavy-group coefficient tensors. All errors are relative Frobenius errors; a perturbative approximation is not assumed exact under the deliberate stress tests.

| Model | Charge-split first-order error / exact | |ΔB_split| / |B_baseline| | Large-mix leading error | Large-mix cubic error |
|---|---:|---:|---:|---:|
| A_-4 | 3.250e-04 | 2.505e-03 | 6.541e-05 | 6.370e-08 |
| A_-2 | 3.250e-04 | 2.505e-03 | 6.541e-05 | 6.370e-08 |
| A_+0 | 3.250e-04 | 2.505e-03 | 6.541e-05 | 6.370e-08 |
| B_-3 | 4.132e-04 | 8.579e-03 | 7.554e-05 | 8.495e-08 |
| B_-1 | 4.132e-04 | 8.579e-03 | 7.554e-05 | 8.495e-08 |
| B_+1 | 4.132e-04 | 8.579e-03 | 7.554e-05 | 8.495e-08 |
| C_-3 | 4.132e-04 | 8.579e-03 | 7.554e-05 | 8.495e-08 |
| C_-1 | 4.132e-04 | 8.579e-03 | 7.554e-05 | 8.495e-08 |
| C_+1 | 4.132e-04 | 8.579e-03 | 7.554e-05 | 8.495e-08 |
| D_-2 | 9.533e-04 | 1.465e-02 | 6.541e-05 | 6.370e-08 |
| D_+0 | 9.533e-04 | 1.465e-02 | 6.541e-05 | 6.370e-08 |
| D_+2 | 9.533e-04 | 1.465e-02 | 6.541e-05 | 6.370e-08 |
| E_-4 | 3.091e-04 | 1.253e-03 | 5.664e-05 | 4.777e-08 |
| E_-2 | 3.091e-04 | 1.253e-03 | 5.664e-05 | 4.777e-08 |
| E_+0 | 3.091e-04 | 1.253e-03 | 5.664e-05 | 4.777e-08 |
| E_+2 | 3.091e-04 | 1.253e-03 | 5.664e-05 | 4.777e-08 |

## 6. Domain of validity and limitations

- The 16 input states use a restricted scalar quadratic: no scalar-B terms, and the same unmixed S1 and S2 masses for components of each multiplet.
- With identical charge-sector kernels and compatible heavy-generation contractions, representation dependence can collapse to an overall candidate CG coefficient. Scalar degeneracy alone does **not** establish the necessary tensor identity; Batch 6 checked it numerically.
- The heavy-generation tilt changes three eigenmasses while keeping the spectral projectors; it is a mass-matrix experiment, not a UV renormalization.
- Large scalar mixing resums higher odd powers in δ; the cubic expansion is conditional on small enough |δ| and nondegenerate unmixed masses.
- Charge splitting deliberately violates the multiplet-degenerate restricted mass assumption; it is **not** asserted to arise from an SU(2)-invariant T3 potential.
- Frozen Yukawa differences explain the normalized spectra, but the CG coefficients in this report are not standalone physical amplitudes.
- The source audit uses a separate CG contraction convention; its signs cannot be equated to the candidate Weyl coefficient signs without B1/B2.
- No source physical normalization certification is claimed: **0/16**.

## 7. Scientific synthesis

At matched Yukawas and within the original restricted scalar-sector assumptions, the saved heavy-generation coefficient tensors of all 16 models are nearly proportional, with baseline errors at most 4.530e-10. The candidate CG magnitude pattern is (A,B,D):1, C:1/√3, E:√2 relative to A_-2. Consequently the *normalized* neutrino mass spectrum is controlled by the common Yukawa flavour texture at the checked spectrum. The leading analytically identified violations are charge-dependent scalar mass shifts and nonlinear scalar mixing. This is a **source-conditional numerical-plus-algebraic statement**, not a universal UV theorem or EFT result.

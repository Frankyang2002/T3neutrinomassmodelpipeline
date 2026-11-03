# Track B A/D/E — source-derived Dirac broken-phase UV calculation

## 1. Exact scope and evidence boundary

Models: `A_-4,A_-2,D_-2,D_+0,D_+2,E_-4,E_-2,E_+0,E_+2`. `A_+0` is **excluded**. This work neither changes the repository's production EFT or numerical code nor calls matching, RGE or integration-out procedures. GitHub `main` observed at `d2163e6bb95109c0fac419f7e41f335f0d347b62`; the installed Windows HEAD must be recorded independently. A native CG export made previously with Matchete is a measured tensor, but **not** a fresh native quadratic-action or spinor-amplitude certificate.

The source inspected was `Lagrangian/model/{T3ModelCatalog,T3Fields,SU2Invariants,LagrangianBuilder}.wl`, `Lagrangian/interactions/{T3Topology,ScalarPotential}.wl`, `output/general_observable/lo_direct_one_loop/source_{cg,portal_cg}_all16/`, and the numerical fitted-reference benchmark JSONs. The archived native full-index identity span audit is `output/general_observable/lo_direct_one_loop/source_audits/native_portal_identity_decomposition.json`. Neither provisional input tensors nor `Numerical/direct_one_loop/core.py` enter the independent calculation.

## 2. SU(2) basis and all charged channels

Source hypercharges use `Q=T3+Y`, with `YS1=alpha/2`, `YS2=(alpha+2)/2`, and `YF=(alpha+1)/2`. The scalar triplet is originally in Matchete's Cartesian basis. The saved native generator has

\[
T_3^{(3)}=\begin{pmatrix}0&0&-i\\0&0&0\\i&0&0\end{pmatrix},\qquad
U_3=\begin{pmatrix}-i/\sqrt2&0&i/\sqrt2\\0&1&0\\1/\sqrt2&0&1/\sqrt2\end{pmatrix},\quad
U_3^\dagger T_3^{(3)}U_3=\mathrm{diag}(1,0,-1).
\]

The doublet uses `T3=diag(1/2,-1/2)` and an identity basis transform. The family-dependent internal fermion charges are

\[
Q_{F_1}=\frac{\alpha+2}{2},\qquad Q_{F_2}=\frac\alpha2.
\]

`A` has `S1` singlet and contributes only the second fermion charge; `D` has `S2` singlet and contributes only the first; `E` contributes **both**, including charged fermions and charge-matched scalar components. No neutral-only shortcut or frozen multiplicity is used.

From the **native** Cartesian tensors extracted with `GetCGTensor` (fresh) or the saved actual Matchete `uv_cg_registry` (conditional), for a charge-matched scalar weight pair `(i,j)` define

\[
 c_{1f}=\sum_a (\mathrm{T3Y1CG})_{L_1,f,a}(U_1)_{ai},\qquad
 c_{2f}=\sum_b (\mathrm{T3Y2CG})_{L_1,f,b}(U_2)^*_{bj},
\]
\[
 c_{mf}=\sum_{a,b}(\mathrm{T3MixCG})_{H_2,H_2,a,b}(U_1)_{ai}(U_2)^*_{bj}.
\]

The conjugation here follows the **field objects** in the source operators: `y1 Bar[L] CConj[F] S1`, `y2 Bar[L] F Bar[S2]`, and `lambdaT3 H H S1 Bar[S2]`; it is not fitted from provisional mass matrices. Each H.c. is included separately. The source Lagrangian has original Yukawa coefficient `+1` and Dirac mass `-MF Bar[F] F` (no `1/2`). The general `T3Fields.wl` indexed-mass option `T3_EXPERIMENTAL_FLAVOR_MASSES=1` supports `MF[r]` with a shared flavour index but cannot by itself certify three numerical values. Here the benchmark explicitly supplies diagonal masses `100,110,120 TeV`.

## 3. Scalar action and entire portal-invariant span

The bare scalar quadratic action is `L_free=-MS1^2 S1dag S1-MS2^2 S2dag S2`. The complete native Higgs portal family is

\[
{\cal L}_{\rm portal} = \frac12\operatorname{PlusHc}\left[\sum_{k,n}\lambda_{Hk,\mathrm{Inv}n}
  H^\dagger H S_k^\dagger S_k\, C^{(k,n)}\right]
\]

so at `H=(0,v/sqrt2)` the **mass** diagonal block is

\[
A_k=M_{Sk}^2\mathbf1-\frac{v^2}{2}\sum_n \lambda_{Hk,\mathrm{Inv}n}
 (C^{(k,n)})_{H_2,H_2}.
\]

At the **zero-portal source point** all such portal coefficients are zero, preserving the fitted-reference `MS1Sq, MS2Sq` as the bare squared masses, with the source `lambdaT3` kept nonzero. At the separate **conditional benchmark-identity point**, solve the real-linear **full** tensor equation

\[
\sum_n a_n C^{(k,n)}_{ijab}=\delta_{ij}\delta_{ab},
\quad \lambda_{Hk,\mathrm{Inv}n}=a_n\lambda_{Hk}^{\rm benchmark}.
\]

The recorded native archived tensors give singlet `Inv1` coefficient `1`, triplet `Inv1` coefficient `0` and `Inv2` coefficient `1`. **These assignments are an observed complete-tensor decomposition, not an Inv1/Inv2 naming assumption.** The source benchmark's physical correspondence to a particular RGBeta coupling basis remains conditional and must not be certified from a geometric equality.

The source mixing interaction has coefficient `+1`, **not** `+1/2`:

\[
{\cal L}_{\rm mix}=\lambda_{T3}HH S_1\bar S_2\,\mathrm{T3MixCG}+\mathrm{h.c.},
\qquad A_{S_2,S_1}^{(Q)}=-\frac{\lambda_{T3}v^2}{2}c_{mQ}.
\]

Each relevant charge block is Hermitian,

\[
A_Q=\begin{pmatrix}a_Q&A_{S_2,S_1}^{(Q)*}\\A_{S_2,S_1}^{(Q)}&b_Q\end{pmatrix},\quad
s_\pm=\frac{a_Q+b_Q}{2}\pm\sqrt{\frac{(a_Q-b_Q)^2}{4}+|A_{S_2,S_1}^{(Q)}|^2}.
\]

The program also constructs the **full** `dS1+dS2` Hessian, partitions by electric charge, checks that all off-charge mixing vanishes, and rejects tachyonic scalar sectors.

## 4. Independent one-loop self-energy and spinor qualification

In the explicitly declared (not native-phase-certified) two-component convention `F=(chi,bar(xi))`, `C=+i gamma^2 gamma^0`, `L_Y=-nu(chi*g1*S1bar+xi*g2*S2)+h.c.`,

\[
g_{1,ar}^{(Q)}=-[y_{1,ar}c_{1Q}]^*,\qquad
 g_{2,ar}^{(Q)}=-[y_{2,ar}c_{2Q}]^*,
\]
\[
(m_\nu)_{ab}=\sum_{Q,r}\frac{M_r}{16\pi^2}
 (A_Q)_{S_2S_1}\mathcal D(M_r^2;s_{Q-},s_{Q+})
 [g_{1,ar}^{(Q)}g_{2,br}^{(Q)}+g_{2,ar}^{(Q)}g_{1,br}^{(Q)}].
\]

The **two different** numerical methods are

\[
\mathcal D=\frac{B_0^{\rm fin}(0;M^2,s_+)-B_0^{\rm fin}(0;M^2,s_-)}{s_+-s_-},\quad
 B_0^{\rm fin}=1+\frac{s\log(s/M^2)}{M^2-s},
\]
\[
\mathcal D=-\int_0^1\frac{dx}{s_+-s_-}\log\left(1+
\frac{(1-x)(s_+-s_-)}{xM^2+(1-x)s_-}\right).
\]

The first uses 75-digit `mpmath`; the second uses adaptive SciPy `quad`. Uniform Simpson integration is **not** used as the required independent validation, because `M^2/s~10^4` produces a sharp near-endpoint boundary layer. UV divergence coefficient cancels by exact `sum_i U_{S2,i} U^*_{S1,i}=0` for each charge block; the numerical code verifies this, rather than discarding divergent terms term-by-term. At fixed source mass parameters, `lambdaT3=0` and `y1=0` give an exactly vanishing one-loop result. The full matrix is complex symmetric; its Takagi masses are the singular values of `m_nu`.

Under nonuniform complex Yukawa deformations, globally rotate `y1 -> exp(i t1)y1`, `y2 -> exp(i t2)y2`, `lambdaT3 -> exp(i tL)lambdaT3`. This chosen convention implies `m_nu -> exp(i(tL-t1-t2))m_nu`. A flavour diagonal neutrino rephasing `nu -> P nu` leads to congruence `m_nu -> P m_nu P`, with accompanying source-Yukawa phases. An intentionally wrong `y1` conjugation must change the complex result. Synthetic three-family engine tests are executable without live Matchete; **nine-model real-source complex tests run in the batch runner on Windows**, not in the static archive snapshot.

## 5. Acceptance, tolerances, comparison

Geometry/charge checks use `5e-11` relative full invariant tensor residual and `1e-12` relative charge/pole checks. The independent loop check tolerates at most `5e-10` relative disagreement; adaptive integration at these benchmark masses achieves much better precision. Complex covariance checks use `2e-11` relative norm due to ordinary double-precision matrix multiplication. Takagi singular values are computed using SVD and checked against the Gram spectrum; the `2e-7` secondary residual permits rounding when a singular value is small compared with the largest. All tolerances are set **before** frozen comparisons.

After independent zero/conditional matrices, integrals, complex controls, and Takagi checks complete, the script may optionally read the frozen provisional matrices for a **comparison-only** Frobenius relative residual. No provisional tensors, loop factors or fitted CG factors are ever used as input. A numerically tiny residual cannot certify the unresolved concrete Matchete `GammaCC` phase, fresh native quadratic-action witness, or the physical benchmark-to-native portal/RGBeta coupling dictionary. `certification_eligible` therefore stays `false` for every model until those gates are independently resolved.

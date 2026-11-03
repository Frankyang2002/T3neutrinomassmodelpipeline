# C_-1 — independent Majorana-triplet direct UV mass (G3)

**Scope:** only the restricted broken-phase one-loop UV self-energy. No EFT, RGBeta computation, running, or production-engine mutation. The live Matchete probe must be run on the user's Windows installation; the result below is **computed here independently using current GitHub-native source data**, NOT a newly executed Matchete witness. Source inspected read-only at `2ac58cc3de12c7af11e2e7e214415dadf0138f4f` (2026-11-03).

## 1. Native source and exact component selection

Model `C_-1`: `dS1=2, dS2=2, dF=3, alpha=-1`; hypercharges `Y(S1)=-1/2`, `Y(S2)=+1/2`, `Y(F)=0`. `F` is self-conjugate in the **real Cartesian adjoint** representation. Native assembled `LFree` has `-MF/2 F^T GammaCC F`; Yukawa originals and Hermitian conjugates each have coefficient +1. Their original fields are `Bar[L] F S1` and `Bar[L] F Bar[S2]` respectively. Native CG source file:

`output/general_observable/lo_direct_one_loop/source_cg_all16/T3_dS1_2_dS2_2_dF_3_alpha_m1/uv_cg_registry.json`.

With 1-based native component labels, the electron-neutrino external leg is `L[1]`, the neutral Majorana triplet Cartesian axis is `F[2]`, the neutral first scalar is `S1[1]`, and the neutral second scalar is `S2[2]`. Let `c=3^(-1/4)`:

| Native coefficient | Exact component | Value |
|---|---|---|
| `T3Y1CG` | `L[1],F[2],S1[1]` | `+c` |
| `T3Y2CG` | `L[1],F[2],Bar[S2[2]]` | `-c` |
| `T3MixCG` | `H[2],H[2],S1[1],Bar[S2[2]]` | `-2/sqrt(3)` |

`H[2]=v/sqrt(2)` induces **only** the neutral scalar mixing `S1[1]`↔`S2[2]`. Charged F± states occur in individual Yukawa vertices, but there is no corresponding Higgs-neutral-induced charged-scalar propagator joining the two different Yukawas. Thus the mixed one-loop neutrino mass has one neutral internal triplet channel, **not a factor 3**.

## 2. Triplet Majorana metric and Weyl mass

The self-conjugate source action has `-MF/2` once. No extra Majorana `2` or `1/2` multiplies the loop self-energy. In the **native real Cartesian gauge basis**, the invariant mass metric is `I3`. For the native-generator-audited weight rotation `F_cart = U3 F_weight` with weights `(+1,0,-1)`, the Weyl symmetric mass metric is instead

\[
G=U_3^T U_3=\begin{pmatrix}0&0&1\\0&1&0\\1&0&0\end{pmatrix}.
\]

This pairs `F+` with `F-` and `F0` with itself. The new reviewer independently reads `matchete_generator_analysis.json` and verifies `U†T3U=diag(1,0,-1)` and `U.T@U=G`. This existing generator report is separate supporting evidence, **not a fresh installed generator measurement**. The loop uses the *native Cartesian neutral* `F[2]`, with mass `MF`, and therefore does not incur a charged-component sign or a triplet multiplicity.

## 3. Scalar quadratic action and full portal identity

`LInt` contains `+lambdaT3 H H S1 Bar[S2] T3MixCG + h.c.`. Since `V=-L` and `H0=v/sqrt(2)`, for the neutral scalars, `Phi=(S1^0,S2^0)` with `V2=Phi† A Phi`,

\[
A=\begin{pmatrix}a&A_{21}^*\\A_{21}&b\end{pmatrix},\qquad
A_{21}=-\frac{\lambda_{T3}v^2}{2}\left(-\frac{2}{\sqrt3}\right)
=\frac{\lambda_{T3}v^2}{\sqrt3}.
\]

For source benchmark `v=246.22 GeV`, `lambdaT3=0.1`: `A21=3500.144922716951 GeV^2`. The live probe must measure this coefficient and the actual action prefactor; it must not reconstruct it from the provisional scalar matrix.

The **exact full four-index native portal identity** for each scalar doublet is

\[
I_{ij}I_{ab}=\tfrac12\mathrm{Inv1}_{ijab}+\tfrac{\sqrt3}{2}\mathrm{Inv2}_{ijab},
\]

not `Inv1` alone. The Higgs-neutral projections are `Inv1=diag(0,1)` and `Inv2=diag(2,1)/sqrt(3)`. For the **selected** benchmark identity-portal assignment, `lambdaH1Inv1=lambdaH1/2`, `lambdaH1Inv2=sqrt(3)*lambdaH1/2`, and analogously for `S2`. This *choice* is not a native assertion that the frozen `lambdaH1,lambdaH2` names directly exist in the Matchete Lagrangian. It gives

\[
a=m_{S1}^2-\lambda_{H1}v^2/2,\qquad
b=m_{S2}^2-\lambda_{H2}v^2/2.
\]

The zero-portal point sets all four native portal couplings to zero. Additional symmetry-allowed scalar couplings are restricted to zero for these benchmark comparisons, not assumed absent by symmetry.

## 4. Spinor convention and independently integrated one-loop mass

In the declared left-handed Weyl convention `L_Y=-nu psi (g1 phi1*+g2 phi2)+h.c.`, source Hermitian conjugation gives

\[
g_{1,ar}=-c\,y_{1,ar}^*,\qquad g_{2,ar}=+c\,y_{2,ar}^*.
\]

The original `+1` Matchete action coefficients and single `-MF/2` are source witnessed; the **absolute installed 4×4 `GammaCC` orientation is still unmeasured**. Coherent external-neutrino signs are unphysical; an incorrectly transported *relative* sign is not. No source physical certification is claimed without the full convention map.

For `s±` eigenvalues of the 2×2 scalar A, evaluate without the production `core.py`/`B0` using

\[
\mathcal D_r=-\frac{1}{s_+-s_-}\int_0^1 dx\,\log\left[1+\frac{(1-x)(s_+-s_-)}{xM_r^2+(1-x)s_-}\right].
\]

Then

\[
(m_\nu)_{ab}=\frac{1}{16\pi^2}\sum_r M_r A_{21}\,\mathcal D_r\big(g_{1,ar}g_{2,br}+g_{2,ar}g_{1,br}\big)
=-\frac{\lambda_{T3}v^2}{3(16\pi^2)}\sum_rM_r\mathcal D_r
\big(y_{1,ar}^*y_{2,br}^*+y_{2,ar}^*y_{1,br}^*\big).
\]

The explicit coefficient `1/3` follows from `c^2=1/sqrt(3)` and `Cmix=-2/sqrt(3)` with `v²/2`. It is **not** fitted to the old amplitude.

### Independent numerical outputs (eV)

Three source points are distinguishable:

- **Native-compatible common MF=100 TeV and zero portals:** Takagi masses, descending: `(0.07029331234,0.04202541753,0.03121193929)` eV. Native one-symbol `MF` describes one common mass; three flavours are retained for Yukawas.
- **Three diagonal benchmark masses (100,110,120 TeV), zero portals:** `(0.06228006048,0.03858577091,0.03749401355)` eV. This is a *benchmark extension*, not the Matchete native one-MF action.
- **Three benchmark masses and selected identity portals:** `(0.062285004868,0.038588872805,0.037497256790)` eV, with independent matrix

\[
(m_\nu/\mathrm{eV})=\begin{pmatrix}
-0.015057960963&-0.030197701558&0.018639953177\\
-0.030197701558&-0.015518566954&-0.037755190720\\
0.018639953177&-0.037755190720&-0.030616860936
\end{pmatrix}.
\]

The provisional **was read after** this independent calculation; its conditional relative Frobenius matrix residual is about `8.19e-12`. It supports numerical compatibility but **does not provide a new physical source certificate**.

## 5. New local experiment and acceptance

Files supplied in one ZIP:

1. `scripts/RunPhysicalCm1NativeScalarQuadratic.wl`: fresh isolated Matchete process, native `LFree`/`LInt`, exact real-triplet Majorana `-MF/2`, charged/neutral CGs, full-index portal identity and scalar mixing. No output overwrite.
2. `scripts/ReviewPhysicalCm1NativeLoop.py`: checks recorded native tensor values and source hashes; validates the independent generator/weight metric; independently integrates the loop; refuses prior outputs; compares frozen only on explicit flag and only *after* computing the result.
3. `Numerical/lo_certification/physical_normalisation/CM1_TRIPLET_MAJORANA_NATIVE_LOOP_GATE.md`: this source-to-loop derivation.

**Pass condition for C_-1 G3 subgate:** native Matchete action and registered tensors pass, neutral-channel coefficient matches the exact analytic contraction, source SHA guard passes, triplet metric check passes, independent numerical self-energy is nonzero and stable. Benchmark source identification, single-vs-three mass distinction, and numerical `GammaCC` phase remain explicit limitations. The physical certification count remains **0/16** until resolved.

# A_+0: independent neutral-component UV one-loop amplitude (2026-11-03)

## Evidence status and scope

- **GitHub read-only source inspected at** `1fb8f4696a5568ec80f5f6dd5c8d6c61780ce40c`, one additive commit after the prior handoff's `6c6321f...`. `T3Topology.wl`, `T3Fields.wl`, `SU2Invariants.wl`, `LagrangianBuilder.wl`, `ScalarPotential.wl`, `T3ModelCatalog.wl`, and `core.py` were inspected. The new commit did not edit these files.
- **Native-source numerical witnesses inspected**: `output/general_observable/lo_direct_one_loop/source_audits/component_support_audit.json` (A_+0 entries); `source_cg_all16/T3_dS1_1_dS2_3_dF_2_alpha_p0/uv_cg_registry.json`; `source_audits/matchete_generator_analysis.json`; `source_portal_cg_all16/T3_dS1_1_dS2_3_dF_2_alpha_p0/native_higgs_portal_registry.json`; `output/general_observable/physical_normalisation/b1_action_contract_20261103_103545/native_action_coefficients.json` (native action coefficient witness). These were generated previously on the user's installation, not generated afresh in this deliverable.
- **Independently derived and evaluated in this analysis environment**: the A_+0 charge-selected, two-complex-scalar propagator integral, with a Feynman-parameter integral evaluated at high precision. This does not import `Numerical/direct_one_loop/core.py`, `provisional_inputs`, or the frozen result. The supplied script was smoke-tested locally against separately transcribed source numbers, **not run on the user's Windows Matchete installation**. When run locally it loads the actual JSON source files.
- **Not claimed**: actual native portal coupling identification of the frozen benchmark, installed numerical orientation of Matchete's abstract `GammaCC`, an absolute cross-package CG map to RGBeta, or certification of any of 16 physical models. No EFT or running.

## 1. Exact component and charge selection

`A_+0` is `(dS1,dS2,dF)=(1,3,2)`, with `(YS1,YS2,YF)=(0,1,1/2)` using `Q=T3+Y`. The relevant `SU(2)` triplet in the actual source Cartesian basis has

\[
T_3=\begin{pmatrix}0&0&-i\\0&0&0\\i&0&0\end{pmatrix},\quad
u\text{ is } L_1,\quad F^0=F_2,\quad S_2^0=S_{2,-1},\quad
u_-=\frac{1}{\sqrt2}(i,0,1)^T.
\]

Here `S_{2,cart}=U_3 S_{2,weight}` and `U_3` is taken from the **native generator witness**, not from the provisional input builder. Because the source operators contain `Bar[S2]`, project their CG coefficients with `conj(u_-)`. Define `c=3^{-1/4}`. The **native CG component audit** contains:

\[
\begin{array}{c|c|c}
\text{vertex}&\text{original Cartesian component coefficients}&\text{neutral projection}\\\hline
T3Y1CG&(L_1,F_2,S_1):1& C_1=1\\
T3Y2CG&(L_1,F_2,\bar S_{2,1}):c,\ (L_1,F_2,\bar S_{2,3}):-ic &C_2=-i\sqrt2\,c\\
T3MixCG&(H_2,H_2,S_1,\bar S_{2,1}):-c,\ (H_2,H_2,S_1,\bar S_{2,3}):+ic &C_m=+i\sqrt2\,c
\end{array}
\]

Thus `C2=-Cm`. This is an **explicit physical electric-charge channel contraction**; the six nonzero general tensor entries are **not** a loop multiplicity. For A_+0 only the common `F^0`, `S1^0`, and `S2^0` line contributes to the `y1*y2` neutrino mass.

## 2. Four-spinor to Weyl vertex coefficients

Native A_+0 action coefficients previously measured: `-1 MF Bar[F] F`; `+1 y1` and `+1 y2` for their original terms **and** H.c.; `y1` uses `Bar[L] CConj[F] S1`, `y2` uses `Bar[L] F Bar[S2]`. In the declared Weyl basis `F=(chi,bar(xi))` and `C=+i gamma^2 gamma^0`, the H.c. forms are `+y1* nu chi S1*` and `+y2* nu xi S2` times conjugated native CGs. For the bookkeeping convention `L_Y=-nu chi (g1 S1*)-nu xi(g2 S2)+h.c.`,

\[
g_{1,ar}=-y_{1,ar}^{*}C_1^*=-y_{1,ar}^*,\quad
g_{2,ar}=-y_{2,ar}^*C_2^*=-i\sqrt2\,c\,y_{2,ar}^*,\quad
\mathcal L_M=-M_r\chi_r\xi_r+\mathrm{h.c.}
\]

No artificial `1/2` or `2` is included for a Dirac mass. The chosen `C` convention's overall phase remains unanchored to a literal installed 4x4 numerical Matchete matrix, but an overall coherent field rephasing changes no physical masses. For this real benchmark the resulting matrix is real.

## 3. Physical neutral scalar mass mixing

The source operator is **plus** `lambdaT3 H H S1 Bar[S2] CG + h.c.` with **no explicit prefactor 1/2**. Set `H_2=v/sqrt(2)` and write the complex charge-neutral scalar column `Phi=(S1,S2^0)^T`, with `V2=Phi^dag A Phi`:

\[
A=\begin{pmatrix}a&A_{21}^{*}\\ A_{21}&b\end{pmatrix},\qquad
A_{21}=-\frac{\lambda_{T3}v^2}{2}C_m
=-i\frac{\lambda_{T3}v^2}{\sqrt2\,3^{1/4}}.
\]

`a` and `b` are physical diagonal squared masses, including any **explicitly chosen** Higgs portal assignments; the native bare scalar mass parameters are `MS1, MS2`. Diagonalising gives

\[
s_\pm=\frac{a+b}{2}\pm\sqrt{\left(\frac{a-b}{2}\right)^2+|A_{21}|^2}.
\]

Two reproducible points are evaluated:

1. **Explicit zero-portal source point**: `a=mS1Sq`, `b=mS2Sq`, and all other optional quadratic/portal interactions are set to zero. This does **not** reproduce the frozen numerical benchmark's nonzero `lambdaH1,lambdaH2`.
2. **Source-compatible identity-portal assignment**: take `lambdaH1Inv1=lambdaH1=0.02`, `lambdaH2Inv2=lambdaH2=0.03`, `lambdaH2Inv1=0`. The GitHub native portal registry contains the singlet identity tensor for `lambdaH1Inv1` and the `delta(Hbar,H) delta(S2bar,S2)` tensor for `lambdaH2Inv2`. Under this **explicit chosen coupling identification**, `a=mS1Sq-0.02*v^2/2` and `b=mS2Sq-0.03*v^2/2`. The frozen benchmark uses the labels `lambdaH1,lambdaH2`, whereas the actual builder registers `lambdaH1Inv1, lambdaH2Inv1, lambdaH2Inv2`; that identity **mapping is not automatically established** by the benchmark JSON. Accordingly, the script labels this evaluation conditional.

For both points the charged triplet components have the diagonal `b` mass in this selected identity-portal subspace; they do **not** join the y1 neutral loop. A neutral singlet holomorphic quadratic, permitted in broader model variants, is set to zero by the restricted point.

## 4. Independently derived loop integral and absolute coefficient

At zero external momentum, the chirality-changing propagator is proportional to `M_r`, and the scalar correlator is `<S1* S2>`. For any analytic function `f`, a 2x2 Hermitian matrix obeys

\[
[f(A)]_{21}=A_{21}\frac{f(s_+)-f(s_-)}{s_+-s_-}.
\]

Here the **finite difference** of two scalar one-loop functions is ultraviolet finite. Rather than importing `core.py`'s `B0`, calculate its divided difference directly as a Feynman-parameter integral:

\[
\mathcal D(M^2;s_-,s_+) =
-\frac{1}{s_+-s_-}\int_0^1 dx\,
\ln\!\left[1+\frac{(1-x)(s_+-s_-)}{xM^2+(1-x)s_-}\right].
\]

The physical (chosen-Weyl-convention) Majorana self-energy is

\[
\boxed{(m_\nu)_{ab}=\frac1{16\pi^2}\sum_{r=1}^3
M_r A_{21}\mathcal D(M_r^2;s_-,s_+)
\left(g_{1,ar}g_{2,br}+g_{2,ar}g_{1,br}\right).}
\]

Inserting **actual native neutral CG values**, for real positive diagonal `M_r`,

\[
\boxed{(m_\nu)_{ab}=\frac{\lambda_{T3}v^2}{16\pi^2\sqrt3}
\sum_r M_r\mathcal D_r
\left(y_{1,ar}^*y_{2,br}^*+y_{2,ar}^*y_{1,br}^*\right).}
\]

This explicitly traces the `1/sqrt(3)` coefficient to **the product of native CG tensors**, the `v^2/2` scalar insertion, and the complex triplet normalization. No fitted multiplicity, extra `2`, or guessed CG factor is introduced. Since `D<0` for positive physical masses, the mass is proportional to the indicated signed source product, up to the stated Weyl phase convention. For a complex `lambdaT3`, the matrix transforms consistently with the conjugate source spurion invariant `y1* y2* M lambdaT3`.

As an **independent check of propagator algebra**, set `|A21|^2=0` only inside scalar eigenvalues (`s_\pm -> a,b`), retaining one off-diagonal insertion. The difference from the exact two-scalar result is order `|A21|^2/(a-b)^2` and measures about `1.62e-7` of the full matrix norm for this benchmark. Turning off `lambdaT3` makes the entire amplitude zero.

## 5. Numerical results, obtained independently before frozen comparison

Source couplings from `Numerical/benchmarks/fitted_reference/T3_dS1_1_dS2_3_dF_2_alpha_p0.json`: `v=246.22 GeV`, `M_r=(100,110,120) TeV`, `lambdaT3=0.1`, `mS1Sq=1306122.448979592 GeV²`, `mS2Sq=1000000 GeV²`, and the nine real `y1`, `y2` entries in the benchmark. High-precision (`mpmath`, 55 decimal digits) quadrature implements the equation above.

**Zero portals, independently defined source point:**

\[
A_{21}=-3257.251872642\,i\;\mathrm{GeV}^2,\quad
(s_-,s_+)\simeq(999965.345603,1306157.103376)\;\mathrm{GeV}^2,
\]

\[
m_\nu/\mathrm{eV}=\begin{pmatrix}
0.002634548314&-0.001301923991&0.010752697408\\
-0.001301923991&0.020721197519&0.026437299468\\
0.010752697408&0.026437299468&0.021353870032
\end{pmatrix},
\]

Takagi singular masses `(0.00721371755,0.01099915442,0.04849505273) eV` (ascending).

**Identity portals, conditional matching to frozen benchmark:**

\[
(a,b)=(1305516.206095592,999090.635674)\;\mathrm{GeV}^2,
\]

\[
m_\nu/\mathrm{eV}=\begin{pmatrix}
0.002634773571523&-0.001301982231411&0.010753592992369\\
-0.001301982231411&0.020722820246469&0.026439424399483\\
0.010753592992369&0.026439424399483&0.021355599196051
\end{pmatrix}.
\]

Only **after** obtaining these two independent matrices was the existing frozen provisional `output/general_observable/lo_direct_one_loop/provisional/A_+0_direct_one_loop.json` inspected. Its same-point matrix has relative Frobenius discrepancy `3.05e-12` in this evaluation, arising from floating-point cancellation in the provisional eigenstate `B0` route. This agreement is **conditional corroboration, not an independent Matchete portal-action certificate**. The independent physical contribution is the explicit neutral-channel analytical coefficient and Feynman-parameter matrix above.

## 6. Missing decisive native witness

To complete physical source certification for the **nonzero-portal frozen benchmark**, perform a **native broken-phase scalar quadratic action extraction** for A_+0, not another tensor replay. It must directly witness:

- `H_2 -> v/sqrt(2)` in the assembled A_+0 `LInt` and the physical Cartesian-to-charge rotation of S2;
- the `S1 Bar[S2^0]` coefficient `+i sqrt(2) 3^(-1/4) lambdaT3 v²/2` in the **Lagrangian**, so `A21` has the minus sign shown above;
- `lambdaH1Inv1` and `lambdaH2Inv2` diagonal contributions `+lambda v²/2` in `LInt`, and `lambdaH2Inv1=0`, with source coupling identities explicitly recorded;
- whether any unintended quadratic/source-normalisation term was implicitly added by `FreeLag` or default candidate selection.

The attached script computes the physical charge-selected amplitude **without reusing frozen tensors** and reports `physical_source_certified=false` until that local witness and identification are independently established. It never calls `Match`, overwrites older outputs, modifies GitHub, or performs EFT.

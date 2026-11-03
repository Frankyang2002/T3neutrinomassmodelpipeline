# B_-1: independent Majorana-singlet physical UV mass gate (2026-11-03)

## Evidence taxonomy, scope and status

- **Read-only source inspected** at GitHub commit `f82ee42ea2456e07e831a381fcbf2f2a3dd747d5`:
  `Lagrangian/model/T3ModelCatalog.wl`, `Lagrangian/model/T3Fields.wl`,
  `Lagrangian/interactions/T3Topology.wl`, `Lagrangian/interactions/ScalarPotential.wl`,
  `output/general_observable/physical_normalisation/b1_02/native/B_m1.json`,
  `output/general_observable/lo_direct_one_loop/source_cg_all16/T3_dS1_2_dS2_2_dF_1_alpha_m1/uv_cg_registry.json`,
  `output/general_observable/lo_direct_one_loop/source_portal_cg_all16/T3_dS1_2_dS2_2_dF_1_alpha_m1/native_higgs_portal_registry.json`,
  and `Numerical/benchmarks/fitted_reference/T3_dS1_2_dS2_2_dF_1_alpha_m1.json`.
- **Native B_-1 mass/Yukawa action already observed on user's installation:** the cached `B_m1.json` has the `-MF/2` Majorana mass, original and H.c. Yukawa `+1`, and the same heavy four-spinor in both vertices. The cached source CGs (not provisional tensor inputs) specify the contractions below.
- **Independent calculation already executed in the analysis environment:** direct, 60-digit Feynman-parameter integral with explicit 2-by-2 neutral scalar mixing and 3-by-3 real Yukawa entries. The **new local Wolfram Matchete script has not been executed** here; its Python companion was smoke-tested against a deliberately synthetic, clearly labelled source-shaped fixture, not counted as native physical evidence.
- **No physical source-certified model claimed.** The default Matchete `MF` is one symbol, while frozen fitted benchmarks contain 3 distinct flavour masses. The concrete source and the optional diagonal-flavour extension are reported separately. Installed numerical `GammaCC` is not an absolute matrix sign witness. Zero portal and frozen identity portal are distinct points. No EFT/running.

## 1. Native B_-1 component contractions

For `(dS1,dS2,dF,alpha)=(2,2,1,-1)` the hypercharges are
`Y(S1)=-1/2`, `Y(S2)=+1/2`, `Y(F)=0`; the heavy fermion is a Majorana singlet.

With `Q=T3+Y` and the native (weight-aligned) SU(2) doublet axes:

- `nu=L[1]`, neutral `S1=S1[1]`, neutral `S2=S2[2]`;
- `T3Y1CG[nu,F,S1] = T3Y1CG[1,1] = +1`;
- `T3Y2CG[nu,F,bar(S2)] = T3Y2CG[1,2] = +1`;
- `T3MixCG[H0,H0,S1,bar(S2)] = T3MixCG[2,2,1,2] = -2/sqrt(3)`.

There is **one neutral scalar charge channel**, not two factors of the same SU(2) channel. The original native scalar mixing term has `+lambdaT3` without `1/2`; after inserting *both* Higgs vevs `H[2]=v/sqrt(2)`:

\[
 {\cal L}_{\mathrm{mix}}=-\frac{\lambda_{T3}v^2}{\sqrt3}S_1^0 S_2^{0*}+\mathrm{h.c.}
\]

For `V_2=Phi^dag A Phi`, `Phi=(S1^0,S2^0)`, the off-diagonal mass entry is

\[
 A_{21}=+\frac{\lambda_{T3}v^2}{\sqrt3}.
\]

At the frozen `v=246.22 GeV`, `lambdaT3=0.1`, this is **`+3500.144922716951 GeV^2`**.

The native Higgs-portal doublet tensors `lambdaH1Inv1`, `lambdaH2Inv1` project to `I_2` at `H0` (with their `1/2 PlusHc` source factors correctly combined). Their `Inv2` companions project to `diag(2/sqrt(3),1/sqrt(3))` and are independently allowed, but set to zero in the selected identity-portal benchmark subspace. Thus:

\[
 a=m_{S1}^2-\lambda_{H1Inv1}\frac{v^2}{2},\qquad
 b=m_{S2}^2-\lambda_{H2Inv1}\frac{v^2}{2}.
\]

The source **does not mandate** the benchmark dictionary `lambdaH1Inv1=lambdaH1`, `lambdaH2Inv1=lambdaH2`; that is an explicit identity-portal identification for the benchmark comparison.

## 2. Majorana mass and two-component vertices

The original Matchete action is

\[
 {\cal L}_F=-\frac{M_F}{2}\bar F F,\qquad
 +y_1\bar L F S_1 +y_2\bar L F S_2^\dagger+\mathrm{h.c.},
\]

with the source-specific CG coefficients in section 1. In the explicitly chosen convention `F=(chi,bar(chi))`, the corresponding left-handed vertex convention is

\[
 {\cal L}_Y=-\nu_a\chi_r\left[g_{1,ar}S_1^{0*}+g_{2,ar}S_2^0\right]+\mathrm{h.c.},
\qquad g_{1,ar}=-y_{1,ar}^*,\quad g_{2,ar}=-y_{2,ar}^*.
\]

The overall Weyl charge-conjugation basis phase remains conventional and was not numerically fixed by the installed abstract `GammaCC` symbol. **Do not insert a second `1/2` in the self-energy**: the canonical Majorana `1/2` belongs to the free action only, while the cross term contains `g1*g2 + g2*g1` once.

## 3. Independent physical loop calculation

For `s_- < s_+` the eigenvalues of the physical neutral scalar squared-mass matrix are

\[
 s_{\pm}=\frac{a+b}{2}\pm\sqrt{\frac{(a-b)^2}{4}+|A_{21}|^2}.
\]

Define, directly with a finite Feynman-parameter integral (no production `B0` calls):

\[
 {\cal D}(M^2;s_-,s_+)=
 -\int_0^1 dx\;\frac{1}{s_+-s_-}
 \ln\left(1+\frac{(1-x)(s_+-s_-)}{x M^2+(1-x)s_-}\right).
\]

The one-loop mass is

\[
 \boxed{
 (m_\nu)_{ab}=\frac{1}{16\pi^2}\sum_r M_r A_{21}{\cal D}_r
 (g_{1,ar}g_{2,br}+g_{2,ar}g_{1,br})
 }
\]

or, for positive real Majorana masses in the chosen Weyl convention,

\[
 \boxed{
 (m_\nu)_{ab}=\frac{\lambda_{T3}v^2}{16\pi^2\sqrt3}
 \sum_r M_r {\cal D}_r
 (y_{1,ar}^*y_{2,br}^*+y_{2,ar}^*y_{1,br}^*) .
 }
\]

The antisymmetric external Majorana spinor indices and symmetric flavour matrix carry the standard mass-operator combinatorics. **No extra `2` or `1/2`** multiplies the result. The scalar mixing and Yukawa product are source-derived; the remaining unmeasured absolute `GammaCC` phase affects only the coherent convention for the selected pure topology, not the Takagi spectrum.

## 4. Independent numerical results

Frozen Yukawas from the B_-1 `fitted_reference` JSON, with `lambdaT3=0.1`, `mS1Sq=1306122.448979592 GeV^2`, `mS2Sq=1000000 GeV^2`. Computed with `mpmath` high-precision quadrature before any provisional matrix was read.

**Native-compatible common `MF=100 TeV` and zero portals:**

\[
 m_\nu / \mathrm{eV} =
 \begin{pmatrix}
 -0.002357796512&-0.009505835126&+0.000326069655\\
 -0.009505835126&-0.019684386584&+0.026683254267\\
 +0.000326069655&+0.026683254267&-0.024340534475
 \end{pmatrix}.
\]

Takagi masses in descending order:
`[0.0497509127555, 0.00908669075363, 0.00571849556947] eV`.

**Frozen three-mass extension** `MF=(100,110,120) TeV` with zero portals:

\[
 m_\nu / \mathrm{eV} =
 \begin{pmatrix}
 -0.001231168917&-0.009877517063&+0.000945517892\\
 -0.009877517063&-0.019332801247&+0.027234618386\\
 +0.000945517892&+0.027234618386&-0.024592751592
 \end{pmatrix}.
\]

Largest Takagi mass approximately `0.05043713 eV`.

**Frozen three-mass extension and chosen identity portals:**

\[
 m_\nu / \mathrm{eV} =
 \begin{pmatrix}
 -0.001231240637&-0.009878342257&+0.000945610552\\
 -0.009878342257&-0.019334391099&+0.027236878856\\
 +0.000945610552&+0.027236878856&-0.024594786699
 \end{pmatrix}.
\]

The *later* comparison to `output/general_observable/lo_direct_one_loop/provisional/B_-1_direct_one_loop.json` gives relative Frobenius residual ~`1.84e-11`. It is explicitly **conditional**, not evidence that the input builder independently certifies its own conventions.

## 5. Next native physics gate and installation

Install these **new full files** only:

- `scripts/RunPhysicalBm1NativeScalarQuadratic.wl` (live B_-1 Majorana, portal and mixing extraction)
- `scripts/ReviewPhysicalBm1NativeLoop.py` (independent quadrature from the LIVE witness)
- `Numerical/lo_certification/physical_normalisation/BM1_NATIVE_MAJORANA_LOOP_NEXT_GATE.md` (this derivation)

Run from the repository root, choose a fresh evidence directory, and first execute native Wolfram extraction. The reviewer refuses to proceed without the live witness and checks the exact source hashes against the local native extraction. It then calculates a native-compatible common-mass point and a separately labelled frozen-extension comparison. The live Wolfram script has **not been run in this analysis environment**.

### Certification boundary

The native extraction is the next decisive **G3 Majorana subgate**, not automatic complete certification. Its *pass* establishes native scalar mass, charged channel, Yukawa and Majorana half-factor for B_-1. A fail identifies a specific incorrect physical assumption. Afterwards, the explicit local source `MF` generality and gamma/phase interpretation must be documented before any `1/16` claim. `C_-1` still requires a separate triplet mass metric, not inferred from this singlet calculation.

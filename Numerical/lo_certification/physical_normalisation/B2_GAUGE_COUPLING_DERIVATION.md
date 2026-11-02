# B2 — SU(2) gauge components and cross-package coupling dictionary

**Scope:** restricted direct-UV T3 topology only. **Source inspected:** public GitHub `main`, commit `0b210f18f6a2cfb10fff03059e2ec4a73d2195a4` (3 November 2026, Sydney). Local Windows Git revision and installed Mathematica/RGBeta/Matchete cannot be observed in this environment. New scripts validate local source hashes rather than assuming equivalence with this reference. All absolute physical Matchete↔RGBeta coupling ratios and `physical_crosspackage_complete` remain **unresolved/false**.

## 1. Actual executable definitions (verified from reference source)

| Coupling | Matchete field monomial within `PlusHc` | RGBeta registered ordered fields and invariant |
|---|---|---|
| `y1` | `y1 NCM[Bar[l],CConj[NewFermion]] S1 CG[T3Y1CG]` | `{Bar[S1],l,F}`; `eps` or `delS2` |
| `y2` | `y2 NCM[Bar[l],NewFermion] Bar[S2] CG[T3Y2CG]` | `{S2,l,Fc}` for Dirac; `{S2,l,F}` for Majorana; `eps` or `delS2` |
| `MF` | mass automatically from `DefineField[NewFermion, ..., Mass -> MF]` | `{F,Fc}` Dirac; `{F,F}` Majorana; doublet Dirac `eps`, triplet `del`, singlet `1` |
| `lambdaT3` | `PlusHc[lambdaT3 H H S1 Bar[S2] CG[T3MixCG]]`, **no 1/2** | `{H,H,S1,Bar[S2]}`, `SelfConjugate->False`, `T3RGBetaMixInvariant` |
| `lambdaH1`, `lambdaH2` | real-coupling `1/2 PlusHc[lambdaHiInvN Hbar H Sbar S CG]` | `{Bar[H],H,Bar[Si],Si}`, `delta(Hbar,H) delta(Sbar,S)` |

Matchete `DefineSU2InvariantCG` calls `InvariantTensors[SU[2],algebraReps]`, then `DefineCG[cgName,cgReps,First[tensors]]`. Its field flags, **not** naive tensor indices, fix bar orientations:

- `T3Y1CG`: dimensions `(2,dF,dS1)`, conjugation `(True,True,False)`.
- `T3Y2CG`: `(2,dF,dS2)`, `(True,False,True)`.
- `T3MixCG`: `(2,2,dS1,dS2)`, `(False,False,False,True)`, identical first two Higgs fields symmetric.

RGBeta `T3RGBetaYukawaInvariant[dS,dF]` defines `(2,1)` as `eps[SU2L[fund],#1,#2]`, `(1,2)` as `eps[SU2L[fund],#2,#3]`, `(2,3)` as `delS2[SU2L,#3,#1,#2]`, `(3,2)` as `delS2[SU2L,#1,#2,#3]`. These are symbolic registry expressions **not numerical components** in any declared coordinate system. `T3RGBetaMixInvariant` has cases `(1,3)`, `(3,1)`, `(2,2)`, `(3,3)`, with the `(3,3)` expression using three `delS2` and three `eps` contractions. Their fully normalised values and index orientations are **pending** direct package evidence. The basic RGBeta quartics explicitly use a singlet `delta * delta` with no additional scalar-index multiplier.

RGBeta internally uses `sqrt(2)` in Yukawa tensors and a nontrivial internal MF multiplier reported by the prior representative native registry probes. **Neither is automatically a physical action multiplier.** In particular, one must reconstruct the real-scalar decomposition, the permutation factors of two Weyl fermions, and the overall sign of the RGBeta action before equating its symbols to Matchete symbols. B1 must resolve `GammaCC` and the 4-spinor charge-conjugation map; introducing a guessed `eta = ±1` as a fact is prohibited.

## 2. Common basis and mass metrics

Let the weight-basis component vector be `phi_weight` and the Matchete original vector be `phi_M = U_d phi_weight`. The concrete diagonalizers recorded by the generator audit are

\[
U_1=(1),\quad U_2=\begin{pmatrix}1&0\\0&1\end{pmatrix},\quad
U_3=\begin{pmatrix}-i/\sqrt2&0&i/\sqrt2\\0&1&0\\1/\sqrt2&0&1/\sqrt2\end{pmatrix}.
\]

Then `U_d^† T3_M U_d = diag(j,j-1,...,-j)`, `j=(d-1)/2`. Conjugated field slots transport with `U_d*`, unbarred with `U_d`. This is what `native_tensor()` does for the three known Matchete CG orientations; a numerical RGBeta component basis has **not** been inferred from names like `fund`, `S2` or `delS2`. Weight-basis phase changes `U_d -> U_d diag(e^{i theta_m})` are allowed when all coupling tensors transform coherently.

For two unbarred fundamental doublets, a chosen invariant is

\[
\epsilon = \begin{pmatrix}0&1\\-1&0\end{pmatrix},\qquad
T_a^T\epsilon+\epsilon T_a=0.
\]

Replacing it with the identity fails the Ward identities. The orientation `epsilon_12=+1` is a **convention**, not an independently exported RGBeta value. For Cartesian triplets use `(T_a)_{bc}=-i epsilon_{abc}` and `g_3=U_3^T U_3`:

\[
g_3=\begin{pmatrix}0&0&1\\0&1&0\\1&0&0\end{pmatrix},\qquad T_a^T g_3+g_3 T_a=0.
\]

This is the mass metric associated with Matchete's reported Cartesian-to-weight transform, **not an assertion** that RGBeta `del[SU2L[S2],...]` has these numerical components. A weight rephasing changes `g_3` by congruence and is permissible.

## 3. Charge sectors and alpha transport

From RGBeta's registered representations, `Y_F=(alpha+1)/2`, `Y_S1=alpha/2`, `Y_S2=(alpha+2)/2`, with electric charge `Q=m+Y` and weight `m=j-index` (index zero-based). The necessary component support in the neutral-external-lepton topology is `Q_F=Q_S1=Q_S2`. Its charge-compatible sectors for all 16 frozen representations are:

| model | possible common component charge `Q` |
|---|---|
| A -4, -2, 0 | -2; -1; 0 |
| B -3, -1, +1 | -1; 0; +1 |
| C -3, -1, +1 | -1; 0; +1 |
| D -2, 0, +2 | 0; +1; +2 |
| E -4 | -2, -1 |
| E -2 | -1, 0 |
| E 0 | 0, +1 |
| E +2 | +1, +2 |

These are necessary *charge skeletons*, not counted nonzero physical loop channels. A channel survives only if **all three oriented native CGs** and the heavy-mass contraction are nonzero. The shipped audit consumes the existing native component-support and generator JSON on the user's local machine if and only if source bytes and audit formats pass strict guards; it emits the labelled diagnostic `conj(C_y1) * conj(C_y2) * C_mix`. This diagnostic **is not the full loop amplitude**, and cannot be treated as a signed physical dictionary. Dependence on alpha arises in selecting which charges are present; reduced SU(2) tensor shapes depend on the representation class. Majorana F requires a charge-pairing mass contraction and may not be reduced to a Dirac field-label substitution.

## 4. Spurion invariant — explicitly *conditional* on B1 Weyl action

Assume the candidate left-handed action

\[
\mathcal{L}\supset y_1\nu F^c S_1+y_2\nu F S_2^* + M_F FF^c +\lambda_{T3} HH S_1S_2^*+\mathrm{h.c.}
\]

and field redefinitions `(F,Fc,S1,S2) -> (e^{ia}F,e^{ib}Fc,e^{ic}S1,e^{id}S2)`. Couplings change by

\[
y_1\mapsto e^{-i(b+c)}y_1,\quad y_2\mapsto e^{-i(a-d)}y_2,\quad
M_F\mapsto e^{-i(a+b)}M_F,\quad\lambda_{T3}\mapsto e^{-i(c-d)}\lambda_{T3}.
\]

Thus the rephasing-neutral candidate is

\[
\boxed{\mathcal I=y_1y_2 M_F^*\lambda_{T3}^*.}
\]

For Majorana `Fc=F` rephase with `b=a`, and `M_F` transforms with `e^{-2ia}`; the same spurion is neutral. This proves only the invariance of the **declared candidate monomials**, not the actual Matchete/ RGBeta sign, field ordering, operator conjugation, or mass-symmetry prefactors. The complex phase of `I` must not be changed by an arbitrary field phase or by choosing a sign for `GammaCC` without consistent propagation.

## 5. Mixing and Higgs-portal operators

Under the source builder's `H^0=v/sqrt(2)` and `V=-L` conventions, the native Matchete scalar-Hermitian block is

\[
(A)_{S_2,S_1}=-\frac{v^2}{2}\lambda_{T3}C^{(M)}_{H^0 H^0 S_1\bar S_2},\qquad
(A)_{S_1,S_2}=(A)_{S_2,S_1}^*.
\]

This is a source-to-builder consistency formula, not the desired physical **RGBeta-to-Matchete** `lambdaT3` map. Even an exact 16/16 agreement with frozen `A` would not certify the latter. Native RGBeta `MixInvariant` numerical components and physical quartic action prefactors are required.

For basic portals define the native Hermitian operator tensor

\[
P_k(\bar h,h,\bar s,s)=\tfrac12\left[C_k(\bar h,h,\bar s,s)+C_k( h,\bar h,s,\bar s)^*\right].
\]

Solve **over real coupling coefficients** `sum_k c_k P_k = delta_H delta_S`; the complete solutions are `c = c_pivot + sum_a t_a n_a`, with each `n_a` in the real nullspace. Existing source-linked portal evidence covers 10 sectors (five classes times two scalars), but A/S2, D/S1, E/S1, E/S2 are nonunique; they contain an **entire Hermitian null direction**, not an arbitrary choice of a uniquely fixed coupling. The committed previous native-portal dictionary (Git blob `ce6d2ad3e36bf3b82f018330355c84f5b5ba7acf`) gives these **prior** operator embeddings; they are not newly rerun native values:

| sector | one canonical solution `c_pivot` | real Hermitian nullspace | status |
|---|---|---|---|
| A/S1 | `(1)` | none | unique |
| A/S2 | `(0,1)` | `span{(1,0)}` | nonunique |
| B/S1, B/S2 | `(1/2, sqrt(3)/2)` | none | unique |
| C/S1, C/S2 | `(1/2, sqrt(3)/2)` | none | unique |
| D/S1 | `(0,1)` | `span{(1,0)}` | nonunique |
| D/S2 | `(1)` | none | unique |
| E/S1, E/S2 | `(0,1)` | `span{(1,0)}` | nonunique |

Thus for the four nonunique sectors the most general allowed **operator-equivalent** real coupling vector is `c=(t,1)` for arbitrary real `t`. This null direction does **not** represent a freely adjustable physical effect in the basic Hermitian portal operator: it multiplies a vanishing Hermitian tensor. Physical RGBeta action normalisation is still not certified.
 `prior_portal_dictionary()` includes the full exact vectors only when their previous proof file and every recorded source SHA256 match the local files. Otherwise it blocks rather than inventing nullspace coefficients. Native physical RGBeta action normalisation remains unproved even after a successful operator embedding.

## 6. Status matrix and next evidence

| Task | This batch |
|---|---|
| Read reference source operators and their field order | **Verified from GitHub commit** |
| Source-byte guarded all-16 charge skeleton | **Proved algebraically**; 16/16 |
| Doublet mass epsilon / triplet metric Ward checks | **Passed synthetic exact tests**; chosen basis |
| Matchete CG rotation and diagnostic products | **Conditional on local existing native audits** |
| Previous basic portal nullspaces | **Import only with SHA256 source match**; no unique map in four sectors |
| Installed native RGBeta numeric invariant axes | **Not run/not established here** |
| Full physical RGBeta action signs/factors | **Not established** |
| B1 four-spinor↔Weyl action phase | **Pending separate B1** |
| Physical `y1`, `y2`, `MF`, `lambdaT3` source-to-source ratios | **Unknown**; cannot safely integrate |
| Direct one-loop eigenmasses, RGE, EFT | **Not done by B2** |

**Next minimal evidence:** fresh-kernel RGBeta values (and the real action coefficient) for the four invariant operators in A `(1,3,2)`, B `(2,2,1)`, C `(2,2,3)`, D `(3,1,2)`, E `(3,3,2)` in defined index bases, including numeric `eps`, `del`, `delS2` normalisations. The supplied Wolfram diagnostic records registered expressions/definitions, but intentionally does not claim those numerical values; it should be paired with explicit package tensor component data or a native Ward/action measurement. Then solve the joint complex tensor equations, carrying one common B1 spinor phase `eta` and allowing consistent field rephasings. Until then **`physical_crosspackage_complete=false`, `physical_source_certified=false`**, for all 16.

No existing production code, frozen inputs, CG audits, or GitHub state is modified. No running or EFT is part of this worker. Any later EFT must follow UV → integrate out all BSM fermions together → EFT1 with both BSM scalars → integrate out both BSM scalars together → SMEFT.

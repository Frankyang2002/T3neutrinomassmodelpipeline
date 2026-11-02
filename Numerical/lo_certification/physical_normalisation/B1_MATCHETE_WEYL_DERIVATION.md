# B1 — restricted Matchete four-spinor to Weyl dictionary (2026-11-03)

## Evidence classification and scope

* **GitHub source inspected, read-only**: `main` = `0b210f18f6a2cfb10fff03059e2ec4a73d2195a4`, confirmed by commit lookup and compare. Source file blobs are pinned by `scripts/RunPhysicalB1WeylAudit.py`.
* **Executed in this deliverable's Python environment**: exact exterior-algebra 4x4 identities and independent synthetic tests. These verify algebra **only for the chosen, declared gamma-matrix basis**.
* **Not executed here**: installed Matchete, RGBeta, and the sixteen local benchmark/covariance JSONs. The user's `C:\thesistemp\T3neutrinomassmodelpipeline` working tree is not mounted; `wolframscript` is unavailable. Native action observation therefore remains `false` in this handoff.
* **Physical certificate**: `false`. The shared input builder is provisional, and neither its own replay nor this conditional agreement supplies an installed Matchete ↔ RGBeta SU(2) and charge-conjugation dictionary.

Only the frozen **restricted** topology is intended, with unspecified allowed interactions set to zero. This batch performs neither running nor EFT matching nor numerical result production. Three representatives are A_+0 (Dirac), B_-1 (Majorana singlet), C_-1 (Majorana triplet). There are fourteen Dirac and two Majorana cases across the catalogue.

## 1. Source-executable contracts (GitHub main)

`Lagrangian/interactions/T3Topology.wl` defines

```mathematica
PlusHc[y1[p,r] NCM[Bar[l[i,p]], CConj[NewFermion[...]]] S1 CG[T3Y1CG,...]]
PlusHc[y2[p,r] NCM[Bar[l[i,p]], NewFermion[...]] Bar[S2] CG[T3Y2CG,...]]
PlusHc[lambdaT3[] H[i] H[j] S1 Bar[S2] CG[T3MixCG,...]]
```

The singlet representations suppress their gauge index. `BuildT3Lagrangian` takes `LFree = FreeLag[NewFermion, NewScalar1, NewScalar2] // RelabelIndices`, assembles `LInt`, `LUV`, validates with `CheckLagrangian`, and returns all three in an association. `T3Fields.wl` uses `SelfConjugate -> True` iff the fermion hypercharge vanishes and its SU(2) dimension is odd. `T3MassSpec` defaults to `{Heavy, MF}`; its optional flavour-mass mode describes diagonal eigenvalues, **not an arbitrary off-diagonal Matchete mass**.

`RGE/running/rgbeta/T3RGBetaModel.wl` has `AddYukawa[y1,{Bar@S1,l,F},...]`, `AddYukawa[y2,{S2,l,Fc},...]` for Dirac (or `{S2,l,F}` for Majorana), `AddFermionMass[MF,{F,Fc},...]` for Dirac (or `{F,F}` for Majorana). The native Dirac doublet invariant is `eps` and the Majorana triplet invariant `del`. These are **field and symbolic tensor contracts**, not a measured coefficient/axis mapping. The RGBeta registry witness also reports internal `Sqrt[2]` Yukawa prefactors; that registry factor alone is **not** the physical Lagrangian vertex multiplier.

Previous native `Lagrangian/RunT3YukawaSpinorPhaseNative.wl`, `Lagrangian/ProbeT3BCMajoranaNativeChiral.wl` and `scripts/RunLODirectUVNativeFermionAction.wl` establish abstract package identities and capture actions; they expressly do **not** supply a numerical installed `GammaCC` sign. The new B1 native driver instead isolates original and H.c. terms of both Yukawas from the **assembled** `LInt`, and mass terms from `LFree`; these are exported with full `InputForm` and checks. It does not interpret unexamined `Times` prefactors as absolute physical coefficients. In particular, it does not create numerical CGs or assert a native numerical 4x4 C.

## 2. Explicit spinor and Grassmann convention

Use Minkowski signature `(+,−,−,−)`, `sigma_i` standard Pauli, `epsilon = [[0,1],[-1,0]]`, and

\[
\gamma^0=\begin{pmatrix}0&I_2\\I_2&0\end{pmatrix},\quad
\gamma^i=\begin{pmatrix}0&\sigma_i\\-\sigma_i&0\end{pmatrix},\quad
C= i\gamma^2\gamma^0 =\begin{pmatrix}\epsilon&0\\0&-\epsilon\end{pmatrix}.
\]

Here `C^T = −C`, `C^2 = −I_4`, `C gamma^{mu T} C^{-1} = − gamma^mu`. Define the left-handed scalar contraction `u v := u_0 v_1 − u_1 v_0`, with anticommuting components, and Hermitian conjugation reversing fermion order. For one charged Dirac four-spinor,

\[
F=\begin{pmatrix}\chi_0\\\chi_1\\-\bar\xi_1\\\bar\xi_0\end{pmatrix},\quad
\overline F=(-\xi_1,\xi_0,\bar\chi_0,\bar\chi_1),\quad
F^c=C\overline F^{\,T}=\begin{pmatrix}\xi_0\\\xi_1\\-\bar\chi_1\\\bar\chi_0\end{pmatrix}.
\]

For `L = (l0,l1,0,0)^T`, `bar L = (0,0,bar l0,bar l1)`, direct component multiplication gives

\[
\begin{aligned}
\overline L F^c &= (l\chi)^\dagger, \\
\overline L F   &= (l\xi)^\dagger,\\
\overline F F   &= \chi\xi+(\chi\xi)^\dagger.
\end{aligned}
\]

The exact Grassmann-polynomial evaluator computes these identities using explicit antisymmetric exterior products. With `xi=chi` (self-conjugate Majorana spinor) it also proves

\[
\overline F F=\chi\chi+(\chi\chi)^\dagger,\qquad
\chi\chi=2\chi_0\chi_1.
\]

Thus a canonical free four-spinor Majorana term `−MF barF F/2` corresponds to `−MF chi chi/2 + h.c.`; **no additional factor 1/2 is inserted into the Weyl symmetric mass matrix**. A factor 2 in the symmetric matrix would double the physical Majorana mass. The direction and complex conjugation of a *non-diagonal flavour MF* are conditional matrix-algebra tests, **not native Matchete verification**; the source uses a diagonal-flavour mass parametrisation.

## 3. Yukawa coefficient and field placement

Because the Matchete topology uses `+ y_i barL (...) scalar + h.c.`, its two original vertices become

\[
\mathcal L_Y = y_{1,ar} (l_a\chi_r)^\dagger C_1 S_1
+y_{2,ar}(l_a\xi_r)^\dagger C_2 S_2^* +\mathrm{h.c.}
\]

where `C1,C2` include source gauge-index contractions. Taking the Hermitian conjugates and using the **numerical** convention

\[
\mathcal L_Y=-\nu_a\psi_I\left(Y_{aIs}\phi_s+Z_{aIs}\phi_s^*\right)+\mathrm{h.c.},
\]

we obtain, in this chosen C and scalar basis,

\[
Z_{1;a,\chi_r,s}=-y_{1,ar}^*\,(C_1)^*_{\nu,\chi,S_1},\qquad
Y_{2;a,\xi_r,s}=-y_{2,ar}^*\,(C_2)^*_{\nu,\xi,S_2}.
\]

The `−` signs arise from the convention chosen for the **numerical Lagrangian**, not extra Feynman vertices. With Dirac `psi=(xi,chi)` and complex `3×3` flavour mass `D_{rt}` in the Weyl quadratic action,

\[
M_{\rm Dirac} = \begin{pmatrix}0&D\otimes I_{d_F}\\D^T\otimes I_{d_F}&0\end{pmatrix}.
\]

`Y_2` occupies the **first** (`xi`) block, `Z_1` occupies the **second** (`chi`) block. `D^T` is mathematically obligatory to make the symmetric mass matrix; the native restricted model does not supply a general complex off-diagonal D. A common sign reversal of *both* source-adjusted Y and Z corresponds to `nu -> −nu` and leaves quadratic neutrino masses unchanged. Hence the old provisional builder's `+conjugate(y_i CG_i)` is **not by itself a physical discrepancy**.

For Majorana `xi=chi`, the two vertices use the same heavy Weyl multiplet and no extra multiplicity. The y2 source index is paired through the triplet metric when moving from its conjugate real irrep into the chosen weight basis; the adapter leaves a pending relative phase/gauge-tensor dictionary rather than silently imposing it as a native measurement.

## 4. Triplet mass metric and weight charges

Let the source adjoint Cartesian generators be `(T_a)_{bc}=-i epsilon_{abc}` and let the columns of unitary `U` be the `T3` eigenvectors ordered `m=+1,0,−1`. With a source Cartesian invariant `delta`, the weight-basis mass tensor is **the transpose congruence**

\[
G=U^T U,\qquad G^T=G,\qquad G^\dagger G=I,
\qquad (m_i+m_j)G_{ij}=0.
\]

For the explicit eigenvectors `u_+=(1,i,0)/sqrt2`, `u_0=(0,0,1)`, `u_-=(1,−i,0)/sqrt2`,

\[
G=\begin{pmatrix}0&0&1\\0&1&0\\1&0&0\end{pmatrix}.
\]

Rephasing only `u0` by `i` makes `G_00=−1` (the familiar Condon–Shortley-like sign); that change is a **basis convention** only if the vertices are transformed with the same phase. The mass term is `M=MF \otimes G`, not `MF \otimes I_3` in the ordered weight basis. A diagonal identity metric would pair identical weights instead of opposite electric charges and fails the charge-selection rule. The native C_-1 `delta` contract is source-verified symbolically; B2 must verify the actual imported `U`, its phases and contracted y2 indices across packages.

## 5. Allowed rephasings, physical relative phases

For Dirac `F -> e^{ia} F`, `Fc -> e^{ib} Fc`, `S1 -> e^{ic}S1`, `S2 -> e^{id}S2`, the source monomials carry charge vectors

| Source monomial | `(a,b,c,d)` charge |
|---|---|
| `y1 nu Fc S1` | `(0,1,1,0)` |
| `y2 nu F S2*` | `(1,0,0,−1)` |
| `MF F Fc` | `(1,1,0,0)` |
| `lambdaT3 HH S1 S2*` | `(0,0,1,−1)` |

For Majorana, set `Fc=F`; then the mass has charge `(2,0,0,0)`. In both cases the source-spurion combination

\[
\boxed{y_1 y_2 M_F^*\lambda_{T3}^*}
\]

has zero total field-rephasing charge. Replacing `lambdaT3*` with `lambdaT3` does not. This is a statement about rephasing covariance, **not** the fully contracted loop result before the SU(2)/scalar map is known.

The relative product (and the relative source phases after accounting for fields) can matter. The sign `C -> −C` **alone** changes the y1 action at fixed source couplings, but the same change accompanied by the appropriate heavy-field/mass rephasing is harmless. In an isolated pure `y1*y2` mass amplitude, a single overall sign in the full light Majorana matrix is also unobservable; the hostile test must not call that a physically distinct mass spectrum.

## 6. Loop combinatorics, independent of numerical core

For canonically normalised real mass eigenstates `R_t` and positive heavy Takagi masses `m_k`, with

\[
\mathcal L_Y=-g_{akt}\nu_a N_kR_t+\mathrm{h.c.},
\]

the chirality-flipping one-loop neutrino self-energy at zero momentum, following the single contraction and the `1/2` identical-external-Majorana quadratic-action convention, is

\[
\Sigma_{ab}(0)=\frac{1}{16\pi^2}\sum_{k,t}g_{akt}g_{bkt}m_k
B_0^{\overline{\rm MS}}(0;m_k^2,m_{S,t}^2;\mu^2).
\]

`B0_fin = 1 - [mF^2 ln(mF²/mu²) − mS² ln(mS²/mu²)]/(mF²−mS²)` for unequal positive masses and `−ln(m²/mu²)` for equal masses. There is **no extra Majorana 2** and no extra Majorana `1/2` multiplying this self-energy. For `phi=(R+iI)/sqrt(2)` and one complex Yukawa `h`, the couplings are `h/sqrt(2)` and `i h/sqrt(2)`; therefore

\[
\Sigma=\frac{h^2m_F}{32\pi^2}\left[B_0(m_R^2)-B_0(m_I^2)\right].
\]

This is the canonical real-scalar normalization identity; it is independently tested, and must **not** be confused with verifying the source's scalar `lambdaT3` convention. The latter involves the actual `T3MixCG` and its Higgs/charged scalar map, assigned to B2.

## 7. Hostile controls and pending checks

The isolated tests reject an unaccompanied flip of the C sign at the four-spinor action level; reject the wrong Majorana `1/2` through a doubled heavy mass; reject a wrongly conjugated complex y1 for a generic three-flavour mass matrix; distinguish an isolated Y/Z swap at fixed complex scalar mass mixing from an allowed full scalar conjugation; reject the incorrect identity triplet metric/charge pairing; reject an improper non-symmetric complex Dirac mass block; check invariance under consistent heavy-fermion rephasings; reject false native-certificate promotion and unsupported source drift. All these are **synthetic/exact mathematical controls**.

The three native reports, if generated using `wolframscript`, establish only that source `LFree`, `LInt`, original/H.c. Yukawa terms, projectors and symbolic mass terms were observed in a fresh kernel. They must be inspected before any claim about their precise couplings. The driver is intentionally fail closed if its term extractor cannot distinguish original from H.c. Its output does **not** claim numerical `GammaCC`, off-diagonal source flavour masses, nor live SU(2) basis alignment.

**Final B1 blocker**: Matchete's installed absolute C convention and RGBeta's concrete epsilon/CG index ordering relative to `U.T @ U` are not yet measured together in a single numerical component basis. B2 must supply that dictionary. No benchmark is physically certified by this batch, even if all sixteen conditional tensor residuals vanish.

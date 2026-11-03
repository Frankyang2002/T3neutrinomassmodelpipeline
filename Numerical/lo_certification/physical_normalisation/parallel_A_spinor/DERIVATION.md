# A1 — Independent four-component → Weyl convention

## Scope, orientation and limits

These equations define **one exact, declared** spinor representation, not an experimentally observed numerical `Matchete` `GammaCC` convention. The source's `+y1`, `+y2` and their Hermitian conjugates, and the `-MF`/`-MF/2` mass action coefficients, were verified by existing installed-Matchete witnesses; the current package reads them as *prior evidence*. The bridge of the installed symbolic `GammaCC` to the particular explicit matrix below remains unmeasured.

Follow the **actual Matchete original-term orientation**, not merely a symbolic spurion mnemonic: the original `y1` source multiplies `Bar[L] CConj[F] S1`, and the original `y2` multiplies `Bar[L] F Bar[S2]`. This is important when rephasing complex masses and taking Hermitian conjugates.

Use Minkowski metric `diag(+1,-1,-1,-1)`, `gamma^0=[[0,I],[I,0]]`, `gamma^i=[[0,sigma_i],[-sigma_i,0]]`, `gamma5=diag(-I,I)`, `epsilon_spinor=[[0,1],[-1,0]]`, and `C=+i gamma2 gamma0=diag(epsilon_spinor,-epsilon_spinor)`. Then `C^T=-C`, `C^2=-I`, and `C gamma^{mu T} C^{-1}=-gamma^mu`. With anticommuting two-component spinors, `u v = u0 v1-u1 v0` and `(u v)^dagger` reverses the ordering of every factor.

For a Dirac multiplet define

```
F        = (chi0,chi1,-bar(xi1),bar(xi0))^T
F^c      = (xi0,xi1,-bar(chi1),bar(chi0))^T
bar(F)   = (-xi1,xi0,bar(chi0),bar(chi1))
L        = (nu0,nu1,0,0)^T.
```

Exact four-component multiplication in exterior Grassmann algebra gives

```
bar(L) F^c = (nu chi)^dagger
bar(L) F   = (nu xi)^dagger
bar(F) F   = chi xi + (chi xi)^dagger.
```

For Majorana `xi=chi`: `chi chi = 2 chi0 chi1` and `bar(F) F=chi chi+(chi chi)^dagger`. The canonical source terms consequently read

```
Dirac:    L_M=-M chi xi + H.c.         [source 4-spinor prefactor -1]
Majorana: L_M=-M chi chi/2 + H.c.    [source 4-spinor prefactor -1/2]
```

The Majorana `1/2` avoids double-counting independent Grassmann fields in the **action**. **Do not apply another 1/2 to a Majorana mass-matrix entry or propagator, nor add a factor two for its self-conjugate field.** An off-diagonal complex family mass matrix is a **Weyl continuation**, not a general mass matrix native to the currently restricted Matchete source.

## Original UV couplings, H.c., chirality

For native component contractions `C1,C2` (including SU(2) indices) the original Matchete source terms are schematically

```
+ y1 (nu chi)^dagger C1 S1 + H.c.
+ y2 (nu xi)^dagger C2 S2* + H.c.
```

The neutrino-holomorphic **H.c.** vertices are therefore

```
+ conj(y1*C1) (nu chi) S1* 
+ conj(y2*C2) (nu xi) S2.
```

When the independent numerical convention is `L_Y=-g1 (nu chi) S1* -g2 (nu xi) S2 + H.c.`, the vertices are **`g1=-conj(y1*C1)`** and **`g2=-conj(y2*C2)`**. The minus signs convert a *positive source Lagrangian coupling* to the chosen *negative numerical Lagrangian convention*; they are not two additional diagrammatic `-1` factors. `y1` contains source `CConj[F]` in A_+0, whereas B_-1/C_-1 have a self-conjugate F; `y2` has non-charge-conjugated F in all three. Original Matchete projector `Proj[1]` (its right-projector index) and H.c. `Proj[-1]` are read from the native action; their translation is the explicit component statement above, not an inference from the sign of the projector number alone.

For Dirac source mass `D` with Weyl ordering `(xi,chi)`:

```
M_W = [[0, D ⊗ I_dF], [D^T ⊗ I_dF, 0]],
```

provided D is defined as the **xi–chi** holomorphic mass coefficient; `y2` acts on the `xi` block and `y1` on `chi`. If reordering `F`/`F^c` or using opposite Dirac mass index ordering, transpose D consistently. With Majorana triplet in the charge-weight basis, `M_W = M_flavour ⊗ G`, with `M_flavour=M_flavour^T`.

## SU(2) tensors and triplet Majorana metric

The fundamental `SU(2)` invariant is `epsilon_12=+1`, `epsilon_21=-1`. It is **not** the identity `delta_ij` for two fundamentals. `SU2Invariants.wl` builds `T3Y1CG` using conjugation flags `(True,True,False)` and `T3Y2CG` using `(True,False,True)`. The RGBeta Dirac doublet mass uses `eps`, whereas its neutral singlet and triplet Majorana masses use `1` and `del`, respectively. These are symbolic construction facts, **not** yet a numerical Matchete↔RGBeta phase alignment.

For the source Cartesian triplet generator

```
T3_cart = [[0,0,-i], [0,0,0], [i,0,0]]
```

choose normalized `T3` eigenvectors `u_+=(-i,0,1)/sqrt(2)`, `u_0=(0,1,0)` and `u_-=(i,0,1)/sqrt(2)`; `U=(u_+,u_0,u_-)` converts weight components to Cartesian. In the weight basis the Cartesian Majorana metric `delta` becomes **transpose congruence**

```
G=U^T U = [[0,0,1], [0,1,0], [1,0,0]].
```

This is symmetric and unitary and obeys `T3_weight^T G + G T3_weight=0`. The weight-basis identity matrix fails this charge-selection identity. A rephase `u_0 -> i u_0` instead yields `G_00=-1` with all matching vertices also rephased; `+/-` for the neutral triplet metric is therefore a **complete-basis choice**, not an independently measurable observable. C_-1's neutrino-neutral component is the `m=0` triplet state, contributing one heavy neutral species per flavour, not threefold multiplicity.

The existing installed-native charge-selected neutral CG values, transcribed from previous native scalar/CG witnesses, are:

| Model | C1 | C2 | Cm |
|---|---|---|---|
| A_+0 | 1 | `-i sqrt(2) 3^(-1/4)` | `+i sqrt(2) 3^(-1/4)` |
| B_-1 | 1 | 1 | `-2/sqrt(3)` |
| C_-1 | `3^(-1/4)` | `-3^(-1/4)` | `-2/sqrt(3)` |

The C_-1 triplet metric is to be inserted before tracing source CG indices if a charged/off-diagonal channel is needed. The neutral chosen-weight channel has `G_00=+1`. The native model does **not** fix an independent `i` on each eigenvector unless its charge basis is recorded in the same test.

## Scalar mixing and independent loop normalisation

Set the **restricted**, zero-portal test point `H_2=v/sqrt(2)`, `phi=(S1,S2^0)`, `V2=phi^dagger A phi`, with

```
A = [[a, conj(A21)], [A21,b]],
A21 = -lambdaT3 v^2 Cm/2,  a,b>0.
```

The quartic source coefficient is **plus** `lambdaT3 H H S1 Bar[S2] + H.c.`; the physical mass-squared potential is minus the interaction Lagrangian. Thus `A21` has this sign. The non-trivial `Cm` phase and `1/2` from two Higgs VEVs are essential. No benchmark identity-portal mapping is assumed by this test.

For eigenvalues `s_-<s_+` of A define the finite divided difference

```
D(M^2;s_-,s_+) = - integral_0^1 dx log(1+(1-x)(s_+-s_-)/(x M^2+(1-x)s_-))/(s_+-s_-)
                 = [B0_fin(M^2,s_+)-B0_fin(M^2,s_-)]/(s_+-s_-),
B0_fin(M^2,s;mu^2=M^2) = 1 + s log(s/M^2)/(M^2-s).
```

For **real positive heavy mass eigenvalues**, the chosen-Weyl loop matrix is

```
(m_nu)ab = sum_r [M_r A21 D_r/(16 pi^2)] (g1_ar*g2_br + g2_ar*g1_br).
```

For **a complex holomorphic Weyl mass parameter** `M_r`, the chirality-flipping two-component propagator `<psi psi>` has numerator `M_r^*` and denominator `p^2-|M_r|^2`; the consistent complex generalisation is therefore **`M_r^* D(|M_r|²)`**, not `M_r D(|M_r|²)`. This is an independently derived Weyl statement; the current restricted installed `MF` actions were tested with real masses, **so it is not a native validation of complex `MF` support**. For an ordinary real mass, the difference is invisible. The negative complex-mass control explicitly separates `M` from `M*`.

The external Majorana quadratic action is `-1/2 m_ab nu_a nu_b + H.c.`: sum both distinct cross orientations `(1,2)` and `(2,1)` exactly once. No further `1/2` or factor `2` multiplies the self-energy. The reference one complex scalar `phi=(R+iI)/sqrt(2)` illustrates the same factor: vertices `h/sqrt(2)` and `i h/sqrt(2)` yield `h² M [B0(mR²)-B0(mI²)]/(32 pi²)` for real M (this identity alone is **not** the T3 scalar mixing check).

## Rephasings and exactly what is physical

Keep external `nu`, H fixed. Define *old fields in terms of new fields* `chi_old=e^{ia}chi_new`, `xi_old=e^{ib}xi_new`, `S1_old=e^{ic}S1_new`, `S2_old=e^{id}S2_new`. In the **original Matchete four-component coupling orientation** described above,

```
y1'      = exp[i(-a+c)] y1,
y2'      = exp[i(-b-d)] y2,
M'       = exp[i(a+b)] M,
lambda'  = exp[i(c-d)] lambda.
```

For Majorana `xi=chi`, enforce `b=a`, hence `M'=exp(2ia)M`. With `g1=-conj(y1*C1)` and `g2=-conj(y2*C2)`, the invariant physical one-loop **holomorphic-neutrino** factor is

```
conj(y1) * conj(y2) * conj(M) * lambdaT3
```

(up to the consistently rotated CGs/scalar eigenvectors and the loop mass function). Its overall phase is *not* separately observable: `m_nu -> -m_nu` preserves Takagi masses and the CP-odd invariant built from `m_nu m_nu^dagger`. An internal source-coupling sign reversal that violates the above phase transport is not a valid single-field rephasing **at fixed other couplings**, but a pure isolated one-contribution overall mass sign is still not measurable. The previous source-numerator mnemonic `y1*y2*M^**lambda^*` in `B1_PHASE_TRANSPORT_DERIVATION.md` explicitly defines monomials `nu chi S1`, i.e. a different holomorphic assignment from the *original four-component* `+y_i Bar[L] ... + H.c.` terms. Do not combine the two name orientations without transporting the conjugations. In particular this track's complex-`M` extension is **not** a claim that the old symbolic mnemonic has been natively invalidated.

For arbitrary flavour mass matrices, perform actual flavour-basis congruences and use the inverse chirality-flipping propagator; the simple elementwise diagonal mass formula is restricted to nondegenerate diagonal-flavour benchmark points. Neither a general complex non-diagonal Matchete free action nor a fully numeric SU(2) Matchete↔RGBeta coupling dictionary is certified here.

## Physical source-certification obstruction

Installed Matchete's `GammaCC` is an abstract Clifford symbol. An independent extraction of a convention-sensitive **relative** Feynman rule/operator coefficient after fixing the same external Weyl basis, C orientation, and SU(2) CG axes is still needed. `native_rule_probe.wl` exercises the installed `FeynmanRules` API where available and gives a wrong-H.c.-sign control, but **cannot itself fix an unmeasured numerical `GammaCC`↔Weyl change of basis**. Until that bridge is observed, the `eta1/eta2` relative factors must remain conditional; `0/16` physical source certificates.

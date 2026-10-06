# T3 Implementation-Specific Numerical Choices

This document records implementation choices in the T3 neutrino-mass pipeline that can change numerical results independently of the underlying model physics.

The purpose is to distinguish:

1. **physical model assumptions**;
2. **matching/RGE truncation choices**;
3. **numerical conventions**;
4. **fit/search choices**;
5. **implementation artefacts or failure modes**.

The current implementation was checked against the repository structure and source on the current `main` branch. These notes are intended to explain why two calculations of the same nominal T3 model can produce somewhat different intermediate expressions, fit quality, benchmark points, or low-energy observables.

---

# 1. Lagrangian and model implementation specifics

## 1.1 Hypercharge convention

The project uses

$$
Q=T_3+Y.
$$

Therefore

$$
Y(S_1)=\frac{\alpha}{2},
\qquad
Y(S_2)=\frac{\alpha+2}{2},
\qquad
Y(F)=\frac{\alpha+1}{2}.
$$

This is a convention rather than a physical restriction, but comparisons with literature must account for the factor of two. A literature value of hypercharge cannot be inserted directly unless the convention is translated consistently.

---

## 1.2 Restricted production representation set

The standard production model classes are

$$
\begin{aligned}
A &: (d_{S_1},d_{S_2},d_F)=(1,3,2),\\
B &: (2,2,1),\\
C &: (2,2,3),\\
D &: (3,1,2),\\
E &: (3,3,2).
\end{aligned}
$$

Normal production support is restricted to singlet, doublet, and triplet representations.

The topology validator separately requires

$$
d_{S_i}=d_F\pm1
$$

for both scalar Yukawa vertices and requires the scalar product $S_1\otimes S_2$ to contain the $J=1$ channel needed for the Higgs-scalar quartic contraction.

Therefore the production model set is narrower than the mathematically possible set of larger-representation T3 constructions.

---

## 1.3 Neutral-state selection

The model-selection layer can test whether a BSM multiplet contains an electrically neutral component using

$$
Q=T_3+Y.
$$

This means the set of model/$\alpha$ points considered in scans can be filtered before the numerical calculation begins.

This is a model-selection implementation choice rather than a statement about the formal existence of the T3 topology.

---

## 1.4 Shared-scalar/scotogenic branch

Ordinary T3 models use physical fields

```text
F, S1, S2
```

while the shared-scalar branch uses

```text
F, S
```

with one physical scalar `S` filling the formal topology roles `S1` and `S2` up to conjugation.

The formal expansion into `S1,S2` occurs only at the matching boundary.

The currently supported shared-scalar numerical branch requires

$$
d_S=2,
\qquad
d_F=1\ \text{or}\ 3,
\qquad
\alpha=-1.
$$

Therefore the shared-scalar branch should not be interpreted as merely taking the ordinary theory and numerically setting $S_1=S_2$. The physical field counting and matching bookkeeping are different.

---

## 1.5 Heavy-flavour multiplicity

The numerical implementation uses three heavy fermion generations. The T3 Yukawas are therefore $3\times3$ lepton-heavy flavour matrices.

This fixes the available rank and flavour freedom of the generated neutrino-mass matrix.

A one-heavy-generation implementation of the same model would generally have less flavour freedom and can give qualitatively different fit behaviour.

---

## 1.6 Majorana versus Dirac heavy fermions

A heavy fermion is treated as self-conjugate when

$$
Y_F=0
$$

and its SU(2) representation is real, implemented for odd-dimensional SU(2) irreps.

For $$Y(F^c)=-Y(F)$$

This changes:

- the heavy-fermion degrees of freedom;
- mass-matrix symmetry requirements;
- allowed contractions;
- beta-function structure;
- numerical state validation.

This is particularly relevant for $B,\alpha=-1$, where the singlet heavy fermion is Majorana.

---

## 1.7 Yukawa orientation and conjugation

The two T3 Yukawa interactions are implemented schematically as

$$
\bar L\,F^c\,S_1
$$

and

$$
\bar L\,F\,S_2^\dagger.
$$

The implementation therefore fixes where charge conjugation and complex conjugation appear in the final Weinberg coefficient.

Another implementation can use an equivalent field convention with different explicit conjugations. Raw expressions can then look different even when the physical result is equivalent after translating conventions.

---

## 1.8 T3 scalar-mixing quartic

The topology quartic is implemented schematically as

$$
\lambda_{T3}\,H H S_1 S_2^\dagger+\mathrm{h.c.}
$$

with representation-dependent Clebsch-Gordan tensors.

This convention fixes which version of $\lambda_{T3}$, or its conjugate, appears in matched coefficients.

---

## 1.9 Symmetric Higgs-Higgs channel

The two Higgs fields in the $HH S_1S_2^\dagger$ vertex are identical bosons.

The implementation explicitly projects the two-Higgs SU(2) tensor product onto the symmetric subspace before defining the Clebsch-Gordan tensor.

This removes antisymmetric Higgs-Higgs contractions by construction.

---

## 1.10 Pseudoreal SU(2) index conventions

Even-dimensional SU(2) representations are pseudoreal. The Wolfram/Matchete implementation therefore distinguishes the orientation of conjugated and unconjugated indices using `CRep` and barred representation objects.

Real is when our conjugate representation $R^*=R$, while pseudoreal requires an intertwiner for even representation

This is important for doublets and other pseudoreal irreps.

Different tensor conventions can change:

- Clebsch-Gordan signs;
- basis phases;
- intermediate symbolic expressions.

These differences are not necessarily physical, provided the entire calculation uses one convention consistently.

---

## 1.11 Topology vertices require a unique invariant

For topology-defining interactions, the implementation requires Matchete to return exactly one invariant tensor.

If more than one independent invariant exists, the topology vertex currently fails instead of introducing several independent couplings.

This is a significant implementation restriction for larger representations or forced model points.

It means the current code does not represent the most general theory whenever a topology vertex admits multiple independent SU(2) contractions.

We assume that we only have one coupling, where if we have multiple combinations then we would need multiple couplings to represent them which is not what we want. So if we have more than one invariant tensor then we say it is invalid

---

## 1.12 First-valid alternative selection

Interaction candidates can belong to an `AlternativeGroup`.

After Matchete validation, the implementation keeps the first valid candidate in each alternative group.

Therefore, if two mathematically valid alternatives are generated in one group, the selected representative depends on candidate ordering.

For currently supported models this is intended to choose one equivalent contraction representation, but it is an implementation-specific basis choice.

---

## 1.13 Scalar-potential invariant basis

The scalar potential is built by:

1. generating all SU(2) invariant tensors;
2. projecting identical-boson pairs onto symmetric tensor subspaces;
3. flattening tensors;
4. using matrix-rank tests to retain a linearly independent basis.

Therefore the scalar coupling basis is algorithmically chosen.

Equivalent bases can produce different-looking RGEs and different coordinates in numerical parameter space even if they describe the same physics.

---

## 1.14 Real versus complex scalar couplings

Many scalar self-couplings and portal couplings are defined as self-conjugate and therefore real. We decided it ourselves

We know they are self conjugate as we look at the operator, like $H^\dagger H$ which is self conjugate. Many of these are self conjugate and the couplings with that can be treated as real as they have their hermitian conjugate in the lagrangian itself, so the imaginary component does not contribute

The topology coupling $\lambda_{T3}$ is complex as we have HHS1S2^dagger, which is not self-conjugate so it is not real

Some representation-specific quartics are also complex.

This restricts the scalar phase space compared with the most general complex parametrisation.

---

## 1.15 Special $B,\alpha=-1$ scalar sector

For the doublet-doublet case with

$$
(d_{S_1},d_{S_2},d_F)=(2,2,1),
\qquad
\alpha=-1,
$$

the numerical state contains extra quartic structures beyond the generic five-coupling scalar sector.

Several are initialized to zero in the default benchmark but remain part of the RGE state.

Therefore this model has a larger running coupling space than a simplified treatment that retains only

$$
\lambda_{S_1},\lambda_{S_2},\lambda_{H1},\lambda_{H2},\lambda_{12},\lambda_{T3}.
$$

| Code name                    | Type in numerical state | Schematic field content / meaning                                                  | Why it exists                                                                                                             |
| ---------------------------- | ----------------------- | ---------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `lambdaH1Adj`                | Real                    | Alternative/adjoint contraction of $H^\dagger H\,S_1^\dagger S_1$                  | Two SU(2) doublet bilinears can be contracted through a nontrivial triplet/adjoint channel as well as the singlet channel |
| `lambdaH2Adj`                | Real                    | Alternative/adjoint contraction of $H^\dagger H\,S_2^\dagger S_2$                  | Same reason for $S_2$                                                                                                     |
| `lambda12Adj`                | Real                    | Alternative/adjoint contraction of $S_1^\dagger S_1\,S_2^\dagger S_2$              | Additional independent SU(2) contraction                                                                                  |
| `lambdaHHdagS2S2`            | Complex                 | $H^\dagger H^\dagger S_2S_2+\mathrm{h.c.}$                                         | Allowed because $2Y(H^\dagger)+2Y(S_2)=0$                                                                                 |
| `lambdaHHdagS1barS1bar`      | Complex                 | $H^\dagger H^\dagger S_1^\dagger S_1^\dagger+\mathrm{h.c.}$                        | Hypercharge also cancels specifically for $Y(S_1)=-1/2$                                                                   |
| `lambdaS1bar2S2bar2`         | Complex                 | $(S_1^\dagger)^2(S_2^\dagger)^2+\mathrm{h.c.}$                                     | Extra hypercharge-neutral doublet quartic                                                                                 |
| `lambdaS1barS2S2bar2`        | Complex                 | Schematic $S_1^\dagger S_2(S_2^\dagger)^2+\mathrm{h.c.}$                           | Special hypercharges permit another non-self-conjugate quartic                                                            |
| `lambdaS1S1bar2S2bar`        | Complex                 | Schematic $S_1(S_1^\dagger)^2S_2^\dagger+\mathrm{h.c.}$                            | Conjugate-pattern-specific quartic allowed in this representation                                                         |
| `lambdaHHdagS1barS2barCross` | Complex                 | Additional independent contraction involving $H,H^\dagger,S_1^\dagger,S_2^\dagger$ | Same field content can admit more than one independent SU(2) contraction                                                  |

---

## 1.16 Matching truncation

Production matching targets

$$
L=1,
\qquad
d_{\rm EFT}\le5.
$$

Therefore the implementation is a one-loop, dimension-five-truncated EFT calculation.

Effects that first enter through dimension-six operators are intentionally outside production scope.

---

## 1.17 Scalar-first threshold orderings are excluded

If a scalar is integrated out before the fermion, the intermediate EFT can contain a leading dimension-six structure schematically of the form

$$
\frac{y\lambda_{T3}}{M_S^2}LFHHS.
$$

Because production is truncated at dimension five, a scalar-first calculation would omit the leading intermediate EFT operator.

The pipeline therefore rejects scalar-first threshold orderings.

This is an implementation/truncation restriction, not a statement that scalar-first T3 physics is inconsistent.

---

## 1.18 Partially split scalar thresholds are not implemented

The verified hierarchical production path is

$$
F\rightarrow(S_1,S_2)
$$

or, in the shared-scalar case,

$$
F\rightarrow S.
$$

The two ordinary physical scalars are removed together.

Therefore a hierarchy such as

$$
M_F\gg M_{S_1}\gg M_{S_2}
$$

is not treated with three independent matching thresholds.

---

## 1.19 Weinberg orientation and Majorana symmetrisation

The matching layer extracts one holomorphic Weinberg orientation and avoids separately counting its Hermitian conjugate.

The physical Majorana coefficient is then formed by symmetrising the external flavour indices,

$$
C_{5,pq}^{\rm phys}
=
C_{5,pq}+C_{5,qp}.
$$

In one generation this introduces an ordered-to-physical factor of two.

This is an important convention when comparing with papers or packages that define the Weinberg operator with a different symmetrisation or $1/2$ normalization.

---

## 1.20 Flavor lifting of the one-generation coefficient

For common-threshold flavor matching, the one-generation matched coefficient is factorised as

$$
\kappa_{1g}
=
F_{\rm loop}\,
\overline y_1\,
\overline y_2.
$$

The full three-generation coefficient is then built as

$$
C_{5,pq}
=
\sum_r
F(M_r)
\left[
\overline y_{1,pr}\,\overline y_{2,qr}
+
\overline y_{2,pr}\,\overline y_{1,qr}
\right].
$$

This is effectively a heavy-mass eigenbasis treatment.

A calculation that instead evolves a fully non-diagonal heavy mass matrix without diagonalizing into individual heavy eigenstates is not numerically identical at intermediate steps.

---

# 2. RGE and EFT implementation specifics

## 2.1 One-loop running

The UV, intermediate, and final-EFT RGEs are implemented at one-loop order.

Therefore any comparison against two-loop running will differ by genuine higher-order corrections.

---

## 2.2 RGE variable and loop convention

The independent variable is

$$
t=\ln\mu
$$

with beta functions interpreted as

$$
16\pi^2\frac{dX}{d\ln\mu}
=
\beta_X^{(1)}.
$$

A source using $\ln\mu^2$ instead of $\ln\mu$ differs by a factor of two in the derivative convention.

---

## 2.3 Hypercharge gauge-coupling normalization

The final SM+Weinberg numerical stage uses the ordinary SM hypercharge coupling $g_Y$, not GUT-normalized

$$
g_1=\sqrt{\frac53}\,g_Y.
$$

Beta coefficients therefore must be translated before comparison with formulas written in the GUT-normalized convention.

---

## 2.4 Higgs quartic normalization

The final numerical stage uses

$$
V(H)
=
\frac{\lambda_H}{2}(H^\dagger H)^2.
$$

This differs from the common convention

$$
V(H)=\lambda(H^\dagger H)^2.
$$

The two are related by

$$
\lambda_H=2\lambda.
$$

Therefore both the Higgs beta function and the Weinberg-operator beta function must be compared using the same quartic normalization.

---

## 2.5 Final Weinberg RGE

The implemented final-EFT equation is

$$
16\pi^2\frac{dC_5}{d\ln\mu}
=
\left(
2\lambda_H
-
3g_2^2
+
2T
\right)C_5
-
\frac32
\left[
Y_eY_e^\dagger C_5
+
C_5(Y_eY_e^\dagger)^T
\right],
$$

where

$$
T=
{\rm Tr}
\left(
Y_eY_e^\dagger
+
3Y_uY_u^\dagger
+
3Y_dY_d^\dagger
\right).
$$

The coefficient of $\lambda_H$ is tied to the project quartic convention above.

---

## 2.6 Full-flavour SM Yukawa evolution

The numerical final-EFT stage evolves full complex $3\times3$ matrices

$$
Y_e,\qquad Y_u,\qquad Y_d.
$$

The implementation is therefore not a diagonal-Yukawa approximation.

Basis rotation effects can occur during the running.

---

## 2.7 UV and intermediate beta functions come from RGBeta

The UV and scalar-only intermediate renormalisable beta functions are imported from RGBeta-generated payloads.

The numerical result therefore inherits:

- RGBeta field conventions;
- coupling definitions;
- representation basis;
- exported beta-function normalization;
- the chosen generated coupling closure.

A hand-written reduced RGE system can differ even if it retains the same nominal interactions.

---

## 2.8 Scalar-only intermediate EFT

After the fermion threshold, the implemented hierarchical theory is

$$
{\rm SM}+S_1+S_2
$$

or

$$
{\rm SM}+S
$$

for the shared branch.

This is why only the verified fermion-first hierarchy is available in production.

---

## 2.9 LLSS Wilson coefficient

At the fermion threshold the code produces dimension-five scalar-lepton Wilson operators of LLSS type.

These are transported through the scalar-only intermediate interval and are the source of the hierarchical logarithmic contribution to the final Weinberg coefficient.

---

## 2.10 Strict fixed-order treatment of LLSS self-running

The authoritative one-loop final Weinberg coefficient deliberately does not include the one-loop self-running of the LLSS Wilson coefficient followed by its scalar-loop conversion into $C_5$.

The order counting is

$$
C_{\rm LLSS}^{(0)}
\xrightarrow{\text{1-loop self-running}}
C_{\rm LLSS}^{(1)}
\xrightarrow{\text{scalar loop}}
C_5^{(2)}.
$$

Therefore that effect is treated as two-loop order and excluded from the one-loop authoritative result.

This differs from an RG-improved calculation that resums or partially includes such terms.

---

## 2.11 Direct LLSS to Weinberg mixing

The code does include the direct one-loop generation of the Weinberg operator across the scalar-only interval.

This contribution carries the logarithmic dependence on the ratio of fermion and scalar threshold scales.

Another implementation can instead obtain the same logarithm by expanding the full one-loop matching result rather than treating it as explicit EFT running.

---

## 2.12 Hard plus direct decomposition

For the hierarchical path the authoritative coefficient is organized schematically as

$$
C_5(M_S)
=
C_{5,\rm hard}(M_F)
+
\Delta C_{5,\rm direct}(M_F\to M_S).
$$

The split between the two pieces depends on matching scheme and matching scale.

Only the combined fixed-order result should be compared physically.

---

## 2.13 Matching scale set to $M_F$

The renormalized hard matching expression is evaluated at the fermion threshold by setting recognized logarithms of the form

$$
\ln\frac{\bar\mu^2}{M_F^2}
$$

to zero.

Using a different matching scale changes the partition between hard matching and running by higher-order terms.

---

## 2.14 MSbar pole subtraction

The hierarchical final-$C_5$ builder removes the Matchete UV pole from the hard matching expression and uses an MSbar interpretation.

The current implementation recognizes specific textual forms of the Matchete $1/\epsilon$ contribution.

Therefore this step depends partly on the exact symbolic serialization produced by Matchete.

---

## 2.15 Pole/RGE consistency validation

The direct-running logarithm must satisfy the one-generation relation

$$
\text{direct logarithm coefficient}
=
2\times
\text{hard UV-pole residue}.
$$

The final hierarchical construction is blocked if this validation is not satisfied.

This is a strong consistency check, but it also means a symbolic parsing or convention mismatch can stop a point even if the intended physics is otherwise correct.

---

## 2.16 Integrated-out mass parameters are kept real

After a heavy field is removed, its mass parameter can remain inside Wilson coefficients.

The implementation re-registers integrated-out masses such as

$$
M_F,\qquad M_{S_1},\qquad M_{S_2}
$$

as self-conjugate couplings so that later symbolic operations continue to treat

$$
\overline M=M.
$$

This is symbolic bookkeeping rather than new physics, but it affects algebraic simplification and Hermitian conjugation.

---

## 2.17 Majorana $M_F$ is projected onto the symmetric manifold

For Majorana heavy fermions, each UV numerical RHS evaluation reconstructs the physical state using

$$
M_F
\rightarrow
\frac12(M_F+M_F^T).
$$

The saved numerical solution is also projected onto this manifold.

This suppresses antisymmetric numerical drift.

If the generated beta function contains a small antisymmetric numerical component, the implemented trajectory follows the projected symmetric system rather than the raw ODE.

---

## 2.18 Final $C_5$ is symmetrized numerically

After final SM+Weinberg running,

$$
C_5
\rightarrow
\frac12(C_5+C_5^T).
$$

This removes small numerical violations of the Majorana symmetry condition.

Large asymmetries should still be regarded as a bug or consistency failure rather than something that should be hidden by the projection.

---

## 2.19 Fixed electroweak vev

The numerical neutrino mass matrix is constructed using

$$
m_\nu
=
-\frac{v^2}{2}C_5
$$

with

$$
v=246.22~{\rm GeV}.
$$

The vev itself is not evolved as a running quantity.

A scheme using a running electroweak vev can therefore produce slightly different absolute neutrino masses and mass splittings.

---

## 2.20 Charged-lepton basis rotation

At each saved scale, the charged-lepton Yukawa matrix is diagonalized numerically using an SVD.

The singular values are sorted in ascending order and the corresponding left singular vectors are used to rotate the Majorana mass matrix,

$$
m_\nu
\rightarrow
U_e^T m_\nu U_e.
$$

This fixes the charged-lepton flavour basis numerically.

For realistic charged-lepton Yukawas this is stable, but near accidental singular-value degeneracies the basis can become numerically ambiguous.

---

# 3. Threshold-scale implementation choices

## 3.1 Default fermion masses versus fermion threshold

The default benchmark contains

$$
M_F=
{\rm diag}
(100,110,120)\ {\rm TeV}
$$

but uses one fermion threshold

$$
\mu_F=100\ {\rm TeV}.
$$

Therefore all three heavy generations are transferred into the lower EFT at a common scale even though their masses are not identical.

This is a matching prescription rather than exact sequential decoupling of the three mass eigenstates.

Residual finite logarithms associated with

$$
\ln\frac{110}{100},
\qquad
\ln\frac{120}{100}
$$

are therefore treated differently from a generation-by-generation threshold calculation.

---

## 3.2 Default scalar masses versus scalar threshold

The default ordinary benchmark contains

$$
m_{S_2}=1.000~{\rm TeV}
$$

and

$$
m_{S_1}
=
\sqrt{1.306122448979592\times10^6}\ {\rm GeV}
\approx1.143~{\rm TeV}.
$$

Both physical scalars are removed at

$$
\mu_S=1.000~{\rm TeV}.
$$

Therefore scalar decoupling is also performed at one common scalar threshold rather than one threshold per physical scalar mass.

This can produce percent-scale differences relative to a fully sequential threshold implementation, depending on couplings and hierarchy.

---

# 4. Numerical implementation specifics

## 4.1 ODE solver

The numerical RGE system uses SciPy `solve_ivp` with the high-order explicit solver

```text
DOP853
```

for the UV, intermediate, and final SM+Weinberg stages.

A different integrator such as RK45, Radau, or BDF can produce slightly different trajectories at finite tolerances.

---

## 4.2 ODE tolerances

The default tolerances are

$$
{\tt rtol}=10^{-8},
\qquad
{\tt atol}=10^{-11}.
$$

These set the numerical integration accuracy.

Because the state contains couplings and matrix elements with different magnitudes, componentwise error control can affect very small entries differently from order-one couplings.

---

## 4.3 Maximum step size

The UV and intermediate runners default to

$$
{\tt max\_step\_log}=\infty.
$$

Therefore the adaptive solver is free to take large steps in $t=\ln\mu$ when the local error estimate permits it.

This normally has negligible effect with the adopted tolerances, but it is still an implementation choice.

---

## 4.4 Complex quantities are integrated as real vectors

Every complex scalar or matrix is converted into separate real and imaginary degrees of freedom.

For example,

$$
A
\rightarrow
\left(
{\rm vec}\,{\rm Re}A,\,
{\rm vec}\,{\rm Im}A
\right).
$$

The ODE solver therefore applies its error control independently to real and imaginary components.

---

# 5. Neutrino-observable implementation specifics

## 5.1 Custom Takagi factorisation

The neutrino Majorana matrix is factorized using a custom realified symmetric eigenproblem.

For

$$
M=A+iB,
$$

the implementation constructs

$$
\mathcal M_{\rm real}
=
\begin{pmatrix}
A & -B\\
-B & -A
\end{pmatrix}
$$

and solves

$$
\mathcal M_{\rm real}
\binom{x}{y}
=
\sigma
\binom{x}{y}.
$$

The complex Takagi vector is then

$$
u=x+iy.
$$

This differs numerically from an SVD-based Takagi implementation, especially in nearly degenerate singular-value subspaces.

---

## 5.2 Input mass matrix symmetry

Before the factorisation, the mass matrix must satisfy approximately

$$
M_\nu=M_\nu^T
$$

using a numerical tolerance.

It is then explicitly projected as

$$
M_\nu
\rightarrow
\frac12(M_\nu+M_\nu^T).
$$

Small antisymmetric floating-point errors are therefore removed.

---

## 5.3 Takagi branch selection

The realified $6\times6$ eigenproblem produces eigenvalues in approximately $\pm\sigma_i$ pairs.

The implementation selects the three largest eigenvalues as the physical positive Takagi branches.

This is robust for ordinary non-degenerate spectra but can become delicate when the lightest singular value is extremely small or two singular values are nearly degenerate.

---

## 5.4 Takagi phase fixing

After constructing the Takagi vectors, each column phase is adjusted using

$$
U_i
\rightarrow
U_i
\exp\left[
-\frac{i}{2}
\arg
\left(
(U^TM_\nu U)_{ii}
\right)
\right].
$$

This makes the diagonal masses non-negative.

It affects phases but not the masses or $|U_{\rm PMNS}|$ away from numerical degeneracies.

---

## 5.5 Takagi ordering

The singular values are sorted from smallest to largest.

For normal ordering this is interpreted as

$$
[m_1,m_2,m_3].
$$

For inverted ordering, the ascending Takagi output is interpreted as

$$
[m_3,m_1,m_2]
$$

and then relabelled to

$$
[m_1,m_2,m_3].
$$

Near level crossings the eigenstate labels can change discontinuously.

---

## 5.6 Automatic ordering heuristic

When the ordering is `AUTO`, the code compares

$$
m_2^2-m_1^2
$$

and

$$
m_3^2-m_2^2.
$$

The spectrum is labelled normal if the lower adjacent gap is smaller, otherwise inverted.

This is an algorithmic ordering prescription rather than a fit over both hierarchy hypotheses.

---

## 5.7 Production benchmark ordering

The default benchmark configuration explicitly uses

```text
"ordering": "NO"
```

so the benchmark search is performed against the normal-ordering target.

An inverted-ordering solution is not allowed to compete in that scan.

---

## 5.8 Takagi residual cut

The reconstructed factorization must satisfy a relative residual below approximately

$$
10^{-7}.
$$

A point exceeding this threshold is treated as a numerical failure.

Therefore some failed model points can be failures of the chosen decomposition/residual criterion rather than physical exclusions.

---

## 5.9 Matrix symmetry tolerances

Different numerical objects use slightly different `allclose` tolerances.

Typical examples are of order

$$
{\tt rtol}\sim10^{-9}\text{--}10^{-10}
$$

with different absolute tolerances for $M_\nu$, $C_5$, and $M_F$.

Therefore a very small matrix can pass or fail symmetry validation depending on its absolute scale.

---

# 6. Default benchmark choices

The default numerical benchmark uses

$$
\mu_{\rm UV}=10^7~{\rm GeV},
$$

$$
\mu_F=10^5~{\rm GeV},
$$

$$
\mu_S=10^3~{\rm GeV},
$$

and

$$
\mu_{\rm low}=100~{\rm GeV}.
$$

The default SM couplings at the UV input scale are approximate benchmark values,

$$
g_Y=0.36,
\qquad
g_2=0.65,
\qquad
g_3=0.6,
\qquad
\lambda_H=0.25.
$$

The SM Yukawa matrices are initialized as diagonal approximate numerical matrices.

These are not precision SM boundary conditions obtained by running measured low-energy parameters upward.

Therefore another implementation using a precision SM matching procedure can produce different T3 running even if the T3 beta functions are identical.

---

# 7. Benchmark-search restrictions

## 7.1 Parameters varied in the default search

The default benchmark search varies:

$$
{\rm Re}\,\lambda_{T3}
$$

and the 18 real entries of

$$
y_1,\qquad y_2.
$$

The search range for the topology quartic is

$$
10^{-3}
\le
{\rm Re}\,\lambda_{T3}
\le
1
$$

with logarithmic sampling.

The Yukawa entries are sampled linearly over

$$
-0.5
\le
y_{1,ij},\,y_{2,ij}
\le
0.5.
$$

This gives a 19-dimensional search space.

---

## 7.2 Parameters held fixed

The default benchmark search does not vary:

- imaginary parts of $y_1$;
- imaginary parts of $y_2$;
- ${\rm Im}\lambda_{T3}$;
- heavy fermion masses;
- scalar masses;
- generic scalar quartics;
- representation-specific extra scalar quartics;
- SM input couplings.

Therefore the resulting benchmark is only the best point found in a restricted 19-dimensional slice of the full model parameter space.

It should not be interpreted as the global best fit of the complete T3 model.

---

# 8. Sobol benchmark-search specifics

The default automatic search uses a scrambled Sobol sequence with

$$
N_{\rm Sobol}=256,
$$

target

$$
\chi^2_{\rm target}=10,
$$

and retains the best

$$
N_{\rm keep}=10
$$

points.

The default randomization seed is

$$
20260927.
$$

Changing any of these changes the benchmark that is found without changing the physical model.

---

## 8.1 Early stopping

The Sobol stage stops when a point satisfies

$$
\chi^2\le10.
$$

Therefore different models may receive different numbers of full pipeline evaluations.

For example, one model can stop at

$$
\chi^2=9.9
$$

after relatively few evaluations while another continues long enough to find

$$
\chi^2=2.
$$

Therefore the stored best $\chi^2$ values cannot automatically be interpreted as a fair model-ranking statistic unless equal search budgets are enforced.

---

# 9. Local optimizer specifics

After the Sobol stage, the search can run local least-squares refinement.

The default refinement uses up to five starting seeds and progressively optimizes only the most locally sensitive

$$
5,\qquad8,\qquad12
$$

parameters.

The corresponding default evaluation budgets are approximately

$$
30,\qquad40,\qquad60
$$

function evaluations.

The full search has 19 varying parameters, so the default local stages do not optimize all 19 simultaneously.

This can materially affect the final reported benchmark.

---

## 9.1 Unit-coordinate parameterization

Every bounded parameter is mapped to

$$
z_i\in[0,1].
$$

Linear parameters use

$$
x
=
x_{\min}
+
z(x_{\max}-x_{\min}),
$$

while log-scaled parameters use

$$
\ln x
=
\ln x_{\min}
+
z
\left(
\ln x_{\max}
-
\ln x_{\min}
\right).
$$

Therefore local sensitivity rankings depend on the chosen parameter coordinates and bounds.

A parameter can appear more or less important simply because it is represented linearly versus logarithmically.

---

## 9.2 Finite-difference sensitivity step

The default local sensitivity step is

$$
\Delta z=0.02.
$$

This controls the estimate of which parameters are most influential.

A different step size can change the selected active-parameter subset and therefore the local optimum found.

---

# 10. Numerical-failure penalty in optimization

When a full model evaluation throws an exception during local optimization, the optimizer receives a fixed residual penalty.

The default residual assigned to every fitted observable is approximately

$$
r_i=10^4.
$$

For five fitted observables this corresponds to

$$
\chi^2_{\rm penalty}
\approx
5\times10^8.
$$

Therefore very different failure modes become the same flat high-$\chi^2$ region.

Examples include:

- Takagi failures;
- Majorana $M_F$ symmetry failures;
- ODE failures;
- invalid threshold boundaries;
- non-finite predictions.

The optimizer therefore cannot distinguish slightly numerically unstable regions from fundamentally invalid ones.

---

# 11. Oscillation-fit implementation specifics

The fitted observable vector is

$$
\left(
\Delta m_{21}^2,\,
\Delta m_{3\ell}^2,\,
\sin^2\theta_{12},\,
\sin^2\theta_{13},\,
\sin^2\theta_{23}
\right).
$$

For normal ordering,

$$
\Delta m_{3\ell}^2
=
\Delta m_{31}^2,
$$

while for inverted ordering,

$$
\Delta m_{3\ell}^2
=
\Delta m_{32}^2.
$$

The angles are extracted from $|U_{\rm PMNS}|$ using

$$
\sin^2\theta_{13}
=
|U_{e3}|^2,
$$

$$
\sin^2\theta_{12}
=
\frac{|U_{e2}|^2}{1-|U_{e3}|^2},
$$

$$
\sin^2\theta_{23}
=
\frac{|U_{\mu3}|^2}{1-|U_{e3}|^2}.
$$

---

## 11.1 Gaussianized NuFIT likelihood

The current NuFIT 6.1 target is not the official NuFIT likelihood surface.

The implementation:

1. takes the reported central values;
2. averages the quoted positive and negative one-sigma uncertainties;
3. treats each observable as Gaussian;
4. uses a diagonal covariance matrix;
5. ignores experimental correlations.

Therefore the pipeline $\chi^2$ is a benchmark-quality metric, not an exact reproduction of the NuFIT global $\Delta\chi^2$.

This can change the location and relative quality of numerical best fits.

---

## 11.2 Quantities not currently included in the fit

The benchmark fit does not include constraints on

$$
\delta_{\rm CP},
$$

the Majorana phases,

$$
m_{\rm lightest},
$$

$$
\sum_i m_i,
$$

or

$$
m_{\beta\beta}.
$$

A parameter point with an excellent five-observable oscillation $\chi^2$ is therefore not automatically a global phenomenological best fit.

---

# 12. Persistent benchmark and cache behaviour

The numerical orchestration layer stores generated model-specific numerical configs.

If a compatible generated config already exists, it is reused unless the run explicitly resets numerical configs.

After a successful benchmark search, the selected best parameters can be written back into the generated model config and benchmark searching disabled.

Therefore a later run can reproduce a previously frozen benchmark rather than perform a fresh search.

This means run history matters.

The commands

```powershell
python pipeline.py --full --numerical
```

and

```powershell
python pipeline.py --full --numerical --reset-numerical-configs
```

are not equivalent from the point of view of benchmark selection.

Two users on the same source commit can obtain different benchmark outputs if their local

```text
configs/generated_models/
```

or

```text
output/benchmarks/
```

contents differ.

---

# 13. Approximate importance of implementation effects

The following ranking is useful when interpreting differences between model points or between this pipeline and another implementation.

| Rank | Implementation effect                                                       | Expected importance                               |
| ---: | --------------------------------------------------------------------------- | ------------------------------------------------- |
|    1 | Restricted scan subspace: real Yukawas, fixed masses, fixed quartics/phases | Very high                                         |
|    2 | 256-point Sobol search and early stopping at $\chi^2=10$                    | Very high                                         |
|    3 | Reduced-dimensional local optimization using only 5/8/12 active variables   | Very high                                         |
|    4 | Cached/frozen model-specific benchmarks                                     | Very high                                         |
|    5 | Takagi branch, degeneracy and residual behaviour                            | High for failed/near-degenerate cases             |
|    6 | Strict Majorana $M_F$ symmetry validation and projection                    | High for affected Majorana models                 |
|    7 | Common fermion threshold for non-degenerate $M_{F_i}$                       | Moderate                                          |
|    8 | Common scalar threshold for non-degenerate scalar masses                    | Moderate                                          |
|    9 | Gaussianized and uncorrelated NuFIT likelihood                              | Moderate to high for best-fit location            |
|   10 | Strict fixed-order treatment of LLSS self-running                           | Small to moderate, but conceptually important     |
|   11 | Approximate SM UV inputs and fixed electroweak vev                          | Small to moderate                                 |
|   12 | DOP853 and ODE tolerance choices                                            | Usually small                                     |
|   13 | Explicit numerical matrix symmetrization                                    | Usually tiny unless masking another inconsistency |

---

# 14. Known failure modes that should not automatically be interpreted physically

## 14.1 Majorana $M_F$ symmetry validation

For Majorana models, the implementation requires

$$
M_F=M_F^T
$$

to tight numerical tolerance.

The UV runner also explicitly projects $M_F$ back onto the symmetric manifold.

Therefore a failure caused by a small antisymmetric component should be interpreted first as a numerical/RGE consistency problem.

It should not automatically be interpreted as evidence that the underlying T3 model has no viable solution.

---

## 14.2 Takagi numerical failures

A failed Takagi residual check means that the implemented decomposition did not reconstruct the supplied Majorana mass matrix to the required tolerance.

It does not by itself imply that the T3 model cannot generate a physical neutrino spectrum.

This distinction is especially important for:

- nearly singular mass matrices;
- nearly degenerate singular values;
- very hierarchical matrix entries;
- points close to eigenstate level crossings.

---

# 15. Compact description of the implemented calculation

A concise summary of the implementation is

$$
\boxed{
\begin{array}{l}
\text{one-loop, dimension-five-truncated EFT},\\
\text{Matchete-generated SU(2) invariant/CG basis},\\
\text{common threshold or fermion-first sequential matching only},\\
\text{RGBeta-generated UV and intermediate RGEs},\\
\text{strict fixed-order treatment of LLSS Wilson running},\\
\text{common fermion and scalar decoupling scales},\\
m_\nu=-\dfrac{v^2}{2}C_5,\quad v=246.22~{\rm GeV},\\
\text{full-flavour final SM+Weinberg running},\\
\text{custom Takagi factorisation and charged-lepton basis rotation},\\
\text{Gaussian five-observable NuFIT benchmark likelihood},\\
\text{restricted Sobol plus reduced local parameter optimization}.
\end{array}
}
$$

The physical model determines the field content and interactions, but the numerical result is additionally a function of

$$
R
=
R\!\left(
\begin{array}{c}
\text{operator basis},\\
\text{CG convention},\\
\text{matching scheme},\\
\text{threshold scales},\\
\text{loop/EFT truncation},\\
\text{RGE basis},\\
\text{ODE tolerances},\\
\text{Takagi prescription},\\
\text{fit likelihood},\\
\text{scan bounds},\\
\text{optimizer budget},\\
\text{cached benchmark}
\end{array}
\right).
$$

This separation should be kept explicit when comparing numerical results between T3 models or against other implementations.

---

# 16. Practical interpretation

When quoting a benchmark result from this pipeline, the most precise description is not simply

> “the T3 model predicts this point.”

A more accurate statement is

> “this is the benchmark found for the specified T3 representation within the implemented one-loop, dimension-five EFT treatment, threshold prescription, numerical conventions, restricted parameter scan, and Gaussianized oscillation-fit metric.”

For most ordinary successful points, ODE discretization and floating-point details should be much smaller than the dominant implementation dependence from:

1. threshold prescription;
2. scan subspace;
3. finite search budget;
4. benchmark caching;
5. Takagi handling near singular or degenerate spectra.

These should therefore be checked first whenever two nominally identical calculations give noticeably different numerical results.

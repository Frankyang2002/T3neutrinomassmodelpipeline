# B2: native RGBeta internal coupling tensors versus physical action

## Status and source lineage

The current public project GitHub main was inspected at
`e2a5def5326cc452f668ef593bd81ae13fdff30b` (read-only). Existing
project source blobs are unchanged relative to the earlier B2 inspection;
`b2_rgbeta_internal_action.py` checks those canonical blobs against the
**actual** local Git tree/index/worktree with Git's path-specific clean filters.
It also checks the new native JSON's recorded Git HEAD and the three
source-blobs used by the Wolfram probe.

The upstream RGBeta implementation was read directly at
`aethomsen/RGBeta/Package/FieldsAndCouplings.m` (observed Git blob
`d524fa948fdbcf70ea5a61358b13a78b16a38faa`) and
`Package/GroupsAndIndices.m` (observed blob
`1be60d95d7c6ef6dbde2585c5a811a5284db99d2`). These are **upstream**
source observations; the installed Windows RGBeta build is not cryptographically
identified just by looking at the upstream repository. The live Wolfram probe
is designed to measure its own registry functions.

## Observed RGBeta constructor normalisations (upstream source)

`AddYukawa` specifies for a non-self-conjugate scalar

\[
N_Y=2\times\frac{1}{\sqrt{2}}=\sqrt{2}.
\]

For the restricted source fields
`y1:{Bar[S1], l, F}` and
`y2:{S2, l, Fc}` (Dirac) / `{S2,l,F}` (Majorana), both scalars are
non-self-conjugate. The registered chirality is `Right`, so the native
`Coupling` function multiplies `Matrix[Bar[y_i]]`, and `CouplingBar`
multiplies `Matrix[y_i]`.

`AddFermionMass` uses the internal symmetric-tensor prefactor

\[
N_{MF}=2
\]

with `Chirality->Right` and the same `Coupling` vs `CouplingBar` reversal.
The mass fields are `{F,Fc}` for Dirac or `{F,F}` for Majorana.
No new free-action `1/2` or Dirac↔Majorana physical mass factor is deduced
from this prefactor.

`AddQuartic` specifies an internal `24` multiplied by `1/Sqrt[2]` once
for **each** non-self-conjugate scalar in the four-field vertex. Thus

\[
N_{\text{quartic}}=\frac{24}{(\sqrt{2})^4}=6.
\]

All four fields of `lambdaT3:{H,H,S1,Bar[S2]}` and basic
`lambdaH1/2:{Bar[H],H,Bar[S_i],S_i}` are complex. The source
`lambdaT3` is `SelfConjugate -> False`, whereas portal quartics are
self-conjugate. Consequently the `Coupling` quartic function multiplies
`Matrix[lambda]`, and `CouplingBar` multiplies `Matrix[Bar[lambda]]`.
The native script checks both functions and the exact ratios; the Python
script rejects any mismatch.

**None of the internal multipliers \(\sqrt 2,2,6\) is by itself an
observable physical Lagrangian coefficient or a Matchete–RGBeta coupling
conversion.** They compensate RGBeta's packaging of real/complex fields,
full tensor symmetrization and the conventions of the RG tensor equations.
The physical action must be reconstructed separately, including index
contractions, identical fields, fermionic signs, and Hermitian conjugation.

## What the B2 native script actually measures

For classes A–E, with the representative alpha choices used by the
previous formal probe, it registers the complete restricted RGBeta set
with no matching/RG operations. It reads the live registries:
`$yukawas[y1/y2]`, `$fermionMasses[MF]`,
`$quartics[lambdaT3/lambdaH1/lambdaH2]`. It evaluates each native
`Coupling` and `CouplingBar` function on symbolic flavor indices and
measures the exact ratio to the corresponding bare symbolic `Matrix`
coefficient. All six operator names, fields, chirality and non-self-conjugate
flags are validated in both Wolfram and Python. Python transports the
*internal representation-level ratios* to all 16 alpha models, never a
signed physical coupling dictionary.

A separate small probe records the native private `GroupInvBar` definition
if discoverable and tests its behavior on `eps`, `delS2` and `I*eps`.
The upstream implementation only conjugates explicit complex numbers and
transposes `tGen`; this is **symbolic conjugation behavior**, not a
numerical SU(2) representation-axis mapping.

## Existing Ward result and remaining obstruction

B2's previous analytic witness has exact squared Ward residuals `(3,0,3)`
for the `(2,2)` literal delta–delta topology mixing invariant under the
fundamental/antifundamental interpretation, versus `(0,0,0)` for the
conditional delta–epsilon alternative. It has a similar B-class y1
orientation discrepancy. The abstract RGBeta `eps` and `delS2` remain
unevaluated for numerical component axes. They are not assigned physical
signs in this batch. The 10 basic portal sectors and 4 nonunique Hermitian
nullspaces from the earlier native Matchete basis audit remain unchanged.

Still missing for a physical dictionary:

1. A **numerical** or otherwise independently established conversion from
   RGBeta's barred fundamental/triplet indices and `delS2` to Matchete's
   exact `GetCGTensor` axes, including epsilon orientation.
2. A full source-to-source physical Lagrangian coefficient or Feynman-rule
   witness for `y1`, `y2`, `MF`, and `lambdaT3`, not an internal RGBeta
   tensor prefactor.
3. The B1 Matchete four-spinor to left-Weyl spinor phase/sign interface.

No production files or source interactions are modified. This batch does
not resolve the conditional Ward disagreement and does not perform EFT or
one-loop matching.

**Physical source certificates: 0/16. Not safe to integrate as a physical
Matchete–RGBeta map.**

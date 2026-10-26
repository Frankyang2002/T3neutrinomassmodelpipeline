# Stage 27 — SU(2) triplet component-basis witness

**Scope:** New, read-only source-informed diagnostic for the six `CalibratedTripletCGCandidate` scalar quartic slots. No production RGBeta, Matchete, or running parameters are altered. Not a physical coupling certificate.

## Provenance and prerequisites

- Stage 22 native Matchete `matchete_tensors.json` (50 slots, 16 models) and native RGBeta `rgbeta_native_invariants.json` (50 slots).
- Stage 25 live `CouplingTensorMap` audit, then Stage 26 signed-action/HC report. The command-line entry point calls the Stage 26 audit afresh and compares its **entire** result with the saved Stage 26 JSON; missing or stale inputs cause an error.
- Inspected source: `Lagrangian/model/SU2Invariants.wl`, `Lagrangian/interactions/ScalarPotential.wl`, `RGE/running/rgbeta/T3RGBetaModel.wl`, and `scripts/Stage22RGBetaAlignedGeometry.py` from the project, plus RGBeta v1.2.0 `Package/FieldsAndCouplings.m`.

## Independent A/D calibration

The Matchete native quartic tensors for `A_+0`, `A_-2`, `D_+0`, and `D_-2` contain two barred Higgs doublets and exactly one SU(2) triplet (the other BSM scalar being a singlet). Extract the **native Matchete** 2×2 matrices \(Q_a(i,j)\) and fix the positive norm \(N\):

\[
N^2=\sum_{ij}|Q_a(i,j)|^2=\frac{2}{\sqrt3},\qquad C_a=\frac{Q_a}{N}.
\]

All four models must produce *the same* \(C_a\), without independent rephasings. On the recorded native data:

\[
C_1=\frac{1}{\sqrt2}\begin{pmatrix}1&0\\0&-1\end{pmatrix},\quad
C_2=\frac{1}{\sqrt2}\begin{pmatrix}0&-1\\-1&0\end{pmatrix},\quad
C_3=\frac{-i}{\sqrt2}\begin{pmatrix}1&0\\0&1\end{pmatrix}.
\]

These form an orthonormal basis of \(\mathrm{Sym}^2\mathbb C^2\). They are **read from Matchete component tensors**, not accepted from Stage 22's previously calibrated `CG` constant.

## Independently predict the two E-series tensors

Insert that single \(C\) in the contraction defining `T3RGBetaS2TripleInvariant`:

\[
F_{abc}=\sum_{ijklmn}C_a^{ij}C_b^{kl}C_c^{mn}
\epsilon_{ik}\epsilon_{jm}\epsilon_{ln},\quad
R_{pqbc}=\sum_a C_a^{pq}F_{abc}.
\]

Compare against every *native* Matchete tensor component for `E_+0` and `E_-2`. The verified tensor-only scales \(R=kM\) are:

| Four simple A/D slots | Two E slots |
|---|---|
| \(k=3^{1/4}/\sqrt2\) | \(k=1/\sqrt2\) |

Different \(k\)'s are **tensor normalisations, not running-coupling transformations**. The code rejects a changed native contraction, sign, phase, component, or source scaling.

## Exact SU(2) generator intertwining

Choose an *explicit reference*, **not a claimed native RGBeta basis**:

\[
B=(|11\rangle,\, (|12\rangle+|21\rangle)/\sqrt2,\,|22\rangle),\qquad
U_{ab}=\langle B_b|C_a\rangle.
\]

With \(t^a=\sigma^a/2\) acting on the fundamental, the two barred doublets use \(\bar t^a=-t^{a*}\). In the symmetric-product reference basis:

\[
(T_B^a)_{bc}=\langle B_b|\bar t^a\otimes1+1\otimes\bar t^a|B_c\rangle,
\qquad T_C^a=U T_B^a U^\dagger.
\]

The audit checks all three exact Hermitian generators, the commutators \([T^1,T^2]=iT^3\) and cycles, the spin-1 Casimir \(\sum_a(T^a)^2=2I\), \(U^\dagger U=I\), and \(T_C^aU=UT_B^a\).

The resulting unitary candidate and conjugation diagnostic are

\[
U=\begin{pmatrix}
1/\sqrt2&0&-1/\sqrt2\\
0&-1&0\\
-i/\sqrt2&0&-i/\sqrt2
\end{pmatrix},\qquad UU^T=\operatorname{diag}(1,1,-1).
\]

**Critical distinction:** \(UU^T\) is a basis-conjugation diagnostic. It does **not** establish the action of Matchete `Bar` or RGBeta `Bar` on physical spin-1 components. Nor have we independently extracted Matchete and RGBeta native triplet generator component matrices; the two generator sets above come from an explicit mathematical reference and the Matchete-inferred \(C\). Therefore the generator intertwining alone cannot certify a physical inter-package basis change.

## Files and outputs

New files only:
- `scripts/Stage27TripletIntertwiner.py`
- `tests/python/test_general_t3_stage27_triplet_intertwiner.py`
- `docs/GENERAL_T3_STAGE27_TRIPLET_INTERTWINER.md`

From the repository root in PowerShell (after unpacking):

```powershell
python -m pytest -q tests\python\test_general_t3_stage27_triplet_intertwiner.py
python -m scripts.Stage27TripletIntertwiner
```

Report: `output/general_t3_stage22/stage27_triplet_generator_witness.json`.

## Certification status / next numerical step

This tests whether one **Matchete-derived triplet CG basis** works for all six cases and whether its spin-1 generator representation is consistent. It does **not** measure native RGBeta `S2` component matrices, Matchete native generator components, or the cross-package field conjugation rule. Consequently:

- `global_matchete_rgbeta_cartesian_basis_certified = false`;
- `physical_coupling_mapping_certified = false`.

Before promoting the 46 Stage 23/26 conditional signed maps to physical RG input, independently extract and identify the native component representation and barred-field conjugation on **both sides**. A candidate \(U\) without these checks is not a certificate. No further arbitrary per-block reindexing is justified.

# Stage 22 — aligned SU(2) geometry diagnostic (not coupling normalisation)

This diagnostic consumes **unchanged native evidence** from Matchete and RGBeta:

- `output/general_t3_stage22/matchete_tensors.json`
- `output/general_t3_stage22/rgbeta_native_invariants.json`
- `output/general_t3_stage22/rgbeta_native_declarations.json`

The inputs supplied in the October 26, 2026 run match the SHA-256 digests embedded in `stage22_aligned_geometry_diagnostics.json`. This step changes **no production code**, generates no `rgbeta_native_tensors.json`, and cannot set any physical-coupling certification flag.

## Run

From the repository root (PowerShell):

```powershell
python -m pytest -q tests\python\test_general_t3_stage22_explicit_basis.py tests\python\test_general_t3_stage22_aligned_geometry.py
python -m scripts.Stage22RGBetaAlignedGeometry
```

The report is `output/general_t3_stage22/stage22_aligned_geometry_diagnostics.json`. Its per-block `geometric.r_in_m` matrices express `R = M A` for **component tensors only**. They are not coupling transformations. The physical potential mapping requires verified relative Lagrangian signs, prefactors, conjugation and component-basis identification:

`r = diag(1/pR) A^(-1) diag(pM) m` (direct holomorphic orientation).

### Findings for all 16 ordinary assignments

| Diagnostic status | Slots | Blocks | Meaning |
|---|---:|---:|---|
| `DirectGeometryCandidate` | 32 | 30 | RGBeta invariant components match within the previously stated direct Cartesian/epsilon convention |
| `PseudorealReindexCandidate` | 8 | 8 | Full exact tensor equality after a **local** epsilon reidentification; not an established global field basis |
| `CalibratedTripletCGCandidate` | 6 | 6 | Exact tensors with a Matchete-calibrated, unit-normalised SU(2) triplet CG basis |
| `IncompatibleNativePair` | 4 | 2 | B_-1 and C_-1 `110101` two-coupling blocks; no common I/J reidentification of their four active doublets |
| **Total** | **50** | **46** |  |

Two direct E-model `2×2` blocks contain two coupling slots each, so there are 30 direct blocks for 32 directly matched slots.

### Doublet index reidentification

Use `J = [[0,1],[-1,0]]`, with `J[i_M,i_R]`. In `B_-1` and `C_-1`, applying it to each occurrence of `S1` and `S1dag` in the native RGBeta expression gives an exactly matching one-dimensional Matchete tensor for blocks `000112`, `000202`, `001201`, `020200`. Every resulting `r_in_m` is exactly `[[1]]`.

**Important:** This operation is selected *per field-content block*. It is not a verified universal field redefinition. Applying a universal rotation to other assignments or the two-coupling block can spoil otherwise matched contractions. It therefore provides geometric candidates only.

### Triplet Clebsch–Gordan calibration and prediction

Represent `delS2[SU2L,a,i,j]` by an orthonormal symmetric tensor `C[a,i,j]`:

```
C[1] =  1/sqrt(2) * [[ 1,  0], [ 0, -1]]
C[2] = -1/sqrt(2) * [[ 0,  1], [ 1,  0]]
C[3] = -I/sqrt(2) * [[ 1,  0], [ 0,  1]]
```

Check `sum_ij conjugate(C[a,i,j]) C[b,i,j] = delta_ab` exactly. This basis was calibrated from Matchete `A_+0/021010`; its triplet component coordinate choice is not independently sourced from RGBeta.

For `A_+0`, `A_-2`, `D_+0`, `D_-2` (one triplet with two symmetric H† fields), the exact identity is

`M = (sqrt(2)/3^(1/4)) * R_C`.

For the two E models, evaluate the *existing verified source code*

```
T3RGBetaS2TripleInvariant[a_,b_,c_] =
  C[a,i,j] C[b,k,l] C[c,m,n] epsilon[i,k] epsilon[j,m] epsilon[l,n]
```

and contract `C[a,h1,h2] T3RGBetaS2TripleInvariant[a,b,c]` with no new normalisation. Every `2×2×3×3` component then satisfies

`M_E = sqrt(2) * R_C`

for both `E_+0/021010` and `E_-2/020101`. These are **independent predictions from the calibration** but not proof that RGBeta's abstract triplet basis is numerically this Cartesian basis.

### Structural obstruction in B/C(-1) block `110101`

The native RGBeta entries are:

```
lambdaT3GeneralQ110101I1:
    del[fund, #4, #1] * eps[fund, #2, #3]
lambdaHHdagS1barS2barCross:
    eps[fund, #1, #3] * eps[fund, #2, #4]
```

Their Matchete analogues form an independent, two-dimensional invariant basis. The diagnostic exhaustively enumerates all `2^4 = 16` choices of identity or `J` reindexing on the four distinct doublet species `H`, `Hdag`, `S1dag`, `S2dag`. Some choices bring the *first* RGBeta tensor into Matchete's span; different choices bring the *second* one into that span. **No single choice brings both in simultaneously.** Thus the full-rank `2×2` map cannot be extracted with this reindexing ansatz.

An algebraic candidate replacement for the first contraction *under the local `J(S1dag)` alignment* is `del[fund,#4,#1]*del[fund,#2,#3]`; it is linearly independent of the existing cross contraction in the Matchete component space. This is **not** evidence it is gauge covariant in the actual RGBeta field/Bar conventions, so **do not change any production quartic definition based only on this candidate**. Verify RGBeta's SU(2) transformation rules and the native tensor map first.

## What is / is not verified

**Verified on uploaded evidence:** exact Matchete tensors; actual RGBeta invariant expressions; exact symbolic equality conditional on explicitly described component identifications; triplet CG orthonormality; all 16 model counts and full-rank checks for the 44 aligned blocks; exhaustive binary pseudoreal reindexing obstruction for the other 2 blocks.

**Not verified:** native RGBeta **component** basis relative to Matchete; RGBeta's `Matrix`/`Lam` coupling-tensor `6*Matrix[...]` prefactors relative to physical scalar-potential coefficients; whether one field basis is globally consistent across all quartics; signs from Lagrangian vs scalar potential and complex conjugation. Physical certification therefore remains **false for all 50 slots**.

The next production-relevant step is to check the two RGBeta `110101` quartic tensors at the group/field-conjugation level, then obtain native RGBeta coupling-tensor components and signed potential prefactors. Do not populate the production coefficient mapper from this diagnostic alone.

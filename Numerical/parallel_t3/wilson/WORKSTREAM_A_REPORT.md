# Workstream A batch implementation report — 2026-11-09

**Public repository inspected, read-only:**
`Frankyang2002/T3neutrinomassmodelpipeline`, `main` HEAD
`a944e0030b8ef24b23342d99e0d8f33cef3cffd5`.
The user's Windows working tree and uncommitted data were not visible.
No GitHub operations with write permissions were used. No existing shared
sources (including `ScalarPotential.wl`) were replaced.

## Implemented and exercised

- All 16 T3 bases: 45 complex SM left-Weyl components, all `H,S1,S2` real
  scalar components (representation-dependent size), SU3/SU2/U1 generators,
  full three-family complex `C[i,j,a,b]` with both independent pair symmetries.
- Full Weinberg SU2-invariant polynomial embedded into and projected from the
  very same `C` tensor, including arbitrary nonzero `C5(MF)`.
- One coupled matrix-free RGE in `ln(mu)` with time-dependent B-provided
  `Y`, `lambda4` and `g2,gY,g3`. Both C5 and LLSS mix in both directions;
  non-projected tensor components remain dynamical. RK2 and DOP853 modes.
- Three unambiguously labelled ADM variants; formula ambiguity remains and
  **no physical anomalous-dimension formula is certified**.
- Full canonical (complex) tensor NPZ plus JSON metadata and input/output
  SHA256 hashes; strict B interface, field coordinates, corrected A(-4)
  `ScalarPotential.wl`, stale provenance and threshold ordering guards.
- Conditional Matchete↔RGBeta quartic polynomial fitting with grouping of
  reused symbols, matrix ranks and fit/holdout residuals, and a separate
  native RGBeta component probe that fails closed on unsupported CG functions.
- Explicit real-orthogonal scalar and unitary Weyl transform helpers including
  quartic, Yukawa, mass and Wilson tensors; unproved matrices are never
  inferred from a shared coupling name.
- Gauge-singlet synthetic `LLSS`/scalar-portal polynomials via **SU2 Ward
  nullspaces** for all 16 models. Not a Matchete-native phase/normalization.

## Tests run here

`python -m unittest discover -s tests -p test_t3_parallel_a.py -q`:
**12 passed** (9–10 seconds). Includes generic complex unitary-basis covariance
comparison: six-term and corrected-crossed candidates pass; printed
unstarred-C seven-term fails. Includes all16 generator/portal Ward tests,
pair symmetries, Weinberg projection, real-linearity, nonzero C5 hard,
LLSS↔C5 mixing, finite evolution, source/basis blockers, disk-backed
B-style formal boundary and rejected stale source.

`python -m scripts.RunT3ParallelA --synthetic-all16 --steps 1 --solver RK2
--out /mnt/data/_t3a_all16_full_llss_v2`: **16/16 Candidate_Qualitative**,
all with nonzero LLSS and C5. No real fermion boundaries were used.

## Blocked physics / NOT implemented as certified

- `Blocked_ADM_PhysicalNormalization`: Eq. (4.85) unstarred-C crossed term
  in `2512.16890v1` violates complex Weyl basis covariance. `C*` and the
  transposed fermion Casimir repair this *necessary* condition but do not
  establish the actual crossed loop counterterm; separate Fonseca Eq. (22)
  supplies the six-term gaugeless reference only. Published versions are not
  analytically reconciled at loop level. Do not use these results as physical.
- `Blocked_GammaCCWeylSignedNormalization`: no native Matchete four-spinor
  -> canonical two-component Weyl sign/factor certificate; code requires
  explicitly supplied number or opt-in formal identity, never derives it.
- `Blocked_MatcheteRGBetaCompleteQuarticMap`: RGBeta's saved Stage22 file is
  declaration metadata, not normalized numerical tensors. The new Wolfram
  probe has **not** been run here. Neither the ordinary nor the 50 extra
  polynomial slots are fully mapped and source-normalized. `lambdaH2Inv1`
  requires `I` inside `PlusHc`, `Inv2` must survive; the corrected source
  is only read/guarded, not edited.
- `Blocked_NativeOneLoopFHardMatching`: formal boundary and optional hard
  `C5(MF)` are supplied by B; no invented physical F matching.
- `Blocked_FullMassAndScalarEigenbasisForSpecialModels`: physical
  neutral-self-bilinear phases and B_-1/C_-1 scalar eigenbases require B/C.
- Closure argument: other **physical** dim5 classes (`phi^5`, `phi F^2`,
  `psi^2 F`) appear forbidden by the unbroken Z2 and SM chiral/gauge charges,
  with no F fields in EFT1. A comprehensive native symmetry enumeration and
  EOM/redundant off-shell Green-basis analysis were **not** performed.
- No real Workstream B trajectory, no one-loop scalar threshold from C,
  no scalar matching, no SMEFT running and no neutrino mass prediction.

## Immediate integration actions

1. B supplies freshly regenerated A(-4) or all16 formal fermion LLSS sources,
   an independently justified optional hard `C5(MF)`, and canonical full
   running tensors/manifest with hashes matching the real working tree.
2. Run `RunT3ParallelAExportRGBetaComponents.wl` on Windows and resolve
   native normalization/CG frame. Then use `RunT3ParallelATensorDictionary`
   with **fully normalized polynomial records**, not declaration strings.
3. Establish spinor sign/factor and the complete one-loop ADM from native
   Feynman conventions or published corrected loop expressions before
   removing any `physical_certified:false` flag.
4. C consumes `wilson_at_MS.npz:C_full` and `C5_EFT1_MS` at one *joint*
   scalar boundary, without re-adding leading logs.

**Outputs remain conditional even when all numerical tests succeed.**

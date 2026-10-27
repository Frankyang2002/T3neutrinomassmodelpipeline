# Native Yukawa2 SU(2) recoupling (experimental)

B_-1 and C_-1 pass the free Lagrangian and first ten separated baseline
contractions, but y2|CG:T3Y2CG aborts `CheckLagrangian` with
`DefineCG::notInvariant`. The source y2 interaction contains `Bar[S2]`;
the mass eigenstate X/Y is an unconjugated doublet after the epsilon
intertwiner. Leaving original `T3Y2CG` times an epsilon bridge in the
transformed expression does not produce a native Matchete CG for those
external indices.

The new `T3MEBuildNativeYukawa2` computes, using the source's current Matchete
SU(2) basis, the old tensor and new tensor, contracts epsilon on the old
scalar slot, and checks every component of `T_old eps = alphaCG T_native`.
The sign alphaCG is **calculated**; in a temporary native Matchete tensor
kernel it was -1 for the B singlet-fermion variant and +1 for the C
triplet-fermion variant. The formula contains its own explicit - sign
from `Bar[S2] = -epsilon (epsilon Bar[S2])` as in the mass-rotation
convention. The new tensor replaces all source y2 terms, including their
Hermitian-conjugate contributions, after lossless algebraic source separation.
The y1 source tensor needs no epsilon transport and remains in the generic
rotation; it is checked by the existing per-operator validator.

## Validation

- Native Matchete `InvariantTensors`: old and new spaces one-dimensional
  for B (dF=1) and C (dF=3); entrywise tensor identity verified.
- Offline Python and source structure tests check the output certificate,
  source completeness, dimensions, and strict failure handling.
- **Full Matchete CheckLagrangian/HeavyMassBasis not yet live-verified**.
- Does not perform physical matching or certify RGBeta coupling normalization.

## Project paths

- `Lagrangian/RunGeneralT3MassEigenMatchete.wl` (replacement)
- `scripts/RunGeneralT3MassEigenMatchete.py` (replacement)
- `tests/python/test_general_t3_mass_eigen_matchete.py` (replacement)
- `docs/T3_MATCHETE_NATIVE_YUKAWA2_CG_FIX.md` (new)

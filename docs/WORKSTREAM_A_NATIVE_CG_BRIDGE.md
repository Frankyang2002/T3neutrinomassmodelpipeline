# Workstream A — Live registered CG baseline evidence

## Scope

This opt-in audit extracts the **actual registered** `T3Y1CG`, `T3Y2CG`, and
`T3MixCG` Clebsch–Gordan tensors from installed Matchete using its public
`GetCGTensor` accessor, following a full Lagrangian build. All three tensors
must occur in the assembled interaction expression. The audit covers 16
ordinary T3 points, i.e. **48 live registered baseline tensors**. The scalar
mass, final extra-quartic and Z2 tests were performed by previous stages and
are required as inputs, not repeated here.

The imported field conventions are the unchanged repository conventions:

- `Yukawa1`: conjugated lepton, conjugated heavy fermion, `S1`.
- `Yukawa2`: conjugated lepton, heavy fermion, conjugated `S2`.
- `T3Mix`: `H H S1 Bar[S2]`, with the first two positions Bose-symmetric.
- A gauge singlet does not create an SU(2) index on the native CG array.

For the first two factors of `T3Mix` the audit checks

\[ C_{ij\cdots}=C_{ji\cdots}. \]

Across alpha values of the same SU(2)-dimension class the audit tests that
registered arrays agree up to an overall **nonzero** complex scalar:

\[ C^{(a)}_{i_1\cdots i_n}=r\,C^{(b)}_{i_1\cdots i_n},\quad r\ne0. \]

The B/C topology CG is separately compared because B and C share the same
`(dS1,dS2)=(2,2)`, despite different fermion representations. This is an
exact arithmetic check, not a comparison of approximate norms.

## Existing conditional normalization results

`output/general_t3_stage22/stage23_quartic_prefactor_candidates.json` has
46 conditional normalization blocks containing 50 additional quartic slots.
The prior exact-equivalence stage checked their matrix identities. The new
stage requires the source evidence to retain its **conditional** status and
inspects the actual RGBeta baseline operator definitions in
`RGE/running/rgbeta/T3RGBetaModel.wl`.

For the physical comparison of two quartic interaction tensors, one needs a
single common field-component convention, including the conjugate SU(2)
doublet/triplet identification and the complex-to-real-scalar convention.
The existing Stage27 and Stage28 reports explicitly **do not provide a
numerical native RGBeta triplet generator basis or certified Bar-field
intertwiner**. A matrix inverse of assumed coefficient maps cannot supply
that missing physical statement.

Consequently, successful execution **cannot certify**:

- a physical signed Matchete-to-RGBeta coupling dictionary;
- the relative native Yukawa, T3 topology, and scalar CG normalizations;
- symbolic field substitution equality for entire UV Lagrangians in both
  package conventions;
- explicit Z2 registry support in Matchete; or
- physical matching or loop-level equality.

It instead reports a distinct milestone: live Matchete baseline tensors and
their group-theoretic internal consistency, while preserving the physical
common-basis obstruction as an explicit unresolved prerequisite.

## Files and commands

New project files:

- `Lagrangian/RunT3UVNativeBridge.wl`: live Matchete source capture for one model.
- `Lagrangian/RunT3UVNativeBridge.py`: 16-model orchestrator, exact arithmetic
  validation and report generation.
- `docs/WORKSTREAM_A_NATIVE_CG_BRIDGE.md`: this note.

From the repository root, in the existing `.venv`:

```powershell
python Lagrangian\RunT3UVNativeBridge.py --self-test
python Lagrangian\RunT3UVNativeBridge.py --all
```

The Python runner sets all three required experimental switches to `1` and
the additional-quadratics switch to `0` for each fresh Wolfram subprocess.
It does not modify inherited user environment settings.

To re-evaluate **already generated** native JSONs without launching Wolfram:

```powershell
python Lagrangian\RunT3UVNativeBridge.py --from-existing
```

A second run that would replace existing native records requires the explicit
`--overwrite` flag. Back up the output directory before using it.

The generated files are contained in
`Reports/output/workstream_a/native_cg_bridge/`:

- `native/<MODEL>.json`, `native/<MODEL>.log`
- `native_cg_bridge.json`
- `native_cg_bridge.md`

No production RGE, EFT matching, or numerical files are modified. No changes
to the GitHub repository are made.

## Expected conditional success

```
LiveMatcheteBaselineCGsPassedRGBetaPhysicalBridgeUncertified
```

This status means all 48 actual Matchete CG records have been extracted and
validated and all 36 cross-model proportionality comparisons succeeded. It
**does not** mean the physical cross-package normalization is closed.

# T3 mixed scalar native Matchete: internal objects compatibility fix

## Trigger

Both B_-1 and C_-1 exited at `original-interaction-assembly-started` with
`AdditionalQuarticCouplingsAbsentInOriginalSource`, listing the same seven
`lambdaT3ExtraQ...` labels. Live logs already showed all seven independent
Matchete quartic candidates generated.

## Verified cause

Matchete's user-facing shortcuts are converted to internal expressions.
Its documented examples show `phi[]` represented by `Field[phi,Scalar,{},{}]`
and `y[]` represented by `Coupling[y,{},0]`. The old check looked for
`lambdaT3ExtraQ...[___]` as an applied head and therefore failed when the
label appeared as an argument to `Coupling`. The old field transformation also
matched only unevaluated `NewScalar1[i]`/`NewScalar2[i]` shortcuts.

## Changes

Full replacement `Lagrangian/RunGeneralT3MassEigenMatchete.wl`:

- Coupling-presence check searches for the exact symbol *at any depth* in the
  already evaluated Matchete expression. The kernel self-test covers both
  shorthand application and `Coupling[label,{},0]`, plus absent-label controls.
- Old-field detection searches for the exact `NewScalar1` and `NewScalar2`
  labels, including their occurrence inside internal `Field[...]` objects.
- Four scalar rotation rules now match `Field[NewScalar1,...]`,
  `Bar[Field[NewScalar1,...]]` and the analogous `NewScalar2` expressions.
- Representative Matchete field/coupling signatures are printed before the
  SU(2) bridge. If the installed version has a different conjugation object
  or index structure, the runner fails closed and reports these structures.
- Existing UV/EFT1 numerical data, mixing matrix, SU(2) bridge checks,
  physical-matching flags and source-provenance checks are not changed.

## Scope and limitations

Offline tests and a kernel-only symbolic mock check pass. No local Matchete
installation was available for a live run. This change does NOT certify the
new heavy-field mass basis, the Matchete–RGBeta normalization map or physical
Wilson-coefficient matching.

If `OriginalScalarFieldsRemainAfterRotation` occurs, inspect the printed
`Raw original Matchete signatures`; they will identify whether the installed
Matchete encodes conjugate fields differently from the expected representation.
Do not disable source-inventory or heavy-mass-basis checks.

## Installation and execution (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Matchete_Internal_Objects_Fix.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
```

Generated run logs and progress files remain under
`output/general_t3_mass_eigen_matchete/`.

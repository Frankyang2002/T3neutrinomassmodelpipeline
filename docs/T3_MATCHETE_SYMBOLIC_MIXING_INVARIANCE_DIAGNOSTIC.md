# T3 B_-1 / C_-1: exact symbolic mixing in native Matchete check

## Trigger

After fixing the typed bridge indices, both models reached `CheckLagrangian`
but aborted with `DefineCG::notInvariant: The tensor is not invariant under group
transformations.` The failure was observed for the assembled transformed
Lagrangian; it has **not** established which vertex, CG orientation, or numerical
operation is at fault.

## Change and rationale

The prior Wolfram runner directly inserted **machine-precision** values of the
2x2 numerical scalar mixing matrix `U` into Matchete expressions. Matchete
constructs and checks group-invariant Clebsch-Gordan tensors exactly. Floating
coefficients might interfere with that exact CG reduction; the existing
messages alone do **not** prove that this is the cause.

This experimental runner now defines four neutral complex Matchete couplings

    T3MEU11[], T3MEU12[], T3MEU21[], T3MEU22[]

and transports interactions algebraically with them and their `Bar[...]`
conjugates. The numerical `U` is *still* checked for unitarity and scalar mass
squared diagonalisation, using the unchanged input JSON. A new result field
`rotation_spurion_values_from_numerical_U` contains the exact double-precision
boundary values of all four symbolic couplings. Python checks equality against
the input U on any purported successful Matchete run.

This avoids approximate coefficients in the gauge-CG calculation. **The four
couplings are formally independent inside Matchete**: a successful structural
`CheckLagrangian` alone would not prove they obey unitarity identities or
certify the Matchete/RGBeta component basis map. The separately verified
numerical map supplies those constraints at the benchmark point only.

## Failure isolation

If the full `CheckLagrangian` still fails, the Wolfram runner first tests
the combined **11-vertex baseline**, then tests each of the seven additional
quartic operators separately. It uses already registered source CGs; it does
not regenerate them. It prints the first failing group or quartic and records
a diagnostic in
`isolated_source_failure`. This is an **investigation tool**, never a way to
turn a failed full check into a pass. If even the first individual term fails,
the bridge definition and conjugate-index orientation need inspection against
the installed Matchete group convention, not further numeric fixes.

No Wilson matching is done, no EFT1 data regenerated, no repository GitHub
content changed, and no cross-package physical tensor normalization certified.

## Files

Full replacements:

- `Lagrangian/RunGeneralT3MassEigenMatchete.wl`
- `scripts/RunGeneralT3MassEigenMatchete.py`
- `tests/python/test_general_t3_mass_eigen_matchete.py`

New document:

- `docs/T3_MATCHETE_SYMBOLIC_MIXING_INVARIANCE_DIAGNOSTIC.md`

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Matchete_Symbolic_Mixing_CG_Check.zip' -DestinationPath . -Force

python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py

python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
python -m scripts.RunGeneralT3MassEigenMatchete --all
```

If there is still a failure:

```powershell
Get-Content .\output\general_t3_mass_eigen_matchete\B_m1_final_matchete.log -Tail 100
Get-Content .\output\general_t3_mass_eigen_matchete\B_m1_final_matchete.json
```

Expected offline test status for this package: 82 tests passing in the
packaging environment. **Live installed Matchete remains untested by the author
of this package.**

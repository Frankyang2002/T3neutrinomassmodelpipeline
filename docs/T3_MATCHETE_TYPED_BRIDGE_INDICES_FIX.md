# Typed pseudoreal-bridge indices: B_-1 / C_-1

## Observed live failure

Both transformed B_-1 and C_-1 Matchete Lagrangians reached the final check,
but the installed Matchete aborted with:

```
DefineCG::invldIndType: j$32088 is not defined as a representation.
```

The original scalar fields have already evaluated into Matchete's internal
form `Field[NewScalar2,Scalar,{Index[d$$3,T3BSMd2]},{}]`.
The bridge interaction inserted a fresh bare `j$...` in
`CG[T3MEBridgeBB|T3MEBridgeUU,{i,j}]`. Matchete cannot infer a
representation for that index when checking the result.

## Change

The two CG slots now contain indices of the same explicit type:

```
CG[T3MEBridgeBB, {i, Index[j,T3BSMd2]}] Bar[NewScalarX[j]]
CG[T3MEBridgeUU, {i, Index[j,T3BSMd2]}] NewScalarX[j]
```

The `j` in `NewScalarX[j]` is deliberately a bare *short-form field index*:
Matchete converts that field shortcut to an internal `Field[..., {Index[j,T3BSMd2]}, ...]`.
The bridge and new field thus contract the same typed dummy index.
Identical changes are made for `NewScalarY`. Original first index `i` is
already a typed `Index[...]` from the source `Field` object.

An additional preflight rejects bridge contractions with untyped or wrongly
typed CG index slots before the ordinary and detailed `CheckLagrangian` calls.
It does not certify numerical component-normalization alignment.

## Install / run

Extract the ZIP into the repository root, replacing the existing Wolfram runner
and test file; no production source or numerical outputs are changed.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Matchete_Typed_Bridge_Indices_Fix.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
```

If a new live check fails, the progress JSON and execution log in
`output/general_t3_mass_eigen_matchete/` identify the stage and Matchete messages.

## Validation scope

The source and offline tests verify the typed-index code path. The Matchete
kernel itself was **not** available in the packaging environment. A native
`HeavyMassBasis` pass remains unverified until both live runs succeed.
No scalar-threshold Wilson matching has been performed.

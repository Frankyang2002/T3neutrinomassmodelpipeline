# Matchete SM field-label collision fix

## Verified live failure

The B_-1 run reached `baseline-build-started`, then `LoadModel["SM"]`
raised `DefineField::FieldLabel` with the **numerical 2x2 U matrix** as the
field label. The original Wolfram runner assigned `u = ...` before loading
Matchete's Standard Model, overwriting Matchete's standard up-quark label `u`.
This is a code namespace collision, not a scalar mass-basis physics failure.

## Correction

In `Lagrangian/RunGeneralT3MassEigenMatchete.wl`, all 15 uses of the global
`u` matrix variable are renamed to `$T3MEUnitaryMatrix`. The JSON input key
`unitary_U`, formulas for unitarity/diagonalisation and the full Yukawa/quartic
field substitution are otherwise unchanged. Existing numerical outputs remain
valid. The diagnostic checkpoint runner and Python provenance guard are kept.

## Installation and exact commands (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_MassEigen_Matchete_SM_u_Symbol_Fix.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
python -m scripts.RunGeneralT3MassEigenMatchete --all
```

If another live error appears, inspect `output/general_t3_mass_eigen_matchete/B_m1_final_matchete.progress.json`
and the matching `.log`; do not interpret the offline tests as a live
Matchete validation. No Wilson matching is performed by this fix.

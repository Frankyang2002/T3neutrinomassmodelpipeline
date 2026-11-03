# General T3 sequential EFT: tensor-axis correction

This package replaces **complete files**, not diffs. Extract at the root of
`T3neutrinomassmodelpipeline` after installing the earlier general sequential
EFT leading-log package. It does not change any GitHub files.

## Bug fixed

The pre-existing numerical `LLSS` and scalar-matching interface expects

`C[p, q, i, j, a, b]` with shape `(3, 3, 2, 2, nScalar, nScalar)`.

The dynamic nonmixed UV fermion-threshold adapter instead allocated axes
`(2, 3, 2, 3, nScalar, nScalar)`, causing the scalar matcher to report
`Full amplitude must have shape (3,3,2,2,4,4)`.

The B(-1)/C(-1) provisional source reweighting also incorrectly transposed the
reference tensor with interleaved SU2/flavour axes, yielding the NumPy einsum
`Size of label 'q' ... (2) ... (3)` error.

Both computations now use the same flavor-first physical tensor layout. The
nonmixed boundary has direct combined lepton/scalar exchange checks, and the
mixed basis reweighting preserves the same symmetries. No renormalisable beta
functions or matching formulas are changed.

## Files (complete replacements)

- `Numerical/running/GeneralT3DynamicFermionWilson.py`
- `tests/python/test_general_sequential_leading_log.py`

## PowerShell (repository root)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_General_Sequential_TensorAxis_Fix.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_general_sequential_leading_log.py
python -X utf8 -m scripts.RunT3GeneralSequentialFitAblation `
  --source-run thesis16_general_c5_full_v1 `
  --run-id general_seq_03 `
  --spinor-factor 1
```

Use a **new** run ID: the existing script refuses to overwrite `general_seq_02`.
The `thesis16_general_c5_full_v1` source contains 15 model fit outputs; D(+0)
was blocked in that saved campaign. To cover all 16, first combine the five
continuation folders from `continue_general16_v1` into a merged source folder
as previously described, then pass `--source-run thesis16_merged_v1`.

## Validation and scope

- 11/11 local algebraic, axis-contract, and interface tests passed.
- Python sources compiled successfully.
- The full numerical pipeline was not executed here (the user's fitted Matchete
  and RGBeta data are not mounted in this environment).
- This changes only tensor indexing. Existing provisional spinor normalization,
  HHSS basis identification and one-scale matching-anchor assumptions remain.

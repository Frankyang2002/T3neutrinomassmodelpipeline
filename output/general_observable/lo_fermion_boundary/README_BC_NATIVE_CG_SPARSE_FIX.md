# B/C native CG product audit: Matchete SparseArray CSR fix

**Project directory:** `output/general_observable/lo_fermion_boundary/`

This package contains **full replacement files** for:

- `check_bc_native_cg_products.py`
- `test_bc_native_cg_products.py`

and this README. It does **not** modify `native_bc_yukawa_cg_probe.wl`, `bc_majorana_exchange_audit.py`, matching code, or the 16-model boundary release.

## Verified failure and repair

The Mathematica extractors produced successful fresh-kernel B_-1 and C_-1 component reports. The downstream Python audit then failed with `ValueError: Sparse coordinate outside dimensions` while decoding the archived *matched* CG registry.

Matchete's `SparseArray[..., {1, {{pointers}, {coordinates}}, {values}}]` stores the first tensor index in compressed-row pointers and **only the remaining rank-1 indices** in the explicit coordinate list. The old parser incorrectly treated those suffix coordinates as full-rank coordinates.

For example, the B_-1 identity CG is stored as

```
SparseArray[Automatic, {2, 2}, 0, {1, {{0, 1, 2}, {{1}, {2}}}, {1, 1}}]
```

`{0, 1, 2}` encodes the two rows. The suffixes `{1}` and `{2}` encode columns. In the C_-1 four-index tensor, suffix coordinate triples are similarly distributed between two first-index rows via `{0, 3, 6}`.

The new `parse_matchete_sparse_array` implements this representation. It rejects row-pointer inconsistencies, out-of-range suffixes, duplicate entries, unsupported Wolfram syntax, and non-listed numeric values. Only the exact rational values present in the pinned matched CG registries are accepted.

## Run (PowerShell)

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
$dir = 'output\general_observable\lo_fermion_boundary'
$zip = "$HOME\Downloads\T3_BC_Native_CG_Sparse_Parser_Fix.zip"
if (-not (Test-Path -LiteralPath $zip)) { throw "ZIP not found: $zip" }
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

python -X utf8 "$dir\test_bc_native_cg_products.py"
if ($LASTEXITCODE -ne 0) { throw 'Regression tests failed' }

$folder = "$dir\normalization_probes"
$b = "$folder\bc_native_yukawa_cg_B_m1_d8cd56bfeab640129bb497178ae003e6.json"
$c = "$folder\bc_native_yukawa_cg_C_m1_63a28f8afbfd4883834c2aa490356833.json"
if (-not (Test-Path -LiteralPath $b) -or -not (Test-Path -LiteralPath $c)) {
    throw 'Successful fresh-kernel Matchete CG reports not found'
}

python -X utf8 "$dir\check_bc_native_cg_products.py" --report-b "$b" --report-c "$c" --check-only
if ($LASTEXITCODE -ne 0) { throw 'Native B/C CG product comparison failed' }
python -X utf8 "$dir\check_bc_native_cg_products.py" --report-b "$b" --report-c "$c"
if ($LASTEXITCODE -ne 0) { throw 'Native B/C CG product comparison failed' }
```

Mathematica **does not need to be rerun**. Expect **18 tests** from the replacement suite. `--check-only` computes the full comparison without writing output; the second invocation writes a deterministic new JSON in `normalization_probes/`. If a downstream numerical identity genuinely differs, the script fails closed and reports the offending tensor instead of certifying it.

## Limitations

This tests numerical SU(2) CG contractions and pinned native B/C source provenance. A successful result is **not** a fully connected native Majorana Feynman amplitude or proof of the Matchete `GammaCC`-to-master-Weyl phase. All such certification flags remain false. The EFT sequence remains UV -> integrate out all fermions -> EFT1 with both scalars -> integrate out all scalars -> SMEFT.

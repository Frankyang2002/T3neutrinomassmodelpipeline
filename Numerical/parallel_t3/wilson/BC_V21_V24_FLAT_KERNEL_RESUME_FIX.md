# Workstream A v24 — v19 flat SVD null-vector contract (B(-1), C(-1))

## Verified source problem

`Numerical/parallel_t3/wilson/complex_scalar_bridge_v19.py` computes

```python
_, _, vh = np.linalg.svd(b, full_matrices=True)
null = vh[-1].conj()
null = null / np.linalg.norm(null)
...
'kernel_direction_normalized_complex': {
    'real': null.real.tolist(), 'imag': null.imag.tolist()
}
```

As `b` has shape `(6,7)`, the saved `null` has shape `(7,)`. In v21,
`validate_special_matrix` incorrectly required `(1,7)` using `finite_matrix(...,(1,7))`.
The synthetic v21 test fixture used a one-row array, hiding the mismatch.

## Correction

- `finite_v19_kernel` accepts the actual saved `(7,)` real/imag arrays and
  a single-row `(1,7)` legacy synthetic fixture. It rejects all other
  shapes, including `(7,1)`, broadcasting between shapes, missing components,
  NaN/Inf, and nonnumeric data.
- The original `rank(B)=6`, `B@n=0`, and `||n||=1` checks are unchanged.
- `finite_matrix` now checks the separate `real` and `imag` shapes before
  addition; it cannot accidentally broadcast a malformed matrix.
- The primary synthetic fixture now uses v19's actual `(7,)` serialization.
  New negative controls cover invalid shapes and a realistic non-diagonal
  B/C two-coupling source sector after column reordering.

No change to the scalar fields, SU(2) generators, `epsilon` transport,
Matchete/RGBeta coupling coefficients, EFT matching, v19 output, or the v23
source-order permutation. The status remains conditional and no physical
Wilson evolution is authorized.

## Install and run (PowerShell)

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
$zip = "$env:USERPROFILE\Downloads\t3_workstream_A_bc_v24_flat_kernel_fix.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

.\.venv\Scripts\python.exe -m unittest discover `
    -s tests -p "test_t3_parallel_a_bc_mass_bridge_v21_v22.py" -v
if ($LASTEXITCODE -ne 0) { throw "v24 regression failed" }

powershell -ExecutionPolicy Bypass -File `
    scripts\RunT3ParallelABCMassEigenBatchV21V22.ps1
if ($LASTEXITCODE -ne 0) { throw "B/C conditional mass-eigen bridge failed" }
```

Source evidence must already exist from v18–v20 and the experimental
`output/general_t3_mass_basis/B_m1_final_mass_basis.json` and
`C_m1_final_mass_basis.json` files. The existing batch runner performs
source-hash preflight, named-column permutation check, v21 action check,
and v22 readiness gate. No prior witness should be overwritten.

Expected existing output locations after passing:

```
output/parallel_t3/wilson/bc_native_mass_eigen_complex_conditional_v21.json
output/parallel_t3/wilson/bc_quartic_mass_readiness_conditional_v22.json
```

## Validation classification

- Exact producer shape: inspected v19 implementation (read-only source).
- Six v24-specific regressions: passed locally with an isolated, lightweight
  dependency harness using the **actual full replacement v21 module**.
- Local Windows integration with saved v18–v20 native witnesses: not executed
  in this container; user must run the above command.
- The existing v21/v22 suite's physics-dependent tests were not re-executed
  against the user's installed Matchete/RGBeta environment here.

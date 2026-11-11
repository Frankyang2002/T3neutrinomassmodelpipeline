# Worker B: constrained fixed-real-scalar basis search

Project paths:
- `scripts/RunT3WorkerBBasisSearch.py`
- `tests/test_t3_worker_b_basis_search.py`
- `docs/WORKER_B_BASIS_SEARCH.md`

Requires Batch3B cache and previous Worker B batch reconciliation output.
This is a **diagnostic** scan, not a certified Matchete/RGBeta coupling map.

Tests 384 species-preserving transformations per model: complex component
permutations within S1/S2, uniform species phases from {1,i,-1,-i}, and
species-wide complex conjugation. The Higgs basis stays fixed. Does not test
arbitrary continuous SU(2) triplet transformations or independent phase per
component; a nonzero minimum is not proof of an incorrect physical tensor.

For lambdaT3 real/imag columns, reports relative errors before and after
allowing an unconstrained **real 2x2 least-squares coupling map**. The 2x2 fit
is only a diagnostic of column-span overlap and is never applied to the
physics tensor. Also compares all response-direction spans after the best
scalar rotation and rechecks lambdaH2 sign.

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
python -X utf8 -m unittest discover -s tests -p "test_t3_worker_b_basis_search.py" -v
python -X utf8 -m scripts.RunT3WorkerBBasisSearch `
  --models A_-4,D_-2 `
  --response-root "output\audits\batch3b_quartic_fixed_20261111_234059" `
  --reconcile-dir "output\audits\workerb_batch_reconcile_v1" `
  --out "output\audits\workerb_basis_search_v1\report.json"
```

The input reconciliation report must match its original NPZ SHA256 hash.
No existing files are modified. This does not perform joint scalar matching.

# Batch 15 — Matchete T3MixCG SU(2) triplet basis correction

## What was verified

The **archived** Matchete `T3MixCG` has six nonzero entries (common factor `3^(-1/4)`):

- C[1,1,1]=+a, C[1,1,3]=+ia
- C[1,2,2]=-a, C[2,1,2]=-a
- C[2,2,1]=-a, C[2,2,3]=+ia

The SU(2) triplet generators inferred by imposing invariance of **that CG** are
pure-imaginary, Hermitian Cartesian generators, rather than the highest-weight
spherical spin-1 generators used in the Batch 13 RGBeta tensors.

This explains exactly the *SU2* tree-operator residuals observed in Batch 14:
`[1.847759065, 1.847759065, 2]`.  With the CG-derived generators the
same source-level toy tensor has residuals O(10^-16).  The triplet generators
obey `[T1,T2]=i T3`, etc., and an explicit unitary map to highest-weight
spin-1 is constructed from their T3 eigenvectors and raising matrix elements.

## Source-aligned change

The source triplet `z_Matchete = U z_spherical`. The interleaved-real
8x8 map is inserted into:

1. The complete 3-flavour F-tree LLSS tensor before evaluating the generic ADM;
2. **ONLY** the native CG-derived mixing-quartic contribution `lambdaT3`;
3. The archived native-quartic -> RGBeta direction-comparison audit.

Baseline RGBeta portal/self reference columns are **already** spherical and are
not subjected to the new source transformation. For D(-2), the existing
independently specified **neutral singlet** real-split mass rotation is
composed with the triplet map only for the Matchete-sourced terms. The D extra
HHSS and scalar-only **reference** tensors are not double-rotated.

The updated Batch13 run fails closed if the full F-tree tensor still violates
SU2/U1/SU3 Ward invariance; it reports scalar-potential and ADM Ward residuals
without assuming they vanish. The updated Batch14 audit reads the **new**
Batch13 directory and checks the corrected tensor with spherical generators.

## Limitations

This fixes a **verified basis mistake**, not the three original physics gates:

- It does not derive the Matchete GammaCC -> Weyl factor/sign;
- It does not produce full finite MSbar counterterms for all first-threshold
  operators or a complete independent two-threshold matching proof;
- It does not certify the full gauge-complete anomalous dimension or the
  absolute Matchete/RGBeta quartic normalization, especially the D extras.

**No source file on GitHub is modified. No previous outputs are overwritten.**

## Full replacement files in this ZIP

- `Numerical/parallel_t3/complete/native_triplet_basis_batch15.py` (new)
- `Numerical/parallel_t3/complete/two_model_master_wilson_batch13.py` (full replacement)
- `Numerical/parallel_t3/complete/two_model_physical_closure_batch14.py` (full replacement)
- `tests/test_t3_native_triplet_basis_batch15.py` (new)
- `README_T3_SOURCE_TRIPLET_CG_BASIS_BATCH15.md` (new)

## Windows PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\t3_source_triplet_cg_basis_batch15.zip" -DestinationPath . -Force
.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_native_triplet_basis_batch15.py tests\test_t3_master_wilson_batch13.py tests\test_t3_physical_closure_batch14.py
.\.venv\Scripts\python.exe -X utf8 -m scripts.RunT3MasterWilsonBatch13 --all --output output/parallel_t3/thesis_master_wilson_batch13_cgfix_v3 --debug
.\.venv\Scripts\python.exe -X utf8 -m scripts.RunT3PhysicalClosureBatch14 --all --generator-details --output output/parallel_t3/thesis_physical_closure_batch14_cgfix_v2 --debug
```

The final `summary.json` includes separate `tree_Weinberg_parent`,
`candidate_ADM_delta`, and `scalar_potential` Ward maxima. A Ward check is
**necessary but not sufficient** for physical correctness, even if all vanish.

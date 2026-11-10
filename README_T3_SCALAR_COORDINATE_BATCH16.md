# Batch 16 — Matchete ↔ RGBeta scalar-quartic coordinate investigation

## Purpose and strict scope

This is a **local-only, read-only source audit plus conditional EFT1 rerun** for
`A_-4` and `D_-2`. The existing Batch 15 triplet CG basis is retained.

The Batch 14 one-to-one comparisons do not determine a scalar invariant map.
Batch 16 instead constructs fourth-derivative `lambda_abcd` tensors in one
common real scalar basis and computes simultaneous native-to-reference linear
coordinate transformations within each gauge-invariant field-content sector.

**No GitHub changes. No automatic modification of existing Batches 9–15.**
The run only writes under `output/parallel_t3/thesis_scalar_coordinate_batch16_v1/`.

## Installed full files (project root-relative)

- `Numerical/parallel_t3/complete/two_model_scalar_coordinate_batch16.py`
- `scripts/RunT3ScalarCoordinatesBatch16.py`
- `tests/test_t3_scalar_coordinates_batch16.py`
- `README_T3_SCALAR_COORDINATE_BATCH16.md`

These are **new, complete files**, not source diffs. Keep the existing
`two_model_master_wilson_batch13.py` and all Batch 15 files unchanged.

## Mathematics

Let `N_i` denote the exact archived Matchete coupling-direction tensors after
the verified Batch 15 source-to-spherical triplet rotation and `R_a` the
candidate RGBeta-potential fourth-derivative reference tensors. For each
field-content sector separately, solve

    N_i = sum_a R_a A[a,i]

by an overdetermined least-squares fit of the **entire** tensor, not just a
single component or previously guessed alias. The matrix `A` converts the
coefficient vector of **Matchete potential monomials** into the reference
potential-coordinate vector:

    g_reference = A @ lambda_native

An inverse is unique **only** if both bases span the same space with matching
independent dimensions. An underdetermined or missing archived native direction
is reported as incomplete and is NOT substituted with zero. Numerical
invariance of the reference quartic is tested by reconstructing complete
reference sectors with the source tensors.

`lambdaT3` is not passed off as one of the Hermitian baseline quartic
invariants. Its real and imaginary orientations are expanded directly from the
native source including `Bar[Coupling[...]]` and checked for real-linearity.
This is **internal source consistency**, not independent RGBeta phase evidence.

For the fixed-order EFT1 check, fully reconstructable sectors are expressed in
the Matchete native basis; incomplete sectors retain the original explicit
RGBeta reference tensor (rather than being guessed). The resulting quartic is
inserted into the **same candidate** seven-channel ADM used by Batch 13. The
calculation uses `ln(muS/muF)/(16 pi^2)` once and replaces the old pole-log
partition, exactly as Batch 13 did. The code checks that the result reproduces
saved Batch 13 to relative tolerance `1e-8`; an unexpected change is a blocker.

## PowerShell — installation and execution

Run from:

    C:\thesistemp\T3neutrinomassmodelpipeline

Download `t3_scalar_coordinate_batch16.zip`, then:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive `
  -LiteralPath "$HOME\Downloads\t3_scalar_coordinate_batch16.zip" `
  -DestinationPath . -Force

.\.venv\Scripts\python.exe -m pytest -q `
  tests\test_t3_scalar_coordinates_batch16.py

.\.venv\Scripts\python.exe -X utf8 `
  -m scripts.RunT3ScalarCoordinatesBatch16 --all --preflight

.\.venv\Scripts\python.exe -X utf8 `
  -m scripts.RunT3ScalarCoordinatesBatch16 --all --debug
```

To skip the final conditional SMEFT observable evaluation (while still doing
the EFT1 calculation), add `--skip-smeft` to the last command.

## Outputs

    output/parallel_t3/thesis_scalar_coordinate_batch16_v1/
      summary.json
      A_-4/result.json
      A_-4/recomputed_EFT1_delta.npz
      D_-2/result.json
      D_-2/recomputed_EFT1_delta.npz

Read `archived_native_vs_RGBeta_coordinate_sectors` in each `result.json` for
native/reference basis ranks, full coordinate transformation matrix, inverse
(where identifiable), directional residuals, missing directions and the
unresolved `lambdaT3` phase convention. The fields
`batch13_reproduction_relative`, `master_EFT1_log_change_relative`, and
`scalar_tensor_reconstruction_relative` quantify fixed-order equivalence.

## Limitations — important

1. These are transformations **relative to the candidate RGBeta reference
   potential definitions presently coded in the repository**. A source-level
   identification of every actual RGBeta generator invariant with this
   reference is still a distinct requirement.
2. The archived Matchete seed does not contain every independently evolved
   RGBeta scalar direction. Incomplete source coverage and the additional
   `D_-2` generalized quartics remain explicitly conditional; no invented
   native couplings are introduced.
3. Internal complex `lambdaT3` phase-linearity does not establish an
   independent absolute Matchete ↔ RGBeta phase/sign convention. Complex
   orientations of all D additional invariants are not independently proven.
4. The **GammaCC → Weyl absolute normalization**, full finite one-loop
   matching, and **independently proven gauge-complete EFT1 ADM** are not
   established by this batch. `physical_C5_certified` stays `false`.
5. No new Yukawa fit is performed. The neutrino observables remain
   conditional benchmark outputs; a pure basis change should not modify them.

## Tests and validation performed before packaging

The eight standalone linear-algebra tests passed locally. Python syntax was
compiled locally. **The full A/D repository integration has not been executed**
in this packaging environment (the locally archived NPZ/source datasets were
not supplied), and must be validated by running the above commands on the
Windows project checkout. A failed source or rerun consistency check is a
blocker, not a reason to force a coefficient mapping.

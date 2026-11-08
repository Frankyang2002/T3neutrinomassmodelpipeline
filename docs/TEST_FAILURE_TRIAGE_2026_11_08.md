# Test-suite triage after historical-test archival

Basis: user-provided Windows pytest log, 8 November 2026, working tree at
`f6baccfd9748b8b0f05d6f50f183663c0719643b` when the retirement plan was applied.
This document is a diagnostic inventory, not an assertion that every retained
unit or integration test is currently correct.

## Verified run results

- 110 test modules archived under `C:\thesistemp\T3_external_archives\test_cleanup_20261108T071306028000Z`.
- Core `scripts/RunCoreTests.ps1`: 189 passed.
- `scripts/RunCoreTests.ps1 -Full`: 4,453 passed, 73 failed, 51 skipped,
  1 collection error; 665 passing subtests (as pytest reported).
- The 14 maintenance tests in the v2 archive tool previously failed in three
  Windows-only test cases and passed in 11.

## Fixes in the accompanying package

These are **test-only** replacements; none changes calculations, matching,
reported numerical constants, the repository archive, or Git history.

1. `tests/python/test_reduce_historical_tests.py`: initialize temporary
   Git repositories with explicit LF blobs and convert newlines to CRLF exactly
   once during test simulation. This avoids a Windows-only test-fixture error.
2. `tests/test_intermediate_scalar_vector.py`: check for the actual
   `IntegratedField` rejection instead of requiring obsolete error-message
   wording. Failure still must raise `ValueError`.
3. `tests/test_threshold_boundary_naming.py`: require the two entry-point
   functions and scalar-only state constructors, not an obsolete comment.
   The boundary remains a state projection; no general Matchete matching is
   certified by this test.
4. `tests/test_t3_bc_native_y_llhh_kernel_v1.py`: normalize paths to `/`
   before checking project-relative `output/audits` locations on Windows.

The new maintenance tests were exercised in an isolated Linux test sandbox;
project-dependent tests must be run in the user's repository to confirm their
behavior on Windows.

## Separate unresolved groups (do not blanket delete)

### Missing saved research artifacts

- `tests/python/test_general_t3_two_threshold_numerical_c5.py`: approximately
  19 failing cases when the B/C two-threshold raw-LLHH inputs or audit JSON
  certificates are absent in `output/general_t3_mass_eigen_matchete/`.
- `tests/test_t3_bc_110101_mass_basis_map.py`,
  `tests/test_t3_bc_110101_potential_map.py`,
  `tests/test_t3_bc_live_110101_action_bridge.py`, and
  `tests/test_t3_bc_live_matchete_110101_registered_cg.py`: many failures share
  absent `output/general_t3_stage22/stage22_110101_wolfram_native_covariance.json`.

These are saved-output integration tests rather than standalone unit tests.
They should be moved to an opt-in integration category **after** testing their
actual file dependencies. Do not fabricate certificates or change missing-data
fail-closed behavior to pass them.

### Outdated test interfaces

- `tests/python/test_general_fermion_threshold_scan.py` still references
  `heavy_mass_basis_at_mu_f` and `mixed_boundary_fractions` even though the
  checked scan script now defines `_nonmixed_heavy_mass_basis` and
  `_mixed_fitted_scalar_fractions`.
- `tests/python/test_general_notation_index.py` assumes a legacy label index
  mapping and helper `_general_quartic_labels` absent from current adapter.
- `tests/python/test_pipeline_general_symbolic_dispatch.py` uses a stub module
  without the current `DEFAULT_DATA_ROOT` import contract.

These need source-based replacements or retirement as presentation-only legacy
contracts; changing production code to add obsolete APIs is not warranted.

### Collection issue

`tests/test_physical_b1_phase_transport.py` imports a function literally named
`test_flavour_covariance` from a production module. Pytest then collects the
import as a test requiring unavailable `y1`/`y2`/`mf` fixtures. The correct
source fix is an alias that does not begin `test_` **without** changing the
actual class method that calls the function. This is still pending.

### Numerically/physically significant failures

- `tests/python/test_general_t3_numerical_c5.py`: example result
  `137.71335734420055` versus expected `148.42297259166764`.
- `tests/python/test_general_t3_reference_c5.py`: example result
  `56.78838881868898` versus expected `61.2038969873076`.
- `tests/python/test_general_t3_quartic_action_identity.py` correctly
  fails closed on a changed Matchete/quartic normalization convention.
- `tests/python/test_general_eft1_integration.py` and
  `tests/python/test_general_two_stage_running.py`: mismatched EFT1 metadata.
- `tests/test_lo_couplings_integration.py`: local RGBeta model source differs
  from the audited historical source baseline; check the current model rather
  than overriding the source guard.

**Do not update expected numerical values until the underlying physics and
convention changes have been verified.**

## Recommended verification

```powershell
$repo = "C:\thesistemp\T3neutrinomassmodelpipeline"
Set-Location $repo
python -X utf8 -m pytest -q tests/python/test_reduce_historical_tests.py
python -X utf8 -m pytest -q tests/test_intermediate_scalar_vector.py tests/test_threshold_boundary_naming.py tests/test_t3_bc_native_y_llhh_kernel_v1.py
.\scripts\RunCoreTests.ps1
```

For full-suite diagnostics without the accidental imported helper (temporary
workaround, not a full test pass):

```powershell
python -X utf8 -m pytest -q --tb=line -k "not test_flavour_covariance" --maxfail=10
```

The full-suite outcome remains an informative integration diagnostic until
historical data-dependent tests are explicitly separated and core numerical
benchmarks are physically reconciled.

## EFT correctness condition

The supported physical order remains:

`UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT`.

No test reduction should claim full general matching or permit scalar-first
production solely to make a test pass.

# Windows baseline test repair — commit aa4dd2b7

## Scope

This is an **offline, local-only test update**. No production code, scientific
inputs, fitted Yukawa data, matching expressions, saved RGE runs, or figures are
modified. The replacement files target exactly repository commit
`aa4dd2b7b6d2c31494e45290a88004822a8c43c5`, verified read-only on GitHub.

The preceding baseline was one import error on full collection, followed by
`2120 passed, 21 failed, 9 skipped` when the obsolete C12 parsing test was
excluded. This number is a **pre-fix observation**, not a claim about the
post-fix suite, which still needs to be executed in the user's Windows checkout.

## Why each change is made

- `tests/test_c12_group_factor_pdf_tables.py`: replace tests for removed
  `_extract_c12_sections` and `transform_group_factor_tex` APIs with meaningful
  checks of the existing `Reports.C12GroupFactorPDF` table and section interface.
  The more extensive `tests/test_c12_group_factor_pdf.py` remains untouched.
- `tests/test_t3b_exploratory_coupled_benchmark.py`: test the **net** diagonal
  Weinberg coefficient `g5 - 3(Ye Ye^dag)ii`. The previous `g5 - (Ye Ye^dag)ii`
  did not include both flavour-leg corrections; this changes a test, not an RGE.
- `tests/python/test_general_t3_quartic_convention_c5.py`: explicitly require that the
  historical two-coupling source checker **reject** the post-Stage24 native
  RGBeta source rather than falsely calling the old map certified. Its
  independent synthetic source-contract tests are still retained.
- `tests/python/test_general_t3_tensor_aligned_geometry.py`: preserve exact
  50-slot/46-block geometry checks, now accepting the registered post-Stage24
  `CorrectedBCReusedInvariantCandidate` status in addition to historical phases.
  Continue to require `physical_coupling_map_certified == False`.
- `tests/python/test_general_t3_tensor_explicit_basis.py`: the JSON report
  omits zero-count status fields; use zero defaults and an exhaustive slot sum
  rather than indexing `ExplicitGeometryMismatch` when it is absent.
- `tests/python/test_general_t3_prefactors.py`: accept the documented
  post-Stage24 classification while retaining the source-normalisation,
  physical-map uncertainty, and signed two-dimensional map tests.
- `tests/python/test_general_t3_bc_covariance.py`: the historical
  *pre-four-invariant-fix* synthetic reconstruction only applies to that
  source phase. Skip just this phase-dependent test if saved source expressions
  are already corrected, while keeping independent covariance and fail-closed
  tests active. **This is not a substitute for a live post-correction audit.**
- `tests/python/test_t3_general_scalar_numerical_campaign.py`: patch both
  `sys.modules` and the `scripts` package attributes in synthetic test fixtures,
  because `from scripts import ...` can return the previously imported genuine
  module. The production same-checkout protection is kept unchanged.
- `scripts/RunBaselineTests.ps1`: enable Python UTF-8 mode for Windows.
  This addresses the observed `cp932` fixture I/O failures without monkeypatches
  or changes to production Unicode files.

## PowerShell installation

Run these commands from the project root. **Do not overwrite a different
local commit or edited versions of the eight replacement test files.**

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$expected = 'aa4dd2b7b6d2c31494e45290a88004822a8c43c5'
if ((git rev-parse HEAD).Trim() -ne $expected) { throw 'Different commit: inspect before installation' }
$targets = @(
    'tests/test_c12_group_factor_pdf_tables.py',
    'tests/test_t3b_exploratory_coupled_benchmark.py',
    'tests/python/test_general_t3_quartic_convention_c5.py',
    'tests/python/test_general_t3_tensor_aligned_geometry.py',
    'tests/python/test_general_t3_tensor_explicit_basis.py',
    'tests/python/test_general_t3_prefactors.py',
    'tests/python/test_general_t3_bc_covariance.py',
    'tests/python/test_t3_general_scalar_numerical_campaign.py'
)
if (git status --porcelain -- $targets) { throw 'Local test modifications detected: back up and inspect first' }
$backup = 'C:\thesistemp\T3_test_baseline_pre_fix_aa4dd2'
foreach ($path in $targets) {
    $dest = Join-Path $backup $path
    New-Item -ItemType Directory -Force -Path (Split-Path $dest) | Out-Null
    Copy-Item -LiteralPath $path -Destination $dest -Force
}
Expand-Archive -LiteralPath 'C:\path\to\T3_TestBaselineFix_aa4dd2.zip' -DestinationPath . -Force
.\scripts\RunBaselineTests.ps1
```

Update the ZIP path to the downloaded archive. If Windows PowerShell's execution
policy blocks running local scripts, instead run the test directly:

```powershell
$env:PYTHONUTF8 = '1'
python -m pytest -q --tb=line
```

## Acceptance criteria

- Pytest collects all tests without using `--ignore`.
- No encoding-related `cp932` failures remain.
- Historical symbolic diagnostics must not silently promote their uncertified
  native tensor/coupling maps to physically certified calculations.
- The T3-B exploratory regression passes with the correct two-leg coefficient.
- No previously passing production-pipeline test becomes a failure.
- Retain any new test failure output for a narrow correction; do not delete
  scientific validations just to obtain a green test suite.

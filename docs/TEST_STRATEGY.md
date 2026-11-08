# T3 testing strategy

This page is about tests, not a chronological progress report.

## What constitutes the supported physics workflow

1. Generate symmetry-allowed **general UV** interactions, including permitted Yukawa, scalar quadratic and quartic couplings. Do not infer complete numerical generality from a structural inventory.
2. Integrate out **all heavy fermions** to obtain an EFT with the heavy scalar modes still present.
3. Integrate out **all heavy scalars** (including mass-mixed eigenstates where applicable).
4. Match and evolve the SMEFT, retaining Weinberg-operator normalization and flavour structure.

Keep regression checks for group theory and charges, scalar mass matrices, Yukawa/mass-basis phases, fermion-first threshold ordering, complete heavy-field decoupling, EFT Wilson coefficients, RGE normalization, and core numerical/neutrino observables. An "offline/synthetic" test of a certificate is not equivalent to a live Matchete or RGBeta computation.

## Running tests

In PowerShell from repository root:

```powershell
.\scripts\RunCoreTests.ps1
.\scripts\RunCoreTests.ps1 -Full
```

The first command runs 18 source-based modules representing the minimal sanity selection; it **does not** certify the full physics. The `-Full` run covers all still-retained pytest tests and should be run after archival and before larger code changes. Neither command rewrites data, GitHub or the Lagrangian.

## Retirement policy

`python -X utf8 scripts/maintenance/ReduceHistoricalTests.py` writes a **read-only** plan at `output/cleanup_audit/historical_test_retirement_plan.json`.

The initial plan uses a fixed Git HEAD and conservative name patterns for presentation regressions, frozen/reduced historical campaigns, batch-labelled pilots, obsolete documentation checks, and four PDF tests that import `scripts/WriteT3GeneralSymbolicPDFs.py` (not tracked at the audited Git commit). It protects core tests and additional general-theory and covariance families.

A filename is not proof of irrelevance. The tool scans retained source/config files for literal references to selected tests; such files are *held*. Before archiving it compares every candidate with its tracked Git blob, then copies and SHA-256 verifies the entire archival set. Only then does it remove the originals. Archives go to a sibling `T3_external_archives` directory, not `git` or GitHub.

Apply after reviewing the plan and verifying the core tests:

```powershell
python -X utf8 scripts/maintenance/ReduceHistoricalTests.py --apply --confirm ARCHIVE_HISTORICAL_TESTS
```

To restore:

```powershell
python -X utf8 scripts/maintenance/ReduceHistoricalTests.py --restore "C:\thesistemp\T3_external_archives\test_cleanup_YYYYMMDDTHHMMSSffffffZ"
```

The archive is intentionally reversible. The remaining test suite may have pre-existing runtime failures, which must be reported rather than reclassified as passes or solved by changing physical conclusions. Do not delete `tests/reference`, `tests/fixtures` or `tests/wolfram` in this pass.

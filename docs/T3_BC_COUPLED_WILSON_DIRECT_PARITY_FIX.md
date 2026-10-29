# B(-1)/C(-1): direct-source 0.5 parity failure

## Failure

The coupled runner stopped before the Wilson ODE:

- `B_-1: master direct source disagrees with saved hierarchical beta: 0.5`
- `C_-1: master direct source disagrees with saved hierarchical beta: 0.5000000000000001`

The guard is deliberately retained at its original relative tolerance (`3e-4`). No outputs from these blocked runs are certified.

## Cause verified from source

`Numerical/experimental/T3BCCoupledWilsonTensor.py:flavor_pairs` previously built the cross-Yukawa flavour factor as

```
F12_old = [conj(Y1) diag(1/M) conj(Y2)^T
         + conj(Y2) diag(1/M) conj(Y1)^T] / 2
```

However the repository's existing `RGE/matching/WeinbergFlavorMatching.py:build_majorana_c5_flavor_matrix` and earlier `scripts/RunT3DirectWilsonSequential.py:pair_matrix` construct the full-flavour symmetric mixed contribution as an **unaveraged two-orientation sum**:

```
F12_prior = conj(Y1) diag(1/M) conj(Y2)^T
          + conj(Y2) diag(1/M) conj(Y1)^T
```

The new coupled runner therefore started with half the mixed `LLSS` tensor used in the saved hierarchical direct C5 baseline. The resulting direct-source beta was exactly half the reference value, producing the observed `0.5` mismatch. This discrepancy was not numerical integration error.

## Correction

Only the `F12` numerical flavour lift is changed to the unaveraged symmetric sum. The `F11` and `F22` contractions are unchanged. The seven-term Wilson master, quadrature, coupling maps, boundary hard C5, and strict direct/extra source parity checks are unchanged.

This aligns the **conditional candidate coupled calculation** with the previously established project flavour convention; it is **not an independent proof** that the historical full-flavour matching normalization is physically correct. The Stage17 extra-quartic normalization is still independently checked and must pass.

## Project files included

- `Numerical/experimental/T3BCCoupledWilsonTensor.py` — full replacement, preserves the parser and mass-symbol fixes
- `tests/test_t3_bc_coupled_wilson.py` — full replacement, adds exact ordered-flavour and mixed-tree regressions
- `docs/T3_BC_COUPLED_WILSON_DIRECT_PARITY_FIX.md` — this document

No production code outside the new coupled Wilson module is modified. No previous output is overwritten intentionally. GitHub is not modified.

## PowerShell

In `C:\thesistemp\T3neutrinomassmodelpipeline`, extract the ZIP, overwriting the two replacement files. Then:

```powershell
python -X utf8 -m pytest -q tests/test_t3_bc_coupled_parser_api.py tests/test_t3_bc_coupled_wilson.py
python -m scripts.RunT3BCCoupledWilsonSequential --model B_-1 --run-id continuous_v1 --resume
python -m scripts.RunT3BCCoupledWilsonSequential --all --run-id continuous_v1 --resume
```

If a prior **successful** run with different source hashes already exists, supply a fresh `--output-id bc_coupled_wilson_v2`. The strict hash guard should not be removed.

## Validation

In the isolated environment: 28 tests passed, 2 repository-dependent checks skipped; Python compilation passed. A complete B/C numerical calculation against the user's local benchmark data is **not performed here** and must be validated from the next user run. The expected immediate result is that the direct parity no longer reports `0.5`. The Stage17 extra parity and final Weinberg-subspace check may still reveal further issues; do not bypass them.

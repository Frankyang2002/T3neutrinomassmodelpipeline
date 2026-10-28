# Batch 12 — continuation of the exploratory B(-1)/C(-1) general-scalar fits

## Purpose

Continue the **already completed** Batch 11 Yukawa optimisations rather than repeat a fit from the historical restricted-potential seeds. The source fits provided by the user are:

| Case | Model | Previous run | Previous chi-squared | General scalar scenario |
|---|---|---|---:|---|
| `small_uv_frozen_5e-9` | `B_-1` | `b_small_fit_v1` | 118.57915264725226 | Seven additional quartics at 0.1\|lambdaT3\|, zero quadratic input |
| `large_uv_frozen_0p5` | `C_-1` | `c_large_fit_v1` | 116.30993286412242 | Same boundary rule, with lambdaT3=0.5 |

Both prior fits are marked `ExploratoryPoorFit`. The goal is chi-squared <= 5.86. The Batch 12 continuation is **not guaranteed to attain it**.

The historical fitted results, previous Batch 9–11 runs, Matchete symbolic coefficient, scalar parameters, and original fit target remain unchanged. All outputs are **new-only**.

## Actual numerical path

Every Yukawa trial uses the **existing Batch 11** `CandidateEvaluator` without altering its physics:

1. Set all 18 *real* Yukawa parameters on a fresh copy of the original historical configuration.
2. Integrate candidate UV RGBeta (when the original case enables UV running) and scalar-only EFT1 RGBeta with the same additional-quartic/quadratic boundary as the source Batch 11 run.
3. Determine heavy-fermion mass-basis Yukawas at `muF=1e5 GeV`, including Takagi/SVD when necessary.
4. Insert endpoint scalar mass mixing, diagonal scalar masses and two direct quartics into the saved **one-scale** Matchete coefficient with the Stage 23/26 **conditional** map `lambda_M = -lambda_RGBeta`.
5. Run SM+Weinberg from `muS=1e3 GeV` to `mu_low=100 GeV`; calculate the unchanged five-observable NuFIT chi-squared.

The optimiser also reuses the existing `Numerical.fitting.SobolBenchmarkSearch.refine_sobol_candidates` implementation, with **all 18 Yukawa coordinates** active simultaneously (instead of the previous staged 5/8/12). This is a numerical optimisation adjustment, not a new matching calculation.

### Important numerical change — local variable scaling

The previous best B(-1) Yukawas reach 1.54 while the previous best C(-1) Yukawas reach 3.09e-4. Using a common optimiser interval `[-3,3]` makes the dimensionless derivative steps excessively large for C(-1). Batch 12 uses individual finite local intervals, centred on **each saved best parameter**:

```
half_width_j = max(relative_window * abs(y_j),
                   floor_ratio * max_k(abs(y_k)), absolute_floor)
lower_j = max(-bound, y_j - half_width_j)
upper_j = min(+bound, y_j + half_width_j)
```

Defaults: `relative_window=0.5`, `floor_ratio=0.02`, `absolute_floor=1e-8`, `bound=3`. It computes the exact source seed once, then launches a single all-18 local least-squares stage with `max_nfev=60`. Optional additional local Sobol search is disabled by default (`--sobol-points 0`) to avoid unnecessary integrations. Set `--sobol-points 8` for a separate broader study only if needed.

## Install

Extract the ZIP at the project root:

`C:\thesistemp\T3neutrinomassmodelpipeline`

The ZIP has **only these three new files**:

- `scripts/ContinueT3GeneralScalarBCExploratoryFit.py`
- `tests/python/test_t3_general_scalar_bc_continuation.py`
- `docs/T3_BC_GENERAL_SCALAR_EXPLORATORY_REFIT_CONTINUATION_BATCH12.md`

The original Batch 11 runner remains untouched.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_BC_General_Scalar_Fit_Continuation_Batch12.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_general_scalar_bc_continuation.py

# First continue B(-1), small scalar, UV frozen:
python -m scripts.ContinueT3GeneralScalarBCExploratoryFit `
  --case small_uv_frozen_5e-9 `
  --model B_-1 `
  --previous-run-id b_small_fit_v1 `
  --run-id b_small_fit_v2 `
  --max-nfev 60

# Then continue C(-1), large scalar, UV frozen:
python -m scripts.ContinueT3GeneralScalarBCExploratoryFit `
  --case large_uv_frozen_0p5 `
  --model C_-1 `
  --previous-run-id c_large_fit_v1 `
  --run-id c_large_fit_v2 `
  --max-nfev 60
```

Optional broader local search if a continuation fails the fit screen: use **another** new run ID, with `--sobol-points 8`; do not overwrite `v2`.

Read the two summaries:

```powershell
Get-Content output\studies\general_scalar_extended\small_uv_frozen_5e-9\b_small_fit_v2\summary.md
Get-Content output\studies\general_scalar_extended\large_uv_frozen_0p5\c_large_fit_v2\summary.md
```

### New outputs

Each run writes four files only to its new study directory:

- `<B_m1/C_m1>_exploratory_fit.json` — full results, parameter vectors, optimisation trace, source hashes and diagnostic status;
- `<B_m1/C_m1>_exploratory_fitted_config.json` — complete distinct fitted configuration, provenance and 18 updated real Yukawas;
- `summary.json`;
- `summary.md`.

The command rejects an existing destination directory, a source run not marked as a completed fit, a mismatch between the saved best 18-vector and previous fitted config, a changed source physics-input hash, a different UV regime or scalar boundary, or a nonreproducible previous fit chi-squared. It does not edit the source file.

## Verification and interpretation

- Offline unit tests: 16 passed in the build environment; includes mocked full 18-parameter continuation, immutable source checks, prior-run consistency, adaptive scaling, and parsing of both user-supplied saved fit JSONs (the latter test auto-skips on Windows where that uploaded analysis ZIP is not present).
- Both actual source runs were checked directly from the user-uploaded ZIP: 18 parameters, expected statuses, and fitted-config values matched the saved best point exactly.
- **Not run here:** full RGBeta/Matchete/SM-Weinberg numerical continuation on the user's Windows repository. The local runs determine whether the fit screen is reached.
- All results remain **exploratory, not physically certified**. The Matchete/RGBeta common component convention, new-operator normalisations and full sequential mixed-scalar matching remain unresolved; the fitted neutrino observables are not independent predictions.

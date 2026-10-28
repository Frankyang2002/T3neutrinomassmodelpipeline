# Batch 17 — constrained B(-1)/C(-1) UV-running Sobol and local optimisation

## What is actually implemented

A new **conditional exploratory** fitter reuses all unchanged Batch 11/14/16
interfaces and saved, input-fingerprint-verified historical and exploratory
source configurations. Each trial updates **all 18 real Yukawa coordinates**
and **one additional real coordinate** `log10_lambdaT3` at 10^7 GeV. The
seven candidate additional scalar quartics start at `0.1*abs(lambdaT3)`;
`mu2_fraction=0` at the UV boundary and its RG-generated value is retained.

On every trial:

1. Run the actual saved candidate RGBeta UV state (10^7 -> 10^5 GeV).
2. Require positive mixed-scalar running-mass eigenvalues and
   `max(mScalar(muF))/min(MF(muF)) <= 0.5`.
3. For a valid UV hierarchy only, run the existing EFT1 candidate scalar RG
   state (10^5 -> 10^3 GeV) and require both scalar masses in `[0.5,5]` TeV.
4. Only if both mass gates pass, use the saved one-scale Matchete C5 with
   **conditional** Stage23/26 quartic mapping, heavy mass-basis rotation,
   and SM+Weinberg running to 100 GeV to compute the five NuFIT observables.
5. Recover the whitened five-component residual from the exact experimental
   target covariance; independently verify its squared norm against the
   returned C5 chi-square.

No scalar-first or sequential EFT matching is generated. These are mass
**running-parameter** inequalities, not a pole-mass decoupling certificate.

## Search method

- Inputs are the **two existing Yukawa textures** (historical and previous
  poor fit) with Batch15-successful common scales: B `0.005,0.01,0.015,0.02`;
  C `0.005,0.01,0.015`. Each is only a seed, not the final fit.
- Use the existing SciPy scrambled-Sobol generator from
  `Numerical.fitting.SobolBenchmarkSearch`, sampling locally around **all 18**
  Yukawa entries and a two-decade window for `log10(lambdaT3)` (clipped to
  global bounds). The `lambdaT3=5e-9/f^2` relation sets only the seed centre.
- Keep screened candidates, rank by **actual recomputed** chi-square, and
  call `scipy.optimize.least_squares` in dimensionless unit coordinates for
  **all 19 variables** from up to two best seeds.
- Ineligible candidates have **no assigned NuFIT chi-square**. The local
  solver instead sees an explicit nonzero penalty residual for a failed mass
  gate. The penalty is only an optimisation aid and is never reported as
  an oscillation statistic. Numerical exceptions likewise cannot become
  accepted candidates.
- Screened acceptance requires both mass conditions and `chi2 <= 5.86`.
  This procedure does **not** guarantee that a screened fit will be found.

## Project files

- `scripts/FitT3BCUVHierarchyConstrained.py` — **new, complete** runner.
- `tests/python/test_t3_bc_uv_hierarchy_constrained_fit.py` — unit and
  mocked full-run tests.
- `docs/T3_BC_UV_HIERARCHY_CONSTRAINED_FIT_BATCH17.md` — this document.

Other project modules, including Batch 16 and all production RG equations,
are unchanged.

## PowerShell commands

Run in `C:\thesistemp\T3neutrinomassmodelpipeline` after extracting the ZIP
into the project root:

```powershell
python -m pytest -q tests\python\test_t3_bc_uv_hierarchy_constrained_fit.py

# First pilot: one texture, selected mass-safe seed factors, modest budget.
python -m scripts.FitT3BCUVHierarchyConstrained `
  --model B_-1 --run-id b_uv_constrained_pilot_v1 `
  --seed-source historical_fit --seed-factors "0.005,0.01,0.02" `
  --sobol-points 8 --seed-count 1 --max-nfev 20

python -m scripts.FitT3BCUVHierarchyConstrained `
  --model C_-1 --run-id c_uv_constrained_pilot_v1 `
  --seed-source historical_fit --seed-factors "0.005,0.01,0.015" `
  --sobol-points 8 --seed-count 1 --max-nfev 20

Get-Content -Encoding UTF8 output\studies\general_scalar_extended\small_uv_running_5e-9\b_uv_constrained_pilot_v1\summary.md
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\small_uv_running_5e-9\c_uv_constrained_pilot_v1\summary.md
```

If the pilots remain poor, use **new run IDs** to widen the search. For
example, drop `--seed-source` and `--seed-factors`, use `--sobol-points 16`,
`--seed-count 2`, `--max-nfev 60` (significantly slower). Never reuse a run ID.

## Outputs and provenance

For each model, independent immutable directory:

`output/studies/general_scalar_extended/small_uv_running_5e-9/<run-id>/`

- `fit.json`: complete candidate search inventory counts, best screened fit,
  fitted Yukawa coordinates, UV topology coupling, predicted oscillation
  observables, running-mass diagnostics and source hashes.
- `fitted_config.json`: new complete input configuration **only if** a
  screened candidate was evaluated, including modified `lambdaT3` and
  Yukawas. Not a certified physical model configuration.
- `summary.md`: concise acceptance/limitations report.

Old campaigns, Batch 11–16 runs, historical fits and code are never overwritten.

## Scientific limitations

A successful Sobol+least-squares `chi2 <= 5.86` fits five observables using
19 free real coordinates, so it demonstrates candidate numerical fit
existence only, not predictivity. The physical Matchete–RGBeta map is not
certified. The existing one-scale coefficient does not include the complete
renormalised general-scalar sequential EFT matching. Positive running scalar
masses at endpoints do not prove vacuum stability or physical pole ordering.
This is a numerical feasibility benchmark, suitable for a thesis only with
explicit qualifications and caveats.

## CLI fix (Batch 17.1)

The first Batch 17 ZIP incorrectly referenced `a.scalar_min` and `a.scalar_max`
inside `main()` while argparse defines `--scalar-min-tev` and
`--scalar-max-tev`. This resulted in an `AttributeError` **before numerical
running or output creation**. The replacement CLI passes
`a.scalar_min_tev` and `a.scalar_max_tev` to `run_model`. Two new regression
tests execute the CLI dispatch with default and explicitly set scalar windows.
No physics functions, beta equations, matching, or fit parameters were changed.

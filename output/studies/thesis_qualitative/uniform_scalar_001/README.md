# Uniform 0.01 BSM scalar quartic benchmark — qualitative T3 refits

## Status and source

Standalone additional driver `scripts/RunT3UniformScalarBenchmark.py`, verified
against the current GitHub `main` versions of the accepted merged 16-model
source, general nonmixed evaluator, and B(-1)/C(-1) conditional scalar backend
on 2026-11-04. It imports the existing
`scripts/RunT3QualitativeScalarRefits.py` for accepted source selection, fit
optimisation, and output conventions. It **does not edit** any existing physics
code or input data, and does not overwrite the old scalar-refit study.

All included tests are isolated unit tests; **no complete numerical neutrino
fit has been executed with this package**. The user's full project checkout
and canonical RGE/Matchete exports are required for physical runs.

## UV boundary at 10^7 GeV

- **Set all allowed dimensionless BSM scalar quartics other than `lambdaT3`
  to real `+0.01`**, including `lambdaS1`, `lambdaS2`, `lambdaH1`,
  `lambdaH2`, `lambda12`, other ordinary representation-dependent quartics,
  and all additional full-general quartic operators in the UV beta layout.
- **Set any allowed `y1Conj`, `y2Conj` matrices to zero** at the UV boundary.
  Ordinary `y1`, `y2` are separately refitted as 18 real parameters.
- Scan `lambdaT3` = **0.0005**, **0.005**, **0.5**; it is *not* fixed to 0.01.
- Keep **SM** `lambdaH=0.25`, gauge couplings, SM Yukawas, fermion mass matrix,
  scalar mass-squared parameters and all threshold scales at saved accepted-fit
  values. In particular, do not overwrite SM `lambdaH` with 0.01.
- Retain the extra mixed **dimensionful scalar bilinear**, where allowed, at
  the saved fraction of `min(mS1Sq,mS2Sq)` (typically 0.01). It is not a
  dimensionless quartic and must not be set to 0.01 GeV^2 by analogy.
- The 16 accepted merged general-fit configurations provide initial `y1/y2`
  seeds and model metadata, **not already-accepted fits under the new
  quartic boundary**. Their original 5-observable fit values are not claimed
  for this new scalar sector.

B(-1)/C(-1) originally encode seven extra quartics as
`quartic_fraction * abs(lambdaT3)` with `quartic_fraction <= 1`. The new
standalone evaluator supplies **absolute** `+0.01` initial quartics in the
existing general scalar state directly. It uses the existing physical beta
payloads, UV/EFT1 integrators, C5, and low-energy fit. It does not clamp
quartics, reduce 0.01 for smaller `lambdaT3`, or modify existing source.

## What is computed

For each model and scanned `lambdaT3`, the driver starts from the accepted
original Yukawas scaled by `sqrt(lambda_original/lambda_new)` as a **seed
only**, and runs a bounded SciPy least-squares fit over 18 real Yukawa entries.
It records the 5 oscillation observable predictions, chi-square, fitted
Yukawa matrices and Frobenius norms, optional low neutrino masses, RGE and
matching diagnostics, original source paths and SHA256 hashes, and a qualitative
norm-vs-lambda slope. Default acceptance is `chi2 <= 5.86`. A small
regulariser selects a nearby branch because 18 parameters cannot be uniquely
constrained by five oscillation observables.

`--accept-seed` may avoid an expensive optimisation **only if the rescaled
seed already passes** the chi-square cutoff. Such points are explicitly
marked `AcceptedRescaledSeed`, **not** `AcceptedApproxFit`; omit this flag
if you require an independent optimiser call at every point.

The general matching is still a **conditional one-scale C5 construction**;
it is **not** a certified full sequential matching calculation
`UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT`.
For B/C the two direct quartics retain the backend's conditional signed
Matchete/RGBeta convention. Quartic equality is imposed in the existing
RGBeta UV coupling basis, not proven as a basis-independent physical relation.

## PowerShell commands (from the repository root)

Copy the contents of `t3_uniform_scalar_001_benchmark.zip` into the root of
the local repository:

```powershell
Expand-Archive -LiteralPath .\t3_uniform_scalar_001_benchmark.zip -DestinationPath . -Force

# Isolated boundary-condition and summary tests:
python -X utf8 -m unittest discover -s tests/python -p 'test_t3_uniform_scalar_benchmark.py' -v

# Read-only source and coupling inventory for A(-4), B(-1):
python -X utf8 -m scripts.RunT3UniformScalarBenchmark --models pilot --preflight

# Start with one model to control CPU cost:
python -X utf8 -m scripts.RunT3UniformScalarBenchmark `
  --models A_-4 --run-id uniform001_A_m4_v1 --max-nfev 5

# Mixed B(-1); unlike the earlier lambda-only fractional parameterisation,
# all three lambdaT3 points now have independent quartics = 0.01:
python -X utf8 -m scripts.RunT3UniformScalarBenchmark `
  --models B_-1 --run-id uniform001_B_m1_v1 --max-nfev 5

# Only if the pilot is useful; run all 16:
python -X utf8 -m scripts.RunT3UniformScalarBenchmark `
  --models all --run-id uniform001_all16_v1 --max-nfev 5
```

Optional fast exploratory mode: append `--accept-seed` to skip the optimiser
only for already-accepted rescaled seeds, with the explicit unoptimised status.
Do not represent those as independently refitted.

The output directory is:

`output/studies/thesis_qualitative/uniform_scalar_001/<run-id>/`

Each new run ID must be unique; overwriting an existing run directory is
refused. Summaries: `summary.json`, `comparison.csv`; full point outputs:
`<model-tag>/lambda_<value>/result.json`, `fitted_config.json`, and
`general_uv_overrides.json` for nonmixed models.

A failed point is written as `Failed`, with the exception recorded. `PoorFit`
means the numerical calculation returned but chi-square exceeded cutoff.
An exit code 2 signals at least one rejected/blocked model; inspect
`comparison.csv` rather than assuming all points were accepted.

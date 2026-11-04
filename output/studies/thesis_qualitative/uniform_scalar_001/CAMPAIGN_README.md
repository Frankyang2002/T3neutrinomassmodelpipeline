# T3 uniform-0.01 benchmark: new reference fits first, then lambdaT3 scan

## Purpose

These are **new independent Yukawa fits** under new UV scalar boundaries, NOT
reinterpretations of previously accepted fits. Existing accepted merged files are
used only as representation-specific `Y1,Y2` seeds and as provenance.

Benchmark at `muUV=1e7 GeV`:

- `lambdaS1`, `lambdaS2`, `lambdaH1`, `lambdaH2`, `lambda12`, and every other
  **symmetry-allowed BSM scalar quartic** except `lambdaT3`: **+0.01**, real,
  independent of `lambdaT3`.
- Allowed `y1Conj`, `y2Conj`: **0 at the UV boundary**. Their subsequent
  running is delegated to the actual RGEs; the statement is *not* that they
  are constrained to zero at every scale.
- `lambdaT3` stays at the *old accepted fit's value* for a new reference fit:
  `0.1` for A/B/C/D, `0.05` for E.
- Optional **subsequent** scan: `lambdaT3=0.0005, 0.005, 0.5`, keeping
  all other BSM quartics fixed at +0.01. Only **newly accepted** references
  provide the initial `Y1,Y2` seeds.
- SM `lambdaH=0.25`, gauge inputs, SM Yukawas, scalar/fermion masses,
  thresholds, and dimensionful scalar-mixing-bilinear boundary policy are
  **unchanged**. Dimensionful mass-squared bilinears are not set to 0.01.
- The fitting variables are all 18 **real** components of `Y1,Y2`; imaginary
  components are fixed to zero. Five neutrino oscillation observables are
  fitted to the existing experiment target; regularisation selects one local
  branch because these observables do not uniquely determine 18 parameters.
- Accepted fit means `chi2<=5.86` within the existing pipeline, NOT proof of
  complete phenomenological viability.

### Limits

`NonMixedEvaluator` is used for the 14 non-mixed representations; B_-1/C_-1
use the registry-based B/C full-general running and a local absolute-quartic
initial state. No change to shared physics files. The supported Matchete C5
is still **conditional one-scale**, including the B/C sign map: the desired
strictly sequential `UV -> integrate all fermions -> integrate all scalars
-> SMEFT` is not certified by this code.

## Installation (PowerShell, from repository root)

Download `t3_uniform001_reference_campaign.zip` into the repository root.
Extract the archive:

```powershell
Expand-Archive -LiteralPath .\t3_uniform001_reference_campaign.zip -DestinationPath . -Force
```

Archive project-relative files:

- `scripts/RunT3UniformScalarBenchmark.py`: complete reusable scalar-boundary evaluator from the previous uniform-0.01 package.
- `scripts/RunT3UniformScalarCampaign.py`: new two-stage campaign controller.
- `tests/python/test_t3_uniform_scalar_benchmark.py`: complete scalar-boundary tests.
- `tests/python/test_t3_uniform_scalar_campaign.py`: new regression tests.
- `output/studies/thesis_qualitative/uniform_scalar_001/CAMPAIGN_README.md`: this file.

The repository's existing
`scripts/RunT3QualitativeScalarRefits.py` is required and is not modified.

## Test and read-only preflight

```powershell
python -X utf8 -m unittest discover -s tests/python -p "test_t3_uniform_scalar_campaign.py" -v
python -X utf8 -m unittest discover -s tests/python -p "test_t3_uniform_scalar_benchmark.py" -v
python -X utf8 -m scripts.RunT3UniformScalarCampaign --stage preflight --models pilot
```

The preflight reads accepted merged outputs and validates scalar boundary
inventories. It does not run RGEs, Matchete matching, or fits, and creates no
outputs.

## Stage 1: fit new reference Yukawas

Start with the cheap-to-debug A(-4) case:

```powershell
python -X utf8 -m scripts.RunT3UniformScalarCampaign `
  --stage reference `
  --models A_-4 `
  --run-id uniform001_A_m4_reference_v1 `
  --max-nfev 5
```

Then test the mixed B(-1) backend separately:

```powershell
python -X utf8 -m scripts.RunT3UniformScalarCampaign `
  --stage reference `
  --models B_-1 `
  --run-id uniform001_B_m1_reference_v1 `
  --max-nfev 5
```

Read the new reference status:

```powershell
Import-Csv .\output\studies\thesis_qualitative\uniform_scalar_001\campaigns\uniform001_A_m4_reference_v1\reference_comparison.csv |
  Format-Table model,lambdaT3,status,chi2 -AutoSize
```

If `PoorFit` or `Failed`, rerun with `--max-nfev 10` and a **new run ID**.
Do not reuse/overwrite earlier results. `--max-nfev 5` is the *solver's*
function-call budget, not a bound on all physical integrations: 18 real
parameters require many finite-difference Jacobian evaluations. No physical
fit is claimed until actual `chi2` output exists.

For all 16 independently:

```powershell
python -X utf8 -m scripts.RunT3UniformScalarCampaign `
  --stage reference `
  --models all `
  --run-id uniform001_all16_reference_v1 `
  --max-nfev 5
```

The optional `--accept-seed` switch checks the seed's actual physics `chi2`
and **skips** optimisation if the seed is already accepted. Such a result is
explicitly `AcceptedRescaledSeed`, not an independently optimised fit. The
switch is OFF by default to meet the request for new optimisations.

## Stage 2: scan only newly accepted reference models

Run with **the same run ID and same model selection**. Every model must have
its new accepted reference saved by stage 1; models without one are blocked.

```powershell
python -X utf8 -m scripts.RunT3UniformScalarCampaign `
  --stage scan `
  --models A_-4 `
  --run-id uniform001_A_m4_reference_v1 `
  --max-nfev 5
```

Combined future all-16 study:

```powershell
python -X utf8 -m scripts.RunT3UniformScalarCampaign `
  --stage scan `
  --models all `
  --run-id uniform001_all16_reference_v1 `
  --max-nfev 5
```

There is also `--stage both`, to do a new reference and then scan in one
invocation, but **reference first** is recommended to see whether the chosen
new quartics permit acceptable Yukawa fits before spending time on 48 scans.

## Outputs and integrity

Output root is
`output/studies/thesis_qualitative/uniform_scalar_001/campaigns/<run-id>/`:

- `reference/<model-tag>/result.json`: original source hashes, five
  observables, `chi2`, fitted matrices, norms and RGE/matching diagnostics.
- `reference/<model-tag>/fitted_config.json`: complete new benchmark UV state
  and **actually fitted** Yukawa matrices (if fit computed).
- `reference/<model-tag>/general_uv_overrides.json`: complete full-general
  UV overrides for non-mixed models (if fit computed).
- `reference_comparison.csv`, `reference_summary.json`: quick 16-model tally.
- `scan/<model-tag>/lambda_<value>/...`: per-scan-point fits and UV states.
- `scan_comparison.csv`, `scan_summary.json`: three-point summaries.

The campaign rejects scan references unless the status meets the threshold,
physical reference lambda is unchanged, UV scalar/SM/representation inputs
match the new uniform policy, the two fitted Yukawa matrices agree across
saved files, and all original accepted-source files still have the recorded
SHA256 hashes. Existing original/other-study outputs are never overwritten.

**Validation in the packaging environment:** 22 isolated unit tests passed
(9 baseline fit logic, 6 scalar boundary, 7 new campaign). The full repository
physics dependencies and generated beta/matching payloads were not present;
**no physical reference fit or new scan was executed** there.

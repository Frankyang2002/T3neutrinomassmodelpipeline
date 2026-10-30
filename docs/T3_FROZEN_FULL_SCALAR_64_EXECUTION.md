# Frozen full-general-scalar 16 × 4: executable candidate pipeline

**Target project:** `Frankyang2002/T3neutrinomassmodelpipeline` (local working copy ONLY).  
**Source inspected:** GitHub `main`, `029ed78f54fac577987afebeed64a7b74cdbba35` (30 October 2026).  
**Installation:** extract the package ZIP at the project root, preserving paths. This package contains only **new** file paths. No original source, numerical data, PDF generator, or original Lambda HTML generator is replaced.

## New project files

- `scripts/T3Frozen64Sources.py` — strict 16-model/4-benchmark catalogue, canonical paths, SHA-256 source checks, manifest validation.
- `scripts/RunT3FrozenFullScalar64.py` — new frozen-UV boundary on the existing full-scalar EFT1 engine, plus optional existing conditional neutrino continuation.
- `scripts/GenerateT3FrozenFullScalarGamma.py` — genuine model-/benchmark-specific T1/T2/T3 perturbation response (no fit).
- `Numerical/plotting/BuildT3FullScalarFrozen64HTML.py` — one HTML report generator reusing original RGE and oval HTML canvases.
- `scripts/RetireSupersededT3HTMLGenerators.ps1` — strictly guarded local archival of three old **HTML-only** entry points after complete reports exist.
- `tests/test_t3_frozen_full_scalar64.py` — manifest, integrity, selection, and empty-report tests.
- `docs/T3_FROZEN_FULL_SCALAR_64_EXECUTION.md` — these instructions.

## Exact physics selection

| Report key | Fitted-Yukawa **seed** | Historical `lambdaT3` case | UV RG | Full scalar EFT1 |
|---|---|---|---|---|
| `optimised_low` | each model's own saved high/low historical fitted Yukawas | `small_uv_frozen_5e-9` | **frozen** | yes, candidate |
| `optimised_high` | each model's own saved historical fitted Yukawas | `large_uv_frozen_0p5` | **frozen** | yes, candidate |
| `shared_a_m2_low` | the A(-2) same-case saved Yukawas inserted into each model | `small_uv_frozen_5e-9` | **frozen** | yes, candidate |
| `shared_a_m2_high` | the A(-2) same-case saved Yukawas inserted into each model | `large_uv_frozen_0p5` | **frozen** | yes, candidate |

For class E, the source historical configurations use half the nominal topology coupling. The numerical code calls `RunT3GeneralScalar48Campaign.historical(...)` **only to verify and read the original scenario-specific Yukawa/config inputs**; it does **not** call the old bulk48 evaluator for the running or neutrino calculation.

**Explicit NEW illustrative boundary assumption:** the existing Stage29 additional quartic values (both its new and pre-existing RGBeta coordinates) are imposed at the frozen **100-TeV** boundary; a nonzero allowed holomorphic quadratic is set to `0.01 * min(mS1²,mS2²)` by default, real phase zero. This is a choice to define one general-scalar numerical example; it is NOT a new independently optimised general-scalar fit, and the original Stage29 values were defined for a different UV boundary. Do not treat this as the uniquely physically preferred general-scalar boundary.

The existing `RunT3FullScalarContinuous.py` verifies the full candidate EFT1 beta inventory, evolves all original and allowed additional scalar quartics and quadratic terms where present, and constructs the scalar mass spectrum. The new driver directly reuses those functions. It writes a flat **UV-boundary guide** (explicitly `uv_ode_evaluations=0`) and integrates only the EFT1 interval from 100 TeV to 1 TeV. It also saves actual *ordinary* quartic samples, distinguished from additional couplings in the HTML selectors.

The existing `RunT3FullScalarNeutrino.py` handles low-energy C5/SM-Weinberg continuation. Its scientific limitations remain in force: six self-bilinear A/D/E C5 controls omit the nonzero mixing correction, B/C mixed results are conditional one-scale, and full sequential Wilson matching / Matchete–RGBeta normalization are not certified. A numerically successful candidate is **not** a fully validated general-scalar neutrino fit. No general-scalar Yukawa re-optimisation is implemented in this package.

## Run from the local project root (PowerShell)

```powershell
# Set this to your actual checkout directory.
cd C:\thesistemp\T3neutrinomassmodelpipeline

# Confirm your local code first; stop to inspect if it diverges from 029ed78f54fa.
git rev-parse --short HEAD
git status --short

# After extracting the ZIP, run the contract tests.
python -X utf8 -m pytest -q tests\test_t3_frozen_full_scalar64.py

# READ ONLY: resolve exact historical fitted configs and Stage29/beta inputs
# for all 64 combinations. Does not create run outputs.
python -X utf8 -m scripts.RunT3FrozenFullScalar64 --plan

# Optional single-case smoke test: new separate run ID, no historical overwrite.
python -u -X utf8 -m scripts.RunT3FrozenFullScalar64 `
    --run-id frozen_general64_smoke_v1 `
    --benchmark optimised_low --models A_+0

# Full 64-case general-scalar candidate EFT1 + conditional neutrino run.
python -u -X utf8 -m scripts.RunT3FrozenFullScalar64 `
    --run-id frozen_general64_v1

# Resume the same exact run/settings after an interruption.
python -u -X utf8 -m scripts.RunT3FrozenFullScalar64 `
    --run-id frozen_general64_v1 --resume

# Genuine T1/T2/T3 gamma scans from the NEW frozen-general-scalar results.
python -u -X utf8 -m scripts.GenerateT3FrozenFullScalarGamma `
    --run-id frozen_general64_v1 --scan-id gamma_scans_v1

# Resume scans only if source hashes, gamma grid, and generator are unchanged.
python -u -X utf8 -m scripts.GenerateT3FrozenFullScalarGamma `
    --run-id frozen_general64_v1 --scan-id gamma_scans_v1 --resume

# Presentation-only: original Lambda RGE/oval canvas templates.
python -X utf8 -m Numerical.plotting.BuildT3FullScalarFrozen64HTML `
    --run-id frozen_general64_v1 --scan-id gamma_scans_v1

Start-Process .\Reports\output\display_general_scalar_frozen_64\frozen_general64_v1\index.html
```

If you need initial scalar plots without neutrino matching, run `RunT3FrozenFullScalar64 --skip-neutrino` with **another run ID**. If any source hash or parameter choice changes, select a new `--run-id`; do not silently resume incompatible input states. If `--plan` fails, fix source provenance before running the ODE campaign. Individual failures are displayed as `BLOCKED` manifest cells rather than silently substituted with older reduced/quartic data.

## Saved data and output locations

```text
output/studies/general_scalar_extended/full_scalar_frozen_64/
  frozen_general64_v1/
    manifest.json                 # 64 strict unique cells, source SHA-256
    optimised_low/
      A_p0/
        continuous_trajectory.json
        A_p0_eft1_quadratic_running.json
      neutrino_v1/
        A_p0_neutrino.json
    optimised_high/...
    shared_a_m2_low/...
    shared_a_m2_high/...
    gamma_scans_v1/
      optimised_low/A_p0.csv      # 3 directions x 11 gamma values
      optimised_low/A_p0.json     # exact source hashes + coverage
      ...
Reports/output/display_general_scalar_frozen_64/
  frozen_general64_v1/
    index.html                    # 16 x 4 status/coverage table
    coverage.csv
    provenance.json
    optimised_low/interactive_rge.html
    optimised_low/sensitivity_ovals.html
    optimised_high/...
    shared_a_m2_low/...
    shared_a_m2_high/...
```

All HTML is standalone and uses the existing Lambda canvas template text. The RGE view **does not** invent C5(mu) or PMNS curves. The sensitivity ovals are only displayed when all 3 directions, the shared gamma=1 centre and the same five neutrino observables are numerically available and reproduce the corresponding saved central case. The historical NuFIT/JUNO contours are *Gaussian display approximations* sourced from `docs/FROZEN_LAMBDA_TARGET_REFERENCE.json`, not official joint likelihood contours.

## Legacy HTML-only generator retirement

**Do not archive the original** `InteractiveModelComparison.py`, `BuildFittedModelEnvelopes.py`, `GeneralScalarDisplayLambdaBridge.py`, `BuildGeneralScalarLambdaReports.py`, or the historical `BuildLambda5e9Displays.py`. The new report **depends on** `BuildGeneralScalarLambdaReports.py`. Keep all PDF generators.

After obtaining **64/64** verified RGE datasets and **64/64** complete full-scalar gamma scans, preview the local archival policy:

```powershell
.\scripts\RetireSupersededT3HTMLGenerators.ps1
```

To move (not delete) `BuildFrozenGeneralScalarLambdaReports.py`, `BuildT3FrozenLambdaNative.py`, and `BuildT3FrozenDisplaySummary.py` into `.cleanup_quarantine`, **only** if coverage and cross-import checks pass:

```powershell
.\scripts\RetireSupersededT3HTMLGenerators.ps1 -RunId frozen_general64_v1 -Apply
```

This archival step is deliberately guarded; do **not** run `scripts/CleanReportOutput.ps1 -Apply` for this task.

## What has and has not been tested here

- **Tested in the package environment:** all Python files compile, and 7 isolated contract tests pass, including a 64-slot HTML skeleton with zero fabricated results.
- **Inspected in the live GitHub source:** the required existing full-scalar numerical methods, input schemas, historical frozen benchmark definitions, reporting templates, and the current saved 16-model candidate continuous run.
- **Not run:** project-local 64-case input preflight, actual SciPy EFT1 integration, C5 evaluation, full gamma scans, and the final HTML against live numerical outputs. The full local repository and physics dependencies are not mounted in this environment. The new numerical integration and report coupling still need an execution test on your machine; on failure, retain the `BLOCKED` output and inspect the exact traceback before proceeding.

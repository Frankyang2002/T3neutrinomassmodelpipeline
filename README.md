# T3 radiative neutrino-mass model pipeline

This repository constructs ordinary T3 models, exports matching and RGE inputs, evolves renormalisable couplings across UV/EFT regimes, and evaluates the Weinberg coefficient and neutrino observables.

## Start here

- **Architecture, research progression, and status:** [`docs/PROJECT_GUIDE.md`](docs/PROJECT_GUIDE.md)
- **Physics and conventions:** [`docs/PHYSICS.md`](docs/PHYSICS.md)
- **Numerical implementation:** [`docs/NUMERICAL.md`](docs/NUMERICAL.md)
- **Saved results and their scope:** [`docs/RESULTS.md`](docs/RESULTS.md)
- **Historical research record:** [`docs/RESEARCH_HISTORY_INDEX.md`](docs/RESEARCH_HISTORY_INDEX.md)
- **Regeneration notes:** [`docs/REGENERATE_CURRENT_REPORTS.md`](docs/REGENERATE_CURRENT_REPORTS.md) (historical saved-data report rebuild; not a verified fresh symbolic rebuild)

## Software entry points

- `pipeline.py` — orchestration of model construction, RGEs, matching and optional numerical stages. Inspect `python pipeline.py --help` for study options.
- `Numerical/fitting/ScanCLI.py` — execute ordinary split-scalar numerical studies from JSON configs.
- `Numerical/orchestration/PipelineNumericalResults.py` — numerical diagnostics and reports inside the production pipeline.
- `scripts/maintenance/SmokeTestFullPipeline.py` — one-point operational test using existing generated artifacts.

### Verified operational smoke test

One existing ordinary T3-B (`alpha=-1`) configuration executed successfully through UV running, intermediate EFT1, matched `C5`, SM + Weinberg running, and oscillation-observable evaluation, returning `PASS_END_TO_END_OPERATIONAL`. This establishes software operation for **one benchmark using existing inputs**; it does not certify the physical matching conventions, all representations, or regenerability of those inputs from scratch.

From a Python environment with project dependencies installed:

```powershell
$env:PYTHONUTF8 = '1'
python scripts/maintenance/SmokeTestFullPipeline.py
```

**If legacy outputs were archived, this command may fail because its saved Matchete/RGBeta JSON inputs are absent.** Regenerate required artifacts first or restore them from the external archive; absence of inputs is distinct from a broken numerical integrator.

## Directory roles

| Directory | Purpose |
| --- | --- |
| `Lagrangian/` | Model, interactions, symbolic construction and matching |
| `RGE/` | RGBeta inputs, symbolic RGE generation and intermediate EFT definitions |
| `Numerical/core/`, `Numerical/running/` | Numerical state and coupled scale evolution |
| `Numerical/fitting/`, `Numerical/orchestration/` | Benchmarks, scans, orchestration and output |
| `physics/` | Mass matrices and neutrino observables |
| `configs/` | Input configurations and benchmark definitions |
| `tests/` | Regression tests |
| `docs/` | Current documentation and consolidated research history |
| `output/`, `Reports/output/` | Generated artifacts; historical contents may have been archived outside the repo |

## Status and scope

The production ordinary split-scalar numerical path is operational for the tested configuration. The separate general-coupling and exploratory paths may require unverified RGBeta–Matchete basis conversions or different matching approximations. Do not interpret their conditional outputs as fully validated physical predictions. See the project guide for details.

## Repository housekeeping

The cleanup and archival operations are **local-only** and should be performed with a preview, an external archive, and a manifest. Never move a file used by source generation, tests or configs solely because its name contains `stage`, `probe` or `experimental`.

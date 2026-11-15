# T3 radiative neutrino mass models

Research tools for one-loop T3 radiative neutrino-mass models. The project covers 16 electroweak field assignments, with symbolic model construction, renormalisation-group evolution, effective-theory matching, numerical benchmarks and neutrino observables.

## Physical calculation

The intended effective field theory (EFT) sequence is

```text
UV theory: Standard Model + heavy fermions + heavy scalars
    | integrate out all heavy fermions
    v
Scalar EFT: Standard Model + remaining heavy scalars + induced operators
    | integrate out all heavy scalars
    v
SMEFT: Standard Model + relevant Wilson coefficients
    | evolve to the low scale
    v
Neutrino mass matrix -> mass eigenvalues, mass splittings and PMNS parameters
```

This defines the physical ordering. Some matching contributions, general interactions and numerical interfaces are not connected into one verified full-general calculation; see [Capabilities and limitations](docs/IMPLEMENTATION_STATUS.md).

## Getting started

Run from the repository root in Windows PowerShell. Python 3.12 is the documented baseline.

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt

python -X utf8 pipeline.py --help
python -X utf8 full16_symbolic.py --help
.\scripts\RunCoreTests.ps1
```

Use an existing functioning environment instead of recreating it. New symbolic calculations require a local Wolfram Language installation with Matchete and RGBeta; `pip` does not install these tools. See [Installation](docs/INSTALLATION.md).

## Calculations

| Entry point | Purpose |
| --- | --- |
| `pipeline.py` | Combined model, matching, running and optional numerical commands. The available physics depends on the selected options. |
| `full16_symbolic.py` | General-interaction symbolic investigation of the 16 ordinary T3 models; **not** a complete all-model numerical neutrino-mass fit. |
| `Numerical/fitting/ScanCLI.py` | Numerical parameter scans using compatible saved coefficient inputs and an explicit configuration. |
| `scripts/RunCoreTests.ps1` | Fast regression checks; `-Full` requests the broader tests. |

For a saved numerical scan configuration:

```powershell
python -X utf8 -m Numerical.fitting.ScanCLI path\to\existing_scan.json --dry-run
```

Replace the example with an actual configuration path.

## Repository guide

| Location | Scientific purpose |
| --- | --- |
| `common/`, `model/` | T3 fields, representations and allowed models |
| `Lagrangian/` | Interaction Lagrangians and symbolic matching |
| `RGE/` | Group theory factors and running equations |
| `Numerical/`, `physics/` | Numerical evolution, fitting, neutrino masses and mixing |
| `configs/`, `input/` | Calculation inputs and benchmark definitions |
| `Reports/`, `display/` | Results and figures |
| `tests/`, `validation/` | Regression and scientific checks |
| `docs/` | Physics methods, computational conventions and reproducibility |

Begin with the [documentation index](docs/README.md), [physical model](docs/MODEL_AND_PIPELINE.md), [numerical method](docs/NUMERICAL.md), and [results](docs/RESULTS.md).

The code supports the November 2026 Honours thesis *Calculating Quantum Corrections in Generalised Scotogenic Radiative Neutrino Mass Models*. The 16-model catalogue and its scientific scope are described in [Thesis scope](docs/THESIS_SCOPE.md).

## Interpretation of results

Numerical fits, symbolic outputs and reported comparisons are valid within their recorded inputs, interaction bases and approximations. A fitted oscillation observable is not an independent prediction. Preserve the input configurations, coefficient sources and saved manifests for reproducibility; some files in `output/` and `Reports/output/` are inputs to later calculations.

See [Capabilities and limitations](docs/IMPLEMENTATION_STATUS.md) and [Reproducibility](docs/REPRODUCIBILITY.md) before treating a calculation as complete one-loop matching.

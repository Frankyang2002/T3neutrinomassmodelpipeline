# T3 radiative neutrino mass models

This repository provides symbolic and numerical tools for studying one-loop T3 radiative neutrino mass models. It covers the ordinary catalogue of **16 electroweak representation and hypercharge assignments**, with model construction, renormalisation-group evolution, effective-field-theory (EFT) matching, neutrino observables and benchmark comparisons.

The calculation is organised by **physical content**, not by the order in which the software was developed:

```text
UV theory: Standard Model + heavy fermions F + heavy scalars S1, S2
    | integrate out all heavy fermions
    v
Scalar EFT: Standard Model + S1, S2 + induced effective operators
    | evolve coefficients; integrate out all heavy scalars
    v
SMEFT: Standard Model + Weinberg operator
    | evolve to the low scale
    v
Neutrino mass matrix -> masses, mass-squared differences and PMNS parameters
```

The physical EFT order is **fermions first, scalars second, SMEFT last**. Distinct mass thresholds within a particle sector require matching at the appropriate scales. This is the intended calculation framework; not every general-interaction and matching option is connected end-to-end in the current code.

## Start here

**Windows PowerShell**, from the repository root:

```powershell
# Set up a Python environment (Python 3.12 is the documented baseline).
py -3.12 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt

# Discover available calculation modes without starting an expensive run.
python -X utf8 pipeline.py --help
python -X utf8 full16_symbolic.py --help

# Run the core Python regression checks.
.\scripts\RunCoreTests.ps1
```

An existing compatible environment need not be recreated. Fresh symbolic generation and matching additionally require a local Wolfram Language installation with **Matchete** and **RGBeta**; these are not installed by `pip`. Details are in [Installation](docs/INSTALLATION.md).

## What each entry point does

| Entry point | Purpose |
| --- | --- |
| `pipeline.py` | Model selection, symbolic matching, RGE and optional numerical workflows; available modes and limitations depend on the selected options. |
| `full16_symbolic.py` | General-interaction symbolic study for the ordinary 16-model catalogue. This is **not** a complete general numerical neutrino-mass fit. |
| `Numerical/fitting/ScanCLI.py` | Numerical scans using compatible, already-generated source coefficients and an explicit JSON configuration. |
| `scripts/RunCoreTests.ps1` | Fast regression selection; `-Full` requests the broader suite. |

For the numerical scan, inspect the supplied configuration and its saved source files before running:

```powershell
python -X utf8 -m Numerical.fitting.ScanCLI path\to\existing_scan.json --dry-run
```

Replace the example path with an existing configuration; do not treat it as a bundled filename.

## Repository guide

| Location | Responsibility |
| --- | --- |
| `common/`, `model/` | Particle representations, model catalogue and threshold definitions |
| `Lagrangian/` | Wolfram model building and matching |
| `RGE/` | Group invariants, beta functions and effective-operator evolution |
| `Numerical/` | Numerical states, RG integration, fitting and analysis |
| `physics/` | Neutrino mass matrices and mixing observables |
| `configs/`, `input/` | Benchmark definitions and calculation inputs |
| `Reports/`, `display/` | Report generation and presentation outputs |
| `tests/`, `validation/` | Regression, algebraic and physics-specific checks |
| `scripts/`, `tools/` | Research drivers and supporting utilities |
| `docs/` | Scientific methods, software interfaces and reproducibility |

Read the [documentation index](docs/README.md) for the physics, code architecture, numerical methods and results.

## Connection to the Honours thesis

This codebase supports the November 2026 Honours thesis *Calculating Quantum Corrections in Generalised Scotogenic Radiative Neutrino Mass Models*. The reader-facing scope is the **16-model T3 catalogue**, the fermion-first/scalar-second EFT prescription, model-dependent mass and mixing, stage-running comparisons, and compensated coupling-sensitivity studies. Start with [Thesis scientific scope](docs/THESIS_SCOPE.md) and the [thesis figure/data map](docs/THESIS_ARTIFACT_MAP.md).

An **external-archive-first** local curation utility is supplied for removing historical development notes and LaTeX build byproducts without touching scientific code or thesis data:

```powershell
python -X utf8 scripts/maintenance/CurateThesisRepository.py plan
```

Its plan identifies protected science, review-only artifacts and safe historical-document candidates. It does not automatically remove numerical exports or source code.

## Interpretation and limitations

The repository contains **working components** for general T3 model construction, symbolic matching and numerical RG studies, alongside restricted or conditional benchmark pathways. The available source does **not** demonstrate a single verified, complete-one-loop, fully general-interaction, fermion-first-to-SMEFT numerical fitting workflow for all 16 models. A small fit residual or a passing software test does not establish that every matching contribution is included.

Always associate numerical results with their input configuration, source exports, approximations and run manifest. Saved files under `output/` and `Reports/output/` can be inputs to subsequent analyses; they should not be removed solely because they look generated.

See [Physics and EFT conventions](docs/MODEL_AND_PIPELINE.md), [capabilities and limitations](docs/IMPLEMENTATION_STATUS.md), and [reproducing results](docs/REPRODUCIBILITY.md).

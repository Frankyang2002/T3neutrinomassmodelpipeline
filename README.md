# T3 radiative neutrino-mass models

<!-- T3_CURRENT_SCOPE_BEGIN: 2026-11-10 -->

**Document scope:** Project overview; current A(-4)/D(-2) task is not finished.

**Current implementation status (10 November 2026):** The requested A(-4)/D(-2) **full-general, complete-one-loop, fermion-first -> joint-scalar -> SMEFT Yukawa fit is not yet an executable end-to-end numerical workflow**. The repository already implements general Matchete/RGBeta components and an actual symbolic `RunSequentialT3Matching` stage engine, but its general Wilson matching outputs are not fully connected to the two-model numerical fit. The current thesis Sobol pilot is **restricted/conditional** (analytic fermion triangle and selected scalar bubble), even if it reaches a low chi-square. See [implementation status](docs/IMPLEMENTATION_STATUS.md), [A(-4)/D(-2) feasibility](docs/T3_A4_D2_GENERAL_PIPELINE_FEASIBILITY_2026_11_10.md) and [implementation batches](docs/T3_A4_D2_FULL_GENERAL_IMPLEMENTATION_BATCHES.md). Earlier numerical results and verification claims below refer to their **dated, historical configurations** unless an individual run manifest establishes otherwise.

<!-- T3_CURRENT_SCOPE_END -->
A symbolic and numerical research pipeline for one-loop Majorana neutrino-mass generation in the T3 topology. The project constructs Standard Model extensions with heavy fermions and scalar multiplets, evaluates renormalisation-group evolution and effective-field-theory (EFT) thresholds, and studies the resulting neutrino masses and mixing.

The ordinary model catalogue covers **16 representation and hypercharge assignments**. The general-interaction inventory and the restricted numerical benchmarks are distinct: the presence of a symmetry-allowed operator in the symbolic model does **not** imply that its physical matching and numerical evolution have been fully validated.

## Calculation overview

```text
UV theory: Standard Model + heavy fermions + heavy scalars
    |
    | integrate out all heavy fermions
    v
Intermediate EFT: Standard Model + heavy scalars + induced operators
    |
    | evolve couplings and Wilson coefficients
    | integrate out all heavy scalars
    v
SMEFT: Standard Model fields + Weinberg operator
    |
    | low-energy running
    v
Neutrino mass matrix, masses and mixing observables
```

The fermion-first ordering is fundamental to this project's target calculation. Scalar-first matching is not a supported production route. A common-threshold calculation is a separate comparison, not equivalent to the sequential fermion-first result.

## Start here

| Topic | Reference |
| --- | --- |
| **Installation, dependencies and first run** | [Installation](docs/INSTALLATION.md) |
| **Implemented, verified and incomplete calculations** | [Implementation status](docs/IMPLEMENTATION_STATUS.md) |
| Model fields, symmetries and EFT sequence | [Model and calculation contract](docs/MODEL_AND_PIPELINE.md) |
| Gauge representations, interactions and conventions | [Physics](docs/PHYSICS.md) |
| Scalar potential and controlled approximations | [Scalar sector](docs/SCALAR_SECTOR.md) |
| Evolution, fitting and numerical methods | [Numerical methods](docs/NUMERICAL.md) |
| Benchmarks, results and their limitations | [Results](docs/RESULTS.md) |
| Documentation by subject | [Documentation index](docs/README.md) |
| Source tree and software entry points | [Repository structure](docs/REPOSITORY_STRUCTURE.md) |

## Running the software

Create a Python environment and install the dependencies as described in [Installation](docs/INSTALLATION.md). From the repository root in PowerShell:

```powershell
python -m pip install -r requirements.txt
python -X utf8 pipeline.py --help
python -X utf8 full16_symbolic.py --help
.\scripts\RunCoreTests.ps1
```

`pipeline.py` coordinates the study workflows. `full16_symbolic.py` is the dedicated general-interaction **symbolic** entry point; it is not a certified general-theory numerical matching pipeline. The core suite is a fast software regression check, not a physical completeness certificate.

The symbolic matching and beta-function workflows require their appropriate Wolfram Language, Matchete and RGBeta environment, as well as the relevant model-dependent inputs. Saved generated data may be required for numerical examples and report regeneration. Consult [Installation](docs/INSTALLATION.md), [Implementation status](docs/IMPLEMENTATION_STATUS.md) and [Numerical methods](docs/NUMERICAL.md) before interpreting a completed run as a fresh computation.

## Mathematical conventions

Electric charge is defined by \(Q=T_3+Y\), with

\[
Y(S_1)=\frac{\alpha}{2},\qquad
Y(S_2)=\frac{\alpha+2}{2},\qquad
Y(F)=\frac{\alpha+1}{2}.
\]

The numerical convention for the dimension-five Weinberg coefficient is

\[
m_\nu=-\frac{v^2}{2}\,C_5.
\]

The charge, flavour, scalar mass-basis and Wilson-coefficient conventions are documented in the physics and calculation references above.

## Scope and verification

The repository includes operational numerical benchmarks and source-based regression tests. However, the complete symmetry-allowed T3 interaction inventory is **not yet a certified, end-to-end general-coupling calculation**: invariant normalisations, mass-basis conversions, one-loop matching and some sequential scalar thresholds still require independent validation.

Results obtained with a restricted scalar boundary or a conditional matching prescription must retain those labels. They should not be presented as independent predictions of the fully general theory. The [results reference](docs/RESULTS.md) separates completed runs, fitted inputs, exploratory calculations and verified conclusions.

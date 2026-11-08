# Repository structure and terminology

This page identifies the **role of each software component**. The names of existing scientific modules and generated data are deliberately preserved because imports, saved source manifests and reproducibility checks may depend on their exact paths.

## Main directories

| Path | Role |
| --- | --- |
| `common/`, `model/` | Model definitions, field identities, threshold plans and study selection |
| `Lagrangian/` | Wolfram Language model construction, symmetry-allowed interactions and matching orchestration |
| `RGE/group_factors/` | SU(2) group algebra, invariant tensors and normalisation checks |
| `RGE/running/` | UV and intermediate-EFT beta-function data and Wilson-coefficient running |
| `RGE/matching/` | Matching and Weinberg-coefficient interfaces |
| `Numerical/core/` | Numerical states, coupling vectors and beta-function evaluation |
| `Numerical/running/` | Scale evolution and heavy-field threshold interfaces |
| `Numerical/fitting/` | Parameter scans, objectives and numerical fits |
| `Numerical/orchestration/` | Combined numerical workflows and result assembly |
| `physics/` | Neutrino mass matrices, diagonalisation and observables |
| `configs/`, `input/` | Study configurations and required inputs |
| `Reports/` | Analytical reports and visualisations |
| `tests/`, `validation/` | Unit, integration, physics and source-validation checks |
| `docs/` | Current technical references and preserved research record |
| `scripts/` | Study-specific commands, export utilities and investigation drivers |
| `output/`, `Reports/output/` | Generated calculations and saved study artifacts; some are required by numerical reruns |

## Entrypoints

- **Main pipeline:** `pipeline.py` — coordinates model selection, matching and running stages. Inspect its command-line options rather than assuming all combinations are supported.
- **General symbolic study:** `full16_symbolic.py` — convenience entry point for the 16-model general-interaction symbolic route. It does not imply that its outputs have certified general numerical matching.
- **Numerical fitting:** `Numerical/fitting/ScanCLI.py` and the workflows under `Numerical/orchestration/`.
- **Fast regression suite:** `scripts/RunCoreTests.ps1`.
- **Full retained test suite:** `scripts/RunCoreTests.ps1 -Full`; some extended tests depend on generated inputs or experimental interfaces and are not necessarily green on a checkout without their artifacts.

## EFT terminology

| Scientific stage | Meaning |
| --- | --- |
| **UV theory** | All Standard Model and T3 heavy fields are dynamical |
| **Fermion threshold** | The heavy fermions are integrated out; heavy scalars remain |
| **Intermediate scalar EFT** | Scalar-containing EFT and fermion-induced operators, with no heavy T3 fermions |
| **Scalar threshold** | All remaining T3 heavy scalars are integrated out |
| **SMEFT / Weinberg running** | Standard Model fields and the dimension-five coefficient \(C_5\) |
| **Low-energy observables** | Masses, squared-mass differences and leptonic mixing parameters |

The older names `EFT1`, `stage 1` and `stage 2` can appear in existing source files and saved outputs. They are historical interfaces, not alternative physical definitions. Follow **field content** and the integration order rather than inferring the EFT from a stage number alone.

## Interpretation of results

- **Symbolic operator inventory:** lists candidate allowed interactions; it does not certify all Clebsch-Gordan normalisations or threshold calculations.
- **Restricted benchmark:** a defined choice of couplings/boundary conditions within the model; it is not a different symmetry principle.
- **Numerical fit:** parameters have been adjusted to reproduce targets; a low objective value is not an independent neutrino prediction.
- **Conditional result:** depends on a stated approximation or matching convention not yet fully validated.
- **Core test pass:** confirms selected software contracts only; it does not certify the complete EFT matching.

## Reproducibility

Do not relocate active Python/Wolfram modules or rename their imported symbols solely to improve their appearance. Likewise, generated filenames may be referenced literally by saved JSON manifests or verification scripts. Older research tools are intentionally kept separate from the principal reader-facing commands, rather than being silently reclassified as production physics.

For further details, see [Model and calculation contract](MODEL_AND_PIPELINE.md), [Numerical methods](NUMERICAL.md), [Testing strategy](TEST_STRATEGY.md) and [Results](RESULTS.md).

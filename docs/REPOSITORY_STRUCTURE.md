# Repository structure

The software is organised by the scientific responsibility of each module. Existing source paths are retained for compatibility with Python imports, Wolfram evaluations, exported coefficients and saved run manifests.

| Path | Responsibility |
| --- | --- |
| `pipeline.py` | Overall pipeline command-line entry point |
| `full16_symbolic.py` | General symbolic study of the ordinary 16-model catalogue |
| `common/`, `model/` | Field definitions, model catalogue and threshold planning |
| `Lagrangian/` | Symmetry-allowed interactions, Wolfram model construction and matching |
| `RGE/group_factors/` | Electroweak group tensors, invariants and checks |
| `RGE/running/` | RGEs for the UV and intermediate effective theory |
| `RGE/matching/` | Matching interfaces and Weinberg coefficients |
| `Numerical/core/` | Numerical states and beta-function evaluation |
| `Numerical/running/` | Differential equation integration and threshold handoffs |
| `Numerical/fitting/` | Parameter optimisation and scans |
| `Numerical/orchestration/` | Coordination of numerical calculations |
| `physics/` | Neutrino mass matrices, diagonalisation and oscillation observables |
| `configs/`, `input/` | Configurations and input coefficients |
| `Reports/`, `display/` | Report generators and presentation material |
| `tests/`, `validation/` | Software, symbolic and physics validation |
| `scripts/` | Reproducible run drivers and specialised calculations |
| `tools/` | Auditing, figure rebuilds and repository utilities |
| `docs/` | Scientific method, interfaces, limitations and results |
| `output/`, `Reports/output/` | Calculated artifacts, stored coefficients and output manifests |

## Working with generated outputs

Many saved results are used as **inputs** by subsequent calculations. Do not move or delete `output/` or `Reports/output/` merely to reduce the apparent size of the repository. Reproduction requires a specific benchmark and its matching source files, not simply a copied figure.

## Scientific terminology

- **UV theory:** Standard Model plus the heavy T3 scalar and fermion fields.
- **Scalar EFT:** the theory remaining after integrating out all heavy fermions; heavy scalars and generated effective operators remain.
- **SMEFT:** the theory after integrating out all heavy scalar multiplets.
- **Weinberg running:** scale evolution of the lepton-number-violating dimension-five coefficient, used to construct the low-energy Majorana mass matrix.

These names describe actual **field content**. Internal filenames or manifest keys should not be used to infer physical validity or completeness.

See [Software architecture](PROJECT_GUIDE.md), [EFT contract](MODEL_AND_PIPELINE.md) and [Capabilities](IMPLEMENTATION_STATUS.md).

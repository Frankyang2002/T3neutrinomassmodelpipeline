# Finding source and data

The repository groups calculations by scientific responsibility. Most readers should begin with the [main README](../README.md), [calculation guide](PROJECT_GUIDE.md) and [results overview](RESULTS_OVERVIEW.md), rather than navigating directly through implementation files.

| Location | What it contains |
| --- | --- |
| `common/`, `model/` | T3 field representations, model definitions and threshold conventions |
| `Lagrangian/` | Symmetry-allowed interactions, model construction and symbolic matching |
| `RGE/` | Group-theory factors, renormalisation-group equations and effective-operator evolution |
| `Numerical/` | Numerical evolution, parameter optimisation and analysis |
| `physics/` | Neutrino mass matrices and oscillation observables |
| `configs/`, `input/` | Benchmark and calculation inputs |
| `Reports/`, `display/` | Numerical tables, figures and presentation material |
| `output/`, `Reports/output/` | Saved coefficients, computed results and reproducibility information |
| `tests/`, `validation/` | Software and scientific checks |
| `docs/` | Model description, methods, assumptions, results and reproducibility |
| `scripts/`, `tools/` | Supporting commands and specialist calculations |

The principal commands are `pipeline.py` for available combined calculations, `full16_symbolic.py` for the 16-model symbolic study, and `Numerical/fitting/ScanCLI.py` for configured numerical scans. Use the documented command-line help to determine which operations and source inputs are supported.

## What should be preserved

Saved outputs can be inputs to subsequent studies. In particular, files holding matched coefficients, generated RGEs, fitted Yukawa values, numerical configurations or figure provenance should not be discarded or renamed simply to simplify the directory listing. Their original paths may be required by import statements or archived run manifests.

Technical filenames do not establish the physics or accuracy of a result. For the scientific definition of the effective theories, see [the T3 convention guide](MODEL_AND_PIPELINE.md), and for the evidence supporting each calculation, see [Capabilities and limitations](IMPLEMENTATION_STATUS.md).

# Research and workflow commands

The main public entry points are in the repository root:

- `pipeline.py`: model selection and the coordinated matching/running workflow.
- `full16_symbolic.py`: general-interaction symbolic model generation and reporting.
- `scripts/RunCoreTests.ps1`: a fast set of project-critical regression tests; `-Full` runs the retained extended suite.

Other scripts provide **model-specific studies**, symbolic export, independent validation, numerical experiments and figure generation. A script being present here does not mean it is part of the fully verified production path.

## Where the substantive implementations live

| Calculation | Directory |
| --- | --- |
| T3 model and invariant interaction construction | `../Lagrangian/` |
| Group theory and renormalisation-group data | `../RGE/` |
| Numerical evolution, fitting and orchestration | `../Numerical/` |
| Neutrino observables | `../physics/` |
| Technical and physical assumptions | `../docs/` |

Source identifiers referring to prior study stages are preserved where a direct rename would break imports, checkpoints, saved provenance or Mathematica commands. Consult the [repository structure](../docs/REPOSITORY_STRUCTURE.md) and [documentation index](../docs/README.md) before running a specialised script.

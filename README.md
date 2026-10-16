# T3 radiative neutrino-mass pipeline

An honours-physics project implementing symbolic construction, matching and renormalisation-group evolution of T3 radiative neutrino-mass models, and numerical threshold running, benchmark studies, fits and diagnostics.

This README is the **project entry point**, not a statement that every phenomenological model is validated. The authoritative source code and saved numerical inputs, rather than old report figures, determine calculations.

## Where to start

| Task | Location |
|---|---|
| Select a model and dispatch a study | `pipeline.py`, `common/`, `model/`, `studies/` |
| Build UV interactions and matching | `Lagrangian/` |
| Group factors, RGEs, threshold matching and final Weinberg coefficient | `RGE/` |
| Numerically evolve, fit, extract and plot | `Numerical/` |
| Low-energy neutrino masses and oscillation observables | `physics/` |
| Generate symbolic/presentation reports | `Reports/` |
| Regression tests and validation | `tests/`, `validation/` |
| Numerical and benchmark configuration | `configs/`, `Numerical/configs/`, `Numerical/benchmarks/` |

## Documentation

- [`docs/PHYSICS.md`](docs/PHYSICS.md): model conventions, Lagrangian, matching, one-loop power counting, group factors and thresholds.
- [`docs/NUMERICAL.md`](docs/NUMERICAL.md): implementation ownership, benchmark evolution, fitting and testing.
- [`docs/RESULTS.md`](docs/RESULTS.md): benchmark provenance, verified numerical findings, limitations and references.

The previous `INFO_*.md`, `README_*.txt`, `docs/archive/` and `Reports/notes/` materials are **not deleted by this package**. They remain source material for checking whether the consolidation has retained everything important.

## PowerShell commands

From the repository root, in the project's working Python environment:

```powershell
python pipeline.py --help
python -m pytest -q

# Regenerate figures and interactive comparisons from saved numerical diagnostics
python -m Numerical.plotting.BuildDisplayResults --output Reports/output/display

# Inspect available fitted sensitivity studies before running them
python -m Numerical.orchestration.RunFittedSensitivityComparison --help
```

The saved-results display builder is presentation-only; it does not run Matchete, fitting or the RGEs. Full `python pipeline.py --full` workflows may invoke external Wolfram/Matchete/RGBeta programs and be costly. Do not run a full symbolic pipeline merely to rebuild figures.

### Outputs

`output/full/` holds raw saved numerical diagnostics. Some numerical study records (including scan JSON/CSV) live under `Reports/output/full/`. Generated PNG/PDF/HTML can be rebuilt from the appropriate saved data when the relevant generator is available; they are not automatically physics sources of truth.

## Interpretation cautions

- The nominal 16-model × four-scenario grid has 64 numerical cases. A successful numerical trajectory is **not** evidence of agreement with oscillation data.
- Fitted reference Yukawa matrices are deliberately calibrated to NuFIT. Agreement of fitted observables is not an independent model prediction.
- Four fixed `small/large` labels describe benchmark inputs. Later feasibility fits need not respect strict regime boundaries unless they are explicitly imposed.
- Standard production calculations use a fermion-first threshold hierarchy; see the detailed documents before changing threshold order or conventions.
- Keep numerical source data, generated matching artifacts and benchmark manifests when cleaning reports.

## Source

Repository: https://github.com/Frankyang2002/T3neutrinomassmodelpipeline

This consolidated documentation was drafted from existing repository notes; **no executable source file is changed by the ZIP**.

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

### Current scientific reference

- [`docs/PHYSICS.md`](docs/PHYSICS.md): model conventions, Lagrangian, matching, one-loop power counting, group factors and thresholds.
- [`docs/SCALAR_SECTOR.md`](docs/SCALAR_SECTOR.md): **authoritative scalar-potential justification** — full operator inventory, leading-one-loop selection, zero-boundary power counting, quadratic mixing, numerical scope, and limitations.
- [`docs/NUMERICAL.md`](docs/NUMERICAL.md): implementation ownership, benchmark evolution, fitting and testing.
- [`docs/RESULTS.md`](docs/RESULTS.md): benchmark provenance, verified numerical findings, limitations and references.
- [`docs/T3_THESIS_PLAN.md`](docs/T3_THESIS_PLAN.md): working thesis narrative and figure plan; individual historical statements must be checked against completed results.

### Historical research and implementation records

Start with [`docs/RESEARCH_HISTORY_INDEX.md`](docs/RESEARCH_HISTORY_INDEX.md) to find an original research note by filename, locate its consolidated section, and check the recorded source hash. Historical records document experiments, corrections and provisional findings; they do **not** supersede the current scientific reference or certify unresolved calculations.

| Topic | Consolidated record |
|---|---|
| SU(2) tensors, quartics and RGBeta conventions | [`RESEARCH_GENERAL_T3_TENSORS_AND_QUARTICS.md`](docs/RESEARCH_GENERAL_T3_TENSORS_AND_QUARTICS.md) |
| UV running and matching | [`RESEARCH_GENERAL_T3_RUNNING_AND_MATCHING.md`](docs/RESEARCH_GENERAL_T3_RUNNING_AND_MATCHING.md) |
| Matchete and mixed scalar mass eigenstates | [`RESEARCH_MATCHETE_AND_MASS_EIGEN.md`](docs/RESEARCH_MATCHETE_AND_MASS_EIGEN.md) |
| EFT1 thresholds and Weinberg operator | [`RESEARCH_EFT1_THRESHOLDS_AND_WEINBERG.md`](docs/RESEARCH_EFT1_THRESHOLDS_AND_WEINBERG.md) |
| Numerical fits, feasibility and benchmarks | [`RESEARCH_NUMERICAL_FITS_AND_BENCHMARKS.md`](docs/RESEARCH_NUMERICAL_FITS_AND_BENCHMARKS.md) |
| General-scalar studies and pilots | [`RESEARCH_GENERAL_SCALAR_AND_PILOTS.md`](docs/RESEARCH_GENERAL_SCALAR_AND_PILOTS.md) |
| Repository cleanup and migration history | [`REPOSITORY_MAINTENANCE_HISTORY.md`](docs/REPOSITORY_MAINTENANCE_HISTORY.md) |

The historical documents have been consolidated by topic. The external archival manifests and Git history preserve original paths and provenance; saved machine-readable scientific audits and numerical results remain separate from these explanatory records. For report regeneration, see [`docs/REGENERATE_CURRENT_REPORTS.md`](docs/REGENERATE_CURRENT_REPORTS.md).

## PowerShell commands

From the repository root, in the project's working Python environment:

```powershell
python pipeline.py --help
python -X utf8 -m pytest -q

# Regenerate figures and interactive comparisons from saved numerical diagnostics
python -m Numerical.plotting.BuildDisplayResults --output Reports/output/display

# Inspect available fitted sensitivity studies before running them
python -m Numerical.orchestration.RunFittedSensitivityComparison --help
```

The saved-results display builder is presentation-only; it does not run Matchete, fitting or the RGEs. Full `python pipeline.py --full` workflows may invoke external Wolfram/Matchete/RGBeta programs and be costly. Do not run a full symbolic pipeline merely to rebuild figures.

### Outputs

`output/full/` holds raw saved numerical diagnostics. Some numerical study records (including scan JSON/CSV) live under `Reports/output/full/`. Generated PNG/PDF/HTML can be rebuilt from the appropriate saved data when the relevant generator is available; they are not automatically physics sources of truth.

The final reduced-scalar numerical report is under `output/studies/reduced_scalar_final_reports/final_reduced_v2/`, with fit inputs in `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/`. The general-scalar exploratory campaigns and conditional zero-boundary pilots are distinct datasets; consult `docs/SCALAR_SECTOR.md` before interpreting or comparing them.

## Interpretation cautions

- The nominal 16-model × four-scenario grid has 64 numerical cases. A successful numerical trajectory is **not** evidence of agreement with oscillation data.
- Fitted reference Yukawa matrices are deliberately calibrated to NuFIT. Agreement of fitted observables is not an independent model prediction.
- Four fixed `small/large` labels describe benchmark inputs. Later feasibility fits need not respect strict regime boundaries unless they are explicitly imposed.
- Standard production calculations use a fermion-first threshold hierarchy; see the detailed documents before changing threshold order or conventions.
- The reduced scalar potential is an explicit **zero-additional-coupling boundary benchmark**, not the unrestricted general gauge-invariant potential or a certified full general-scalar EFT matching result.
- Keep numerical source data, generated matching artifacts and benchmark manifests when cleaning reports.

## Source

Repository: https://github.com/Frankyang2002/T3neutrinomassmodelpipeline

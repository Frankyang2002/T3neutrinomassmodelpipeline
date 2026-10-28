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
- [`docs/SCALAR_SECTOR.md`](docs/SCALAR_SECTOR.md): **authoritative scalar-potential justification** — full operator inventory, leading-one-loop selection, zero-boundary power counting, quadratic mixing, numerical scope, and limitations.
- [`docs/NUMERICAL.md`](docs/NUMERICAL.md): implementation ownership, benchmark evolution, fitting and testing.
- [`docs/RESULTS.md`](docs/RESULTS.md): benchmark provenance, verified numerical findings, limitations and references.
- [`docs/T3_THESIS_PLAN.md`](docs/T3_THESIS_PLAN.md): working thesis narrative and figure plan; individual historical statements must be checked against completed results.

The earlier scalar Batch 1–8 documents and historical stage notes are **retained for now**. `docs/SCALAR_SECTOR.md` consolidates their scalar-relevance argument but does not replace the machine-readable scientific audits. Obsolete documentation will be archived only after reference and dependency verification; this documentation-only change does not delete it.

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

This documentation consolidation is based on a read-only inspection of commit `aa4dd2b7b6d2c31494e45290a88004822a8c43c5`. No executable source file or saved scientific data is changed by the documentation ZIP.

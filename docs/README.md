# Documentation

The documentation is organised here by **scientific purpose**, rather than by the order in which individual calculations were developed. Consult the current overview first; older study notes are retained for reproducibility and should not be interpreted as independently verified results.

## Physics and calculation

- [Model and calculation contract](MODEL_AND_PIPELINE.md) — symmetries, fields, interactions, charge conventions and the target EFT sequence.
- [Physics](PHYSICS.md) — the T3 topology, SU(2) representations, matching conventions and physical scope.
- [Scalar sector](SCALAR_SECTOR.md) — the general operator inventory, mass mixing and controlled restricted-potential benchmarks.
- [General symbolic physics audit](GENERAL_SYMBOLIC_PHYSICS_AUDIT.md) — what has and has not been established for the additional interactions.

## Numerical analysis and results

- [Numerical methods](NUMERICAL.md) — state vectors, RG evolution, threshold interfaces, fits and numerical conventions.
- [Results](RESULTS.md) — existing comparison data, benchmark interpretations and limitations.
- [Regenerating current reports](REGENERATE_CURRENT_REPORTS.md) — rebuilding reports from their saved inputs, distinct from regenerating all matching expressions.

## Software reference

- [Repository structure](REPOSITORY_STRUCTURE.md) — directories, principal entry points and where calculations are implemented.
- [Testing strategy](TEST_STRATEGY.md) — core and extended tests; the difference between unit, synthetic, generated-artifact and live symbolic checks.
- [Project guide](PROJECT_GUIDE.md) — system architecture and a technical record of prior work.

## Historical research record

The subject-organised volumes retain the derivations and investigations behind the current implementation:

- [General tensors and scalar quartics](RESEARCH_GENERAL_T3_TENSORS_AND_QUARTICS.md)
- [General running and matching](RESEARCH_GENERAL_T3_RUNNING_AND_MATCHING.md)
- [Matchete conventions and scalar mass eigenstates](RESEARCH_MATCHETE_AND_MASS_EIGEN.md)
- [Fermion and scalar thresholds, Wilson coefficients](RESEARCH_EFT1_THRESHOLDS_AND_WEINBERG.md)
- [Numerical scans and benchmarks](RESEARCH_NUMERICAL_FITS_AND_BENCHMARKS.md)
- [General scalar exploratory studies](RESEARCH_GENERAL_SCALAR_AND_PILOTS.md)
- [Historical research document index](RESEARCH_HISTORY_INDEX.md)

These records include preliminary analyses and superseded conventions. They are evidence of research progression, **not** a replacement for the current physics and numerical contracts. File identifiers used in scripts, native exports or saved manifests are kept unchanged so that results remain traceable.

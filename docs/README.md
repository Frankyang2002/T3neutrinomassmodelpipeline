# Documentation

The documentation is arranged by scientific question: model definition, calculation, numerical evaluation and results. Start with the physics and then consult implementation details as needed.

## Theory

- [T3 neutrino-mass physics](PHYSICS.md): model motivation, one-loop mechanism and representations.
- [Model and EFT calculation](MODEL_AND_PIPELINE.md): interactions, conventions and the fermion-first, scalar-second matching prescription.
- [Scalar sector](SCALAR_SECTOR.md): scalar operators, permitted mixing and mass bases.

### Technical reference

- [UV models and running](UV_MODEL_AND_RUNNING.md): field content, interactions and symbolic-to-numerical conventions.
- [Fermion matching](FERMION_MATCHING.md): integrating out all heavy fermions and interpreting one-loop matching.
- [Scalar EFT and SMEFT](SCALAR_EFT_AND_WEINBERG.md): scalar-field evolution, scalar matching and Weinberg running.
- [Numerical verification](NUMERICAL_VALIDATION.md): inputs, fitting, observable extraction and tests.

## Computing observables

- [Numerical methodology](NUMERICAL.md): running, matching interfaces, fits and extraction of neutrino observables.
- [Results](RESULTS.md): reported model comparisons and their stated assumptions.
- [Thesis scope](THESIS_SCOPE.md) and [figure/data provenance](THESIS_ARTIFACT_MAP.md): the 16-model study and preserved results.

## Running and assessing calculations

- [Installation](INSTALLATION.md): Python and symbolic-software requirements.
- [Program structure](PROJECT_GUIDE.md) and [directory reference](REPOSITORY_STRUCTURE.md): where scientific components are implemented.
- [Capabilities and limitations](IMPLEMENTATION_STATUS.md): what is supported and what remains conditional or incomplete.
- [Verification](TEST_STRATEGY.md): how checks relate to physical correctness.
- [Reproducibility](REPRODUCIBILITY.md): configurations, source coefficients and run manifests.

**Scientific sequence:** UV theory → integrate out all heavy fermions → scalar EFT → integrate out all heavy scalars → SMEFT → low-energy observables.

Saved numerical studies and specialised research utilities can retain historical identifiers for compatibility. Those identifiers do not define additional physical theories or establish completeness.

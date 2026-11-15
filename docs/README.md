# Documentation

This documentation describes the **scientific model**, the **software that evaluates it**, and the **evidence needed to interpret a numerical result**. The principal interfaces are `pipeline.py` for the combined calculation and `full16_symbolic.py` for the 16-model symbolic comparison.

## Thesis scope and figure provenance

- [Thesis scope](THESIS_SCOPE.md): chapters, sixteen models, benchmark assumptions and EFT conventions.
- [Thesis artifact map](THESIS_ARTIFACT_MAP.md): preserved figure/table sources, fitted Yukawas and known provenance ambiguities.

## Scientific model and calculation

- [Physics and T3 representations](PHYSICS.md): one-loop mass generation, group structure and electroweak representations.
- [Model and EFT contract](MODEL_AND_PIPELINE.md): symmetry-allowed interactions, conventions and the required **UV → scalar EFT → SMEFT** matching sequence.
- [Scalar sector](SCALAR_SECTOR.md): scalar operators, allowed mixing and mass bases.
- [Numerical method](NUMERICAL.md): running, matching interfaces, fitting and observable extraction.

## Using the software

- [Installation](INSTALLATION.md): Python and external Wolfram/Matchete/RGBeta requirements.
- [Software architecture](PROJECT_GUIDE.md): entry points, module responsibilities and data flow.
- [Directory reference](REPOSITORY_STRUCTURE.md): locations of source code and inputs.
- [Capabilities and limitations](IMPLEMENTATION_STATUS.md): implemented components and their validity boundaries.
- [Verification](TEST_STRATEGY.md): what each test type establishes.
- [Reproducibility](REPRODUCIBILITY.md): configuration, saved-source and manifest requirements.
- [Thesis-first repository curation](THESIS_CLEANUP.md): reversible removal of nonessential historical material with safeguards for scientific inputs.

## Results and reports

- [Numerical results](RESULTS.md): numerical studies and their stated assumptions.
- [Regenerating reports](REGENERATE_CURRENT_REPORTS.md): regenerating presentation from existing saved results.

The default reading path is **model → calculation → numerical implementation → limitations → results**. Individual research utilities, experimental derivations and saved studies do not define additional public-facing versions of the pipeline. Their names may still appear in internal code and provenance records for reproducibility.

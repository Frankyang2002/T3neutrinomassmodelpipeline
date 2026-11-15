# Software architecture and calculation flow

The project separates field-theory construction, effective-theory evolution, numerical evaluation and reporting. A single software entry point is not evidence that all interaction bases and matching prescriptions have been combined into a complete physical calculation.

## Physical calculation

Let \(F\) denote all heavy T3 fermions and \(S\) the complete heavy scalar sector. The required hierarchy is

\[
\mathcal L_{\mathrm{UV}}(\mathrm{SM},F,S)
\xrightarrow{\text{integrate out all }F}
\mathcal L_{\mathrm{scalar\ EFT}}(\mathrm{SM},S)
\xrightarrow{\text{integrate out all }S}
\mathcal L_{\mathrm{SMEFT}}.
\]

Running is performed inside each effective theory, with matching at the relevant thresholds. Below the scalar threshold, the Weinberg coefficient is evolved and converted into the light-neutrino mass matrix using the repository's convention

\[
m_\nu=-\frac{v^2}{2}\,C_5.
\]

Masses are extracted by Takagi factorisation, and mixing observables require a consistent charged-lepton flavour basis.

## Modules and data flow

| Responsibility | Principal locations | Inputs and outputs |
| --- | --- | --- |
| Models and field representations | `common/T3Model.py`, `common/T3Fields.py`, `model/` | Representation and hypercharge choices; field inventories |
| UV Lagrangian and contractions | `Lagrangian/model/`, `Lagrangian/interactions/` | Allowed interactions, Clebsch–Gordan factors and conventions |
| Native symbolic matching | `Lagrangian/RunMatching.wl`, `Lagrangian/` | Matchete expressions and effective operators |
| Group factors and RGEs | `RGE/group_factors/`, `RGE/running/` | Representation-dependent beta functions and Wilson evolution |
| Numerical running | `Numerical/core/`, `Numerical/running/` | Scale-dependent coupling states and effective coefficients |
| Parameter optimisation | `Numerical/fitting/` | Candidate Yukawas, fit residuals and diagnostics |
| Neutrino observables | `physics/` | \(m_\nu\), masses, mass splittings and mixing observables |
| Reporting | `Reports/`, `Numerical/analysis/` | Tables, figures and derived comparisons |
| Tests | `tests/`, `validation/` | Explicit numerical and algebraic checks |

The code includes specialised source exporters, experimental calculations and saved-data processing utilities. These should be used according to their documented input requirements, not as interchangeable substitutes for the standard entry points.

## Entry points

Use `python -X utf8 pipeline.py --help` to inspect the combined pipeline options; `python -X utf8 full16_symbolic.py --help` for the symbolic 16-model study. Numerical scans with existing symbolic exports use `python -X utf8 -m Numerical.fitting.ScanCLI <existing-config.json> --dry-run`. The fast regression suite is `./scripts/RunCoreTests.ps1` in PowerShell.

These commands are **interface checks**, not proofs that the general one-loop matching and numerical fits are complete for every model.

## Result provenance

A reported numerical observable must be linked to:

1. The model representation and complete set of active UV interactions, including assumptions setting allowed couplings to zero.
2. The input configuration, fermion/scalar masses, initial scale, matching scales, coupling conventions and fitted Yukawa entries.
3. The symbolic exports used for UV RGEs, fermion matching, scalar EFT running and scalar matching, with their source fingerprints where available.
4. The perturbative order, any omitted Wilson operators and any restricted or candidate matching prescriptions.
5. Solver settings, run status, numerical checks and figure/table generation command.

A completed numerical scan is evidence for the stated approximation only. In particular, the general symbolic operator catalogue is not by itself a fully connected general-EFT numerical computation. See [capabilities](IMPLEMENTATION_STATUS.md), [verification](TEST_STRATEGY.md), and [reproducibility](REPRODUCIBILITY.md).

## Stable scientific interfaces

Keep the existing import paths of Python and Wolfram modules, the identity of saved source exports and the structure of run manifests stable unless a change has been dependency-audited. File names used by internal programs may reflect prior research terminology, but the public model is defined by its **field content, operators and thresholds**, not those names.

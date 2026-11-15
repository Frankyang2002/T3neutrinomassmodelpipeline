# Capabilities and validity boundaries

This page describes what the inspected repository source supports. It distinguishes executable software paths from physically complete matching calculations. A test or numerical run is interpreted only within its specified operator basis, input data and approximation.

## Capability matrix

| Component | Present in the source | Scope of evidence |
| --- | --- | --- |
| Ordinary T3 catalogue | Representation and hypercharge catalogue for 16 ordinary models | Model selection and structural checks; not a neutrino-mass fit |
| General allowed-interaction construction | Wolfram builders and extended scalar/Yukawa candidates | Some interactions require explicit options; full symmetry completeness is not certified by the builders |
| Renormalisable RGEs | RGBeta generators, stored exports and numerical evaluators | Generated source and selected trajectories; verify basis and coupling map for each run |
| Fermion threshold matching | Matchete matching routines and dedicated fermion-first experiments | Different routines cover different field bases or restricted cases; not interchangeable |
| Scalar intermediate EFT | Numerical evolution and candidate effective-operator mixing tools | Some Wilson-anomalous-dimension and boundary prescriptions remain conditional |
| Scalar threshold and SMEFT | Scalar matching components, Weinberg running and neutrino observable extraction | Some paths include selected one-loop contributions rather than the complete general matching |
| Numerical benchmarks and optimisation | Saved configurations, scans and fitted Yukawa routines | A fit establishes agreement for its configured approximation; it is not an independent prediction |
| Regression testing | Core and extended suites | Validates specified implementation and interface properties, not every physical one-loop term |

**Overall boundary:** The current source does not demonstrate one verified executable calculation combining the **full general renormalisable Lagrangian**, all necessary one-loop fermion matching, complete scalar-EFT operator evolution, complete scalar threshold matching and low-energy fitting for every T3 model.

## Required EFT definition

```text
UV: SM + all T3 heavy fermions and scalars
  -> integrate out all heavy fermions
Scalar EFT: SM + remaining heavy scalars + induced operators
  -> integrate out all heavy scalars
SMEFT: SM fields + Weinberg coefficient and relevant operators
  -> evolve to the low scale and extract observables
```

The source may expose alternative threshold plans or one-step modes. Such options must not be interpreted as implementation of the above calculation without checking the actual generated operators and transitions.

## How to assess one calculation

A calculation can be described as **general** only when all symmetry-allowed interactions relevant to the chosen perturbative order have been represented with the correct mass, tensor, charge and flavour conventions. It can be described as **one-loop matched** only when threshold contributions, counterterms, operator mixing, hard/soft subtraction and inherited terms appropriate to that order are included. It is **numerically verified** only for the configurations and checks actually run.

Special care is required for B(−1) and C(−1), which admit additional mass and interaction structures and potentially independent flavour contributions to the Weinberg coefficient. Restricted couplings can be legitimate benchmark choices, but the restriction must be recorded rather than silently described as a general calculation.

For implementation ownership, see [Software architecture](PROJECT_GUIDE.md). For results, see [Numerical results](RESULTS.md) and the per-run JSON manifests. For tests, see [Verification](TEST_STRATEGY.md).

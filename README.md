# T3 radiative neutrino-mass pipeline

This repository constructs gauge-invariant T3 neutrino-mass models, generates symbolic Lagrangians and matching/RGE inputs, and evaluates neutrino masses and mixing parameters. The reference theory is the **complete renormalisable interaction basis allowed by the specified gauge representations and discrete symmetry**. A selection of coupling values is a benchmark, not a different fundamental Lagrangian.

## Technical reference

- [Model, interactions and EFT sequence](docs/MODEL_AND_PIPELINE.md): fields, charge conventions, operator-basis requirements, matching sequence and implementation boundaries.
- [Repository maintenance](docs/REPOSITORY_CLEANUP_POLICY.md): how to audit and archive historical material without breaking source generation or reproducibility.
- `Lagrangian/model/LagrangianBuilder.wl` and `Lagrangian/interactions/`: Wolfram model construction and candidate invariants.
- `common/PipelinePlan.py`, `common/Thresholds.py`, `pipeline.py`: model/EFT orchestration.
- `RGE/`, `Numerical/`, `physics/`: renormalisation, numerical evolution and neutrino observables.
- `configs/` and `tests/`: inputs and regression checks.

## Physical EFT ordering

```text
UV: SM + all T3 fermions + all T3 scalars
   | run full-theory parameters and match at fermion thresholds
   v
Intermediate EFT: SM + all T3 scalars + fermion-induced operators
   | run intermediate-EFT Wilson coefficients and couplings
   | match at scalar thresholds
   v
SMEFT: SM fields + dimension-five Weinberg operator (within the selected truncation)
   | run C5 to the low-energy scale
   v
Majorana neutrino mass matrix and oscillation observables
```

**All heavy fermions are integrated out before any heavy scalar.** Different masses within each group may require multiple ordered thresholds. A simultaneous all-heavy-field threshold is a separately supported comparison, not the same calculation as the hierarchy above.

## Generality and reproducibility

There is a distinction between *constructing* an inventory of symmetry-allowed operators and *certifying* that their physical-basis normalisations, one-loop RGEs, threshold matching and numerical treatment are complete. The current repository has additional-interaction switches and separate general-symbolic/restricted-production routes. Their existence must not be interpreted as a successful general-theory numerical calculation. Consult the technical reference for explicit evidence and open checks.

Useful commands from the repository root (Windows PowerShell):

```powershell
$env:PYTHONUTF8 = '1'
python -X utf8 pipeline.py --help
python -X utf8 -m scripts.RunT3GeneralSymbolicStudy --self-test
python -X utf8 -m pytest -q --tb=line
python -X utf8 scripts/maintenance/ExtremeRepositoryCleanup.py
```

The general symbolic self-test checks the source audit utility; it is **not** a live Matchete/RGBeta or end-to-end physics validation. Model-specific symbolic runs require the installed Wolfram/Matchete/RGBeta environment and the correct physical-basis setup.

Do not silently use restricted-basis output files as if they were obtained from the general Lagrangian. Inspect source manifests and the coupling/threshold content before any quantitative comparison.

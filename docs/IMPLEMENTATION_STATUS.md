# Implementation status and scientific scope

This page separates **source implementations**, **previously reported executions**, **conditional calculations** and **unresolved physics**. It is based on read-only inspection of Git commit `806cdf08608dbe3f5e0dcd8af785cfad16f26af6` (8 November 2026). No new Wolfram calculation, full numerical rerun or independent matching certification was performed to prepare this status page.

## Target EFT calculation

```text
UV: SM + all heavy T3 fermions + all heavy T3 scalars
  -> integrate out ALL heavy fermions
Intermediate EFT: SM + heavy scalars + induced effective operators
  -> run the supported couplings/Wilson structures
  -> integrate out ALL heavy scalars
SMEFT: SM fields + dimension-five Weinberg operator (selected truncation)
  -> low-energy RGE -> Majorana mass matrix -> neutrino observables
```

For multiple masses within a group, their threshold ordering and mass-basis treatment require explicit matching. The ordinary production plan in `common/PipelinePlan.py` permits a common threshold or the supported **fermion-first, grouped-scalar** hierarchy; arbitrary scalar-first paths and fully general mixed-scalar sequential decoupling are not established as production calculations.

The conventions include `Q = T3 + Y`, the one-loop beta numerator `16*pi^2 dX/dln(mu)`, and `m_nu = -v^2 C5/2` in the stated Weinberg convention. See [physics](PHYSICS.md) and [numerical methods](NUMERICAL.md).

## Source-verified software inventory

| Component | What is implemented in the checked source | What the implementation does **not** establish |
| --- | --- | --- |
| Model catalogue and field selection | `common/T3Model.py`, `common/T3Fields.py`, `model/T3Study.py` and Wolfram catalogues describe the ordinary 16-model representation/hypercharge set | Physical completeness of every independent interaction coefficient |
| General UV interaction construction | `Lagrangian/` and `scripts/RunT3GeneralSymbolicStudy.py` construct/check symmetry-allowed additional interactions and export source-fingerprinted symbolic results | A certified Matchete-to-RGBeta action-level coupling dictionary |
| UV and scalar-only renormalisable running | `RGE/running/rgbeta/`, `Numerical/running/UVRunner.py` and `IntermediateScalarRunner.py` provide exporters/evaluators and SciPy ODE evolution | Complete *general* EFT Wilson-coefficient RGEs or all mixed scalar mass-basis mappings |
| Threshold planning | `common/PipelinePlan.py`, `common/Thresholds.py` and `pipeline.py` represent physical heavy-field removal and reject unsupported configurations | A physics-certified calculation merely because a threshold plan can be represented |
| Restricted numerical Weinberg coefficient | `Numerical/running/FinalC5TrajectoryAdapter.py` supplies source-derived matching inputs to the existing final-C5 evaluator | A new independent C5 formula or all-general-coupling sequential matching |
| Low-energy neutrino calculation | `physics/NeutrinoMass.py`, `physics/NeutrinoObservables.py`, and `Numerical/running/SMWeinbergEvolution.py` handle numerical `C5`, mass matrices, running and mixing | Prediction of fitted oscillation targets independent of fitted Yukawa inputs |
| Fitting and diagnostics | `Numerical/fitting/ScanCLI.py`, `ParameterScan.py` and plotting/report scripts provide benchmark scans and presentation | A unique best-fit solution or automatic experimental exclusion of other models |
| Experimental B/C fermion-first probes | `scripts/RunGeneralT3FermionFirstEFT1Tree.py` and `...OneLoop.py` define live tree and raw one-loop fermion-threshold checks with scalar eigenstates retained | A renormalised, certified physical C5 or complete two-threshold general-EFT result |

The general-interaction audit inventories **50 additional quartic model-slots, 8 additional scalar-bilinear model-slots and 6 extra Yukawa interaction types across 16 models**, as encoded in `scripts/AuditFull16GeneralSymbolicPhysics.py` and the saved inventory documentation. These are totals *across different models*, not simultaneous couplings in every model. Inventory agreement and beta-key coverage alone do not prove physical matching normalisation.

## Execution evidence (previously reported, not rerun here)

| Evidence | Reported result | Appropriate interpretation |
| --- | --- | --- |
| Core regression suite, 8 November 2026 | **189 passed** in the user's Windows checkout | Selected software/physics interface regressions, not end-to-end general-theory certification |
| One-point operational numerical smoke test | `PASS_END_TO_END_OPERATIONAL` for ordinary T3-B(-1), using existing exported inputs | Demonstrates a numerical execution path for that benchmark, not clean source regeneration |
| Reduced-potential frozen-boundary study | **32/32 individually refitted cases accepted** in two input scenarios, per saved `docs/RESULTS.md` | Calibration of selected neutrino observables within a restricted and nonuniform fixed-order matching setup |
| General-scalar exploratory comparisons | Saved conditional statuses summarised in `docs/RESULTS.md` | Preliminary numerical studies; not validated general-theory neutrino predictions |
| Last reported extended pytest run | **4,453 passed, 73 failed, 51 skipped, 1 error** before subsequent targeted infrastructure corrections | Evidence of unresolved integration/data/test issues; this result is not replaced by the green core suite |

The last extended-suite result is historical. The current source's full suite was **not** rerun for this page. Some failures are absent native certificates or stale interfaces, while others are genuine numerical/source-convention discrepancies that must be investigated separately; see [test triage](TEST_FAILURE_TRIAGE_2026_11_08.md).

## Explicitly incomplete or conditional physics

1. **General native coupling dictionary:** quartic tensor bases, action prefactors, signs, Hermitian conjugation and fermionic conventions must agree between Matchete and RGBeta in the physical mass basis. A symbolic inventory is not that certificate.
2. **General sequential matching:** integrating out all fermions and subsequently all scalars, including physically mixed B/C scalar eigenstates, is not yet a complete independently certified one-loop general-theory Wilson calculation. The checked B/C fermion-first scripts explicitly mark their tree/raw one-loop outputs as **not physical C5**.
3. **EFT Wilson evolution:** a general renormalisable UV/EFT1 beta inventory does not demonstrate full Wilson-basis anomalous dimensions, finite threshold pieces or pole cancellations at the requested perturbative order.
4. **Source-dependent numerical reference values:** discrepancies between saved numerical results and current evaluator outputs must be understood before relabelling them verified. Test failures must not be removed or output certificates fabricated to assert completeness.
5. **Restricted potential:** setting additional symmetry-allowed quartics/bilinears to zero at a named boundary is a benchmark choice, not a symmetry prohibition or statement about arbitrary nonzero couplings.

## Reading and verification order

- [Install and run the software](INSTALLATION.md).
- [Model and calculation contract](MODEL_AND_PIPELINE.md): physical fields, invariants and EFT ordering.
- [General symbolic physics audit](GENERAL_SYMBOLIC_PHYSICS_AUDIT.md): what the 16-model inventory checks actually cover.
- [Scalar sector](SCALAR_SECTOR.md): operator and boundary assumptions.
- [Results](RESULTS.md): fit statuses, figures, numerical conventions and saved-data provenance.
- [Testing strategy](TEST_STRATEGY.md): core versus extended tests and live symbolic evidence.

**Project interpretation:** the current repository provides substantial symbolic machinery, restricted numerical T3 benchmark studies and targeted general-theory research probes. It should **not** be advertised as a fully validated, general-coupling, two-threshold neutrino-mass prediction code until the matching and validation requirements above are satisfied.

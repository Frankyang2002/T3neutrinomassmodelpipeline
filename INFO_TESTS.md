# T3 Test Suite — File-by-File Guide

This document explains the problem each test file is protecting and what it checks.

The suite has four roles: physics contracts, numerical contracts, architecture contracts, and external Wolfram/Matchete/RGBeta regressions.

## File-by-file map

| File | Problem solved/protected | What it does |
|---|---|---|
| `tests/test_benchmark_sensitivity.py` | Controlled UV scaling directions | Checks T1/T2/T3 rescalings preserve the leading T3 neutrino-mass prefactor. |
| `tests/test_best_fit_diagnostics.py` | Best-fit diagnostic serialization | Checks optimizer-best loading, complex-matrix JSON round trips and logarithmic save scales. |
| `tests/test_beta_vector.py` | Beta dictionary to ODE vector mapping | Checks layout, complex packing, one-loop factor application and beta-key/type validation. |
| `tests/test_default_result_figure_set.py` | Focused default figure set | Ensures primary thesis/result plots remain default while diagnostic-only plots remain optional. |
| `tests/test_eft_foundations.py` | Generic EFT data model | Checks d<=5 truncation, threshold steps, stage metadata and shared-scalar physical/formal field handling. |
| `tests/test_eft_running_intervals.py` | Running-interval construction | Checks common-threshold, fermion-first, scalar-first representability, shared-scalar intervals and adjacency rules. |
| `tests/test_eft_transitions.py` | Physical EFT transitions | Checks integrating out active fields, truncation propagation and consistent ordinary/shared final EFT content. |
| `tests/test_final_c5_bridge.py` | Final-C5 threshold/basis preparation | Checks threshold values, evaluator schema, vectorlike SVD, Majorana Takagi rotation and canonical evaluator delegation. |
| `tests/test_final_c5_contribution_diagnostics.py` | Hard/direct final-C5 decomposition | Checks C5_final = C5_hard + C5_direct and rejects inconsistent diagnostics. |
| `tests/test_final_c5_trajectory_adapter.py` | Canonical adapter callable | Checks the public callable delegates to the canonical trajectory evaluator. |
| `tests/test_final_c5_trajectory_adapter_architecture.py` | Final-C5 adapter ownership | Checks the canonical implementation exists, old bridge is retired and compatibility aliases remain. |
| `tests/test_final_canonical_architecture.py` | Canonical caller ownership | Prevents internal code/docs from reverting to retired wrapper paths. |
| `tests/test_final_refactor_contract.py` | Whole-refactor architecture contract | Checks retired root wrappers, bridge removal, docs, hypercharge convention, scalar-first scope and script placement. |
| `tests/test_final_sm_boundary.py` | Scalar-threshold numerical boundary | Checks scalar masses, SM projection, scale consistency, full Yukawa retention and Weinberg initial conditions. |
| `tests/test_final_weinberg_matching_ownership.py` | Final-C5 ownership | Ensures authoritative final-Weinberg construction lives under RGE/matching, not running. |
| `tests/test_full_flavor_weinberg.py` | Full-flavor final EFT numerics | Checks complex Yukawa packing, diagonal limit, charged-lepton basis rotation and mass normalization. |
| `tests/test_full_flavor_weinberg_stage_architecture.py` | Full-flavor symbolic stage ownership | Checks the descriptive canonical stage and removal of the historical stage. |
| `tests/test_full_study_architecture.py` | 16-model study orchestration | Checks the model grid, four comparison scenarios and forwarding/execution structure. |
| `tests/test_intermediate_eft_running_architecture.py` | Intermediate production dispatch | Ensures dispatch is backend-agnostic, physical-content based and free of validation-only transport. |
| `tests/test_intermediate_eft_validation_architecture.py` | Independent intermediate validation | Checks validation remains separate from production and owns independent transport/pole checks. |
| `tests/test_intermediate_running_interfaces.py` | Physical intermediate interfaces | Checks canonical ownership of renormalisable running, tensor RGE, direct Weinberg and threshold resume. |
| `tests/test_intermediate_scalar_boundary.py` | Full multi-segment numerical path | Checks UV/intermediate/scalar boundaries, diagnostics, save scales and final continuation. |
| `tests/test_intermediate_scalar_runner.py` | Scalar-only ODE runner | Checks scale direction, save points, endpoint reconstruction and nontrivial running. |
| `tests/test_intermediate_scalar_vector.py` | Scalar-only state vector/evaluator | Checks packing, beta-key agreement, derivative shape and stage metadata. |
| `tests/test_intermediate_weinberg_diagnostics.py` | Intermediate direct-Weinberg diagnostics | Checks the sampled diagnostic stores the expected magnitude trajectory. |
| `tests/test_lagrangian_runner_architecture.py` | Python/Wolfram Lagrangian boundary | Protects separation of validation, subprocess execution, result ingestion and stable public Runner API. |
| `tests/test_legacy_eft1_source_cleanup.py` | Retired EFT1 source architecture | Ensures old EFT1 modules/imports remain absent and canonical stack imports cleanly. |
| `tests/test_low_energy_neutrino_architecture.py` | Post-matching neutrino orchestration | Checks pipeline stage order while detailed low-energy implementation remains in physics/. |
| `tests/test_markdown_math_delimiters.py` | Documentation math rendering | Checks repository Markdown uses the chosen dollar-delimiter convention. |
| `tests/test_matched_weinberg_ownership.py` | Matched-Weinberg parser/benchmark ownership | Separates Matchete parsing, one-generation benchmark algebra, stage output and lazy compatibility wrapper. |
| `tests/test_neutrino_mass_architecture.py` | Symbolic C5->mnu ownership | Ensures mass conversion lives in physics/NeutrinoMass.py rather than RGE stages. |
| `tests/test_neutrino_mass_convention.py` | Physical neutrino-mass normalization | Checks m_nu = -(v^2/2) C5, symbolic/numeric agreement and symmetry. |
| `tests/test_neutrino_observables_architecture.py` | Observable ownership | Ensures Takagi/PMNS phenomenology lives under physics/ and historical RGE phenomenology path stays removed. |
| `tests/test_numeric_neutrino_mass_ownership.py` | Numerical C5->mnu ownership | Prevents duplicated numerical mass-conversion implementations outside physics/. |
| `tests/test_numerical_rge_ownership.py` | RGE model vs numerical solver separation | Keeps analytic beta functions in RGE/ and ODE integration in Numerical/. |
| `tests/test_numerical_state.py` | Canonical UV-state validation | Checks ordinary/shared states and representation-dependent optional quartics. |
| `tests/test_one_generation_weinberg_benchmark_architecture.py` | One-generation benchmark role | Ensures the benchmark is analytic/no-I/O and distinct from physical full-flavor running. |
| `tests/test_oscillation_fit.py` | Oscillation chi-square layer | Checks canonical observables, covariance use, target loading and ordering consistency. |
| `tests/test_oscillation_optimizer.py` | Local fit optimizer | Checks coordinate transforms, sensitivity ranking, warm starts and toy-fit improvement. |
| `tests/test_parameter_scan.py` | Generic scan machinery | Checks deterministic grid/random generation, failures, callbacks, best points and JSON output. |
| `tests/test_pipeline_backbone.py` | Central pipeline readability/ownership | Keeps pipeline.py compact, ordered and dependent on delegated configuration/physics stages. |
| `tests/test_pipeline_cli.py` | User-facing run configuration | Checks ordinary/shared modes, production threshold restrictions, study naming and numerical CLI behavior. |
| `tests/test_pipeline_numerical_default.py` | Default integrated numerical config | Checks schema, default scales, Sobol search and sensitivity settings. |
| `tests/test_pipeline_plan.py` | Single physical execution plan | Checks threshold scales, stage labels, d<=5 metadata, shared-scalar naming and production restrictions. |
| `tests/test_pipeline_reports_architecture.py` | Report orchestration boundary | Keeps detailed reporting out of pipeline.py while preserving output order/contracts. |
| `tests/test_plot_optimizer_results.py` | Optimizer-result plotting | Checks expected convergence/pull/sensitivity files are generated from stored JSON. |
| `tests/test_rgbeta_evaluator.py` | RGBeta InputForm numerical parser | Checks traces/matrices/conjugation, gauge convention, payload evaluation and representation matching. |
| `tests/test_rgbeta_metadata_compatibility.py` | RGBeta metadata backward compatibility | Allows legacy ordinary metadata while rejecting explicit shared-state mismatch. |
| `tests/test_running_diagnostics.py` | Running diagnostic formulas | Checks descending save scales and PMNS-magnitude to mixing-angle conversion. |
| `tests/test_running_result_figures.py` | Primary running figures | Checks the standard stored-diagnostic figure set can be generated. |
| `tests/test_scalar_first_safety.py` | Scalar-first production safety | Checks scalar-first remains representable but is rejected for current d<=5 production. |
| `tests/test_scan_cli.py` | Config-driven scan CLI | Checks nested bindings, relative paths, representation/state construction and stage metadata. |
| `tests/test_sm_weinberg_numerical_naming.py` | Final SM+Weinberg canonical naming | Prevents return of ambiguous historical numerical Weinberg wrappers. |
| `tests/test_sobol_benchmark_search.py` | Automatic Sobol benchmark search | Checks reproducibility, bounds, 19-real-parameter profile, retargeting and local refinement. |
| `tests/test_state_vector.py` | UV state pack/unpack | Checks ordinary/Majorana/shared round trips and vector-size validation. |
| `tests/test_t3_field_content.py` | Physical heavy-field identity | Checks ordinary F,S1,S2 versus shared F,S and formal expansion only at matching. |
| `tests/test_t3_matching_architecture.py` | T3 matching boundary | Ensures pipeline uses the explicit matching module and that it is the sole Runner boundary. |
| `tests/test_t3_model_selection.py` | Neutral T3 scan | Checks the expected 16 ordinary neutral points and neutrality/topology behavior. |
| `tests/test_t3_study_architecture.py` | Study selection ownership | Keeps study selection independent of Lagrangian execution and preserves comparison definitions. |
| `tests/test_t3_trajectory.py` | Complete numerical trajectory | Checks threshold diagnostics, save scales, ordering and final Weinberg continuation. |
| `tests/test_thesis_result_figures.py` | Thesis result rendering | Checks expected table sections and final figure/table output files. |
| `tests/test_threshold_boundary_naming.py` | Boundary vs matching terminology | Ensures numerical threshold modules describe projections/boundaries rather than reimplement matching. |
| `tests/test_uv_runner.py` | UV numerical ODE | Checks endpoint, save scales, state reconstruction, direction and nontrivial running. |
| `tests/test_uv_running_architecture.py` | UV stage ownership/order | Checks one named UV production stage and its machine-output contract. |
| `tests/test_weinberg_flavor_matching_architecture.py` | Flavor-lifting ownership | Separates one-generation flavor lifting from hierarchical final-C5 adaptation. |
| `tests/test_weinberg_trajectory.py` | Scale-dependent SM+Weinberg trajectory | Checks saved scales, C5 symmetry, mass normalization, observables and running. |
| `tests/test_weinberg_trajectory_architecture.py` | Trajectory integration vs interpretation | Keeps ODE integration in Numerical/ and physical interpretation in physics/. |
| `tests/python/check_t3_scalar_quartics.py` | Independent scalar-quartic completeness check | Counts SU(2) singlets/quartic monomials and compares with implemented potential structures. |
| `tests/python/test_rgbeta_t3_running.py` | Standalone RGBeta regression | Exercises/checks the external T3 RGBeta running output contract. |
| `tests/python/test_shared_scalar_mode.py` | Standalone shared-scalar regression | Checks dimensions, physical/formal threshold expansion and one-scalar RGE model counting. |
| `tests/python/test_t3_scotogenic_normalization.py` | Scotogenic normalization benchmark | Checks the critical T3-B alpha=-1 low-dimensional normalization limit. |
| `tests/python/test_weinberg_normalization_cleanup.py` | Weinberg/Wilson prefactor regression | Checks canonical term scaling, C12 flavor symmetry and prefactor preservation. |
| `tests/reference/MaUVRGE.py` | Independent Ma/scotogenic numerical reference | Provides separate UV RGE, trajectory, Takagi and C5 matching for comparison. |
| `tests/reference/T3RGEComponentExport.wl` | Exact component-export reference | Serializes exact T3 tensor component factors and invariant-family metadata. |
| `tests/reference/T3RGEPotentialExport.wl` | Scalar-potential tensor reference | Exports reference Higgs/scalar quartic tensor entries. |
| `tests/reference/T3RGETensorExport.wl` | Wolfram-to-Python tensor reference bridge | Exports scalar metadata, quartics, mixing tensor and raw Yukawa invariant components. |
| `tests/wolfram/ProbeT3SU2Factors.wl` | Isolated SU(2) factor diagnostic | Contracts T3 invariant tensors with internal indices open/closed as needed without loop matching. |
| `tests/wolfram/RegressionC5.wl` | End-to-end external C5 regression | Reads smoke outputs, applies convention/normalization checks and writes the final regression report. |
| `tests/wolfram/T3CouplingConventions.wl` | Generalized-vs-legacy coupling conventions | Stores/converts measured normalization and phase factors for low-dimensional benchmarks. |
| `tests/wolfram/TestT3Conventions.wl` | Representation/CG convention regression | Consolidates exact representation, duality, CG and API convention tests. |
| `tests/wolfram/TestT3RGBeta.wl` | Strict external RGBeta regression | Checks the low-dimensional T3 RGBeta model, closure and exported beta structure. |
| `tests/wolfram/TestT3RGEExport.wl` | RGE tensor-export regression | Checks generalized T3 tensor exports for representative simple and triplet-rich models. |
| `tests/wolfram/TestWeinberg.wl` | Direct Weinberg normalization regression | Projects actual matched LLHH terms onto the physical neutrino/Higgs-neutral component. |

## How to read the suite

Tests with `architecture` in the name mainly protect ownership boundaries and retired modules. They are intentionally different from numerical/scientific tests: a physics calculation can still be numerically correct while living in the wrong layer, which makes future maintenance much harder.

The standalone files under `tests/python/`, `tests/reference/` and `tests/wolfram/` provide additional independent or external checks. In particular, the reference implementations are useful because they avoid simply comparing production code with itself.

## Recommended regression levels

Fast Python suite:

```powershell
python -m pytest -q
```

Real external smoke calculation:

```powershell
python pipeline.py --smoke
```

Full external C5 regression:

```powershell
python scripts/run_regression.py
```

## Overall test-suite overview

```text
small unit tests
    catch local algebra/state-layout mistakes

trajectory/fitting tests
    catch numerical composition errors

architecture tests
    stop ownership boundaries from regressing

reference implementations
    provide independent numerical/tensor comparisons

Wolfram tests
    validate the external symbolic/group-theory stack

smoke/regression runs
    verify the real end-to-end external pipeline
```

The central principle is:

> A passing unit test proves that a local interface behaves as expected; a passing external regression proves that the actual Matchete/RGBeta calculation still produces the expected physics.

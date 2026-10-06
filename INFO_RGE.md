# T3 RGE, Threshold Transport and Final C5 — File-by-File Guide

This document explains every implementation file under `RGE/`.

The subsystem solves four distinct problems:

1. represent generic gauge/scalar/fermion structures;
2. derive/validate representation-dependent group factors;
3. run the UV and intermediate EFTs;
4. construct and run the final SM + Weinberg theory.

The numerical ODE implementations are documented separately in `INFO_NUMERICAL.md`.

---

# 1. Conventions and physical flow

The project uses

$$
t=\ln\mu,\qquad Q=T_3+Y,
$$

and

$$
m_\nu=-\frac{v^2}{2}C_5.
$$

Production flow:

```text
UV T3
    -> UV renormalisable RGE
    -> F threshold
    -> scalar-only intermediate EFT
    -> dimension-five Wilson transport
    -> direct LLSS -> O5 mixing
    -> scalar threshold
    -> final C5
    -> SM + Weinberg RGE
```

At fixed one-loop order, self-running of the intermediate LLSS coefficient followed by another scalar loop would first contribute at two loops, so that effect is not promoted into the authoritative one-loop final $C_5$.

---

# 2. `RGE/general/` — reusable mathematical building blocks

## `RGE/general/AnomalousDimensions.py`

**Problem solved:** generic Wilson RGEs repeatedly need scalar/fermion collinear anomalous-dimension contributions.

**What it does:** computes those universal anomalous-dimension pieces and attaches them to tensor-RGE inputs.

## `RGE/general/FermionBasis.py`

**Problem solved:** fermion multiplets with flavor and gauge indices need one explicit basis before gauge generators can act.

**What it does:** defines Weyl-fermion basis blocks, global fermion bases, embedded SU(2)/U(1) generators and gauge sectors.

## `RGE/general/GaugeGenerators.py`

**Problem solved:** arbitrary SU(2) representations require consistent generator normalisation in complex and real bases.

**What it does:** constructs SU(2) generators, real-scalar generators, embedded global generators, U(1) generators, gauge sectors and quadratic Casimir matrices.

## `RGE/general/ScalarBasis.py`

**Problem solved:** complex scalar multiplets need a deterministic global real-component basis for tensor RGEs.

**What it does:** defines scalar basis blocks and the global scalar basis.

## `RGE/general/WilsonTensorRGE.py`

**Problem solved:** dimension-five intermediate operators require a general tensor beta-function engine rather than model-specific component formulae.

**What it does:** builds wavefunction, scalar-pair, mixed gauge-Yukawa and crossed-Yukawa pieces of the Wilson tensor RGE and validates dimensions/permutations.

---

# 3. `RGE/group_factors/core/` — invariant/tensor algebra

## `MixingQuarticTensorAlgebra.py`

**Problem solved:** representation-dependent mixing quartics exported by Matchete/RGBeta need to be projected onto a physical invariant direction.

**What it does:** builds basis tensors, computes tensor inner products, constructs the cross tensor and projects tensors onto physical directions.

## `QuarticTensorAlgebra.py`

**Problem solved:** sparse scalar-quartic tensors need common low-level algebra.

**What it does:** converts polynomials to tensors, evaluates weighted tensor inner products, performs pair maps, sparse contractions and sector restrictions.

## `RepresentationFactors.py`

**Problem solved:** repeated representation factors such as $C_2$ and canonical Yukawa leg factors need one convention owner.

**What it does:** computes SU(2) quadratic Casimirs and canonical Yukawa-leg factors from representation dimension.

## `ScalarQuarticBasis.py`

**Problem solved:** Matchete's representation-dependent quartic naming must be mapped onto the physical singlet/adjoint/cross basis, and the invariant basis must be checked for completeness.

**What it does:**

- constructs physical scalar-quartic tensors;
- determines Matchete sector names;
- solves the basis transformation;
- checks tensor residuals/rank;
- validates the expected number of mixed scalar invariants.

It is a research/validation driver, not the production ODE integrator.

## `YukawaGroupFactors.py`

**Problem solved:** scalar-leg wavefunction factors from the T3 Yukawa tensors must be extracted in exact representation-dependent form.

**What it does:** loads exported invariant tensors, contracts the scalar leg and verifies the result has the expected Schur/identity structure.

---

# 4. `RGE/group_factors/recoupling/` — translate between invariant bases

## `MixingQuarticRecoupling.py`

**Problem solved:** the T3 mixing-quartic coefficient can be represented in different tensor bases by Matchete and RGBeta.

**What it does:** computes the projection/recoupling between those bases and reports exact residuals.

## `PortalQuarticGroupFactors.py`

**Problem solved:** Higgs-scalar and scalar-scalar portal beta functions require representation-dependent singlet/adjoint coefficients.

**What it does:** derives portal quartic beta-function group factors.

## `PortalRecouplingDerivations.py`

**Problem solved:** cross-portal tensor structures involving conjugate representations require explicit dual-map/charge-conjugation algebra.

**What it does:** constructs relevant generators, metrics and cross tensors and decomposes them into the physical portal basis.

## `RGBetaMixingQuarticRecoupling.py`

**Problem solved:** an independent recoupling calculation is needed directly in the RGBeta representation convention.

**What it does:** reconstructs generators/polynomials/tensors and projects RGBeta cross-beta structures onto the T3 physical basis.

---

# 5. `RGE/group_factors/validation/` — independent group-factor checks

## `ValidateDirectWeinbergGroupFactors.py`

**Problem solved:** direct LLSS $\to O_5$ mixing carries an SU(2) group factor that must match the model's invariant tensors.

**What it does:** derives the factor from T3 representation data, parses diagnostic expressions and validates direct-Weinberg rows.

## `ValidateMassGroupFactors.py`

**Problem solved:** fermion/scalar mass beta coefficients from RGBeta need independent representation-level checking.

**What it does:** extracts gauge/Yukawa coefficients from stored beta expressions and compares them with expected group-theory values.

## `ValidateMixingQuarticGroupFactors.py`

**Problem solved:** the $\lambda_{T3}$ beta coefficients need an independent group-theory expectation.

**What it does:** computes expected mixing-quartic factors and checks the RGBeta expression.

## `ValidateRGBetaRecouplings.py`

**Problem solved:** basis conversion can hide sign/normalisation mistakes even when individual expressions look plausible.

**What it does:** constructs an independent basis calculation and checks RGBeta recoupling residuals.

## `YukawaBetaGroupFactors.py`

**Problem solved:** T3 Yukawa beta functions require representation-dependent self/cross/gauge coefficients.

**What it does:** derives and serialises the complete Yukawa group-factor set and exposes a CLI diagnostic.

---

# 6. `RGE/matching/` — Weinberg matching ownership

## `FinalWeinbergCoefficient.py`

**Problem solved:** in a hierarchy, the physical final Weinberg coefficient receives both a hard scalar-threshold contribution and direct intermediate-EFT running.

**What it does:**

- loads hard and direct pieces;
- combines them into the authoritative final $C_5$;
- performs the established pole/RGE consistency bookkeeping;
- checks flavor symmetry;
- serialises `final_weinberg_coefficient.json`.

This is **matching/bookkeeping**, not numerical RGE evolution.

## `FinalWeinbergAdapter.py`

**Problem solved:** the hierarchical final-C5 JSON must become a symbolic physical Majorana matrix.

**What it does:** parses stored hard/running/final pieces and constructs the symbolic symmetric $C_5$ matrix used downstream.

## `MatcheteC5Parsing.py`

**Problem solved:** Matchete `InputForm` Weinberg expressions are not directly SymPy expressions.

**What it does:** provides the text-level Matchete-to-SymPy parser, including coupling/bar handling.

No RGE equation belongs here.

## `WeinbergFlavorMatching.py`

**Problem solved:** one-generation Matchete matching must be lifted to the full three-generation Majorana flavor structure without changing the representation-dependent loop kernel.

**What it does:**

- factors the one-generation result into loop kernel and Yukawa factors;
- introduces symbolic $y_1,y_2$ flavor matrices;
- builds the symmetric physical flavor $C_5$;
- writes flavor-matching artifacts.

## `WeinbergTensorAdapter.py`

**Problem solved:** the generic Wilson-tensor RGE needs Standard-Model Yukawa/gauge structures in its tensor representation.

**What it does:** constructs the SM EFT tensor inputs and symmetric Yukawa insertions used by the Weinberg/Wilson tensor machinery.

## `MatchedWeinbergRGE.py`

**Problem solved:** historical callers still require the old matched-Weinberg entry point.

**What it does:** provides a lazy compatibility wrapper that routes to the current parser/one-generation benchmark stage.

**Important:** this is a real compatibility boundary and should not be removed casually.

---

# 7. `RGE/running/` — production stage dispatch

## `UVRunning.py`

**Problem solved:** the Python pipeline needs one production entry point for the UV one-loop renormalisable RGE stage.

**What it does:** orchestrates RGBeta generation/loading for the current model, runs the UV RGE stage and records machine-readable outputs/status.

The numerical `solve_ivp` implementation lives in `Numerical/running/UVRunner.py`.

## `IntermediateEFTRunning.py`

**Problem solved:** intermediate EFTs should be selected by physical active-field content, not hard-coded stage numbers.

**What it does:**

- receives an `EFTRunningInterval`;
- asks registered backends whether they support the field content;
- runs the appropriate backend;
- returns a common status/result object.

It intentionally contains no detailed scalar-only beta algebra.

---

# 8. `RGE/running/backends/`

## `ScalarOnlyAfterFermion.py`

**Problem solved:** the verified hierarchy needs one production backend for the theory after $F$ is integrated out while the scalar sector remains active.

**What it does:**

```text
renormalisable scalar-only RGE
    -> dimension-five Wilson RGE
    -> flavor Wilson boundary
    -> direct Weinberg generation
    -> scalar-threshold continuation
    -> final C5 construction
```

This is the main production intermediate-EFT backend.

## `__init__.py`

Package marker only.

---

# 9. `RGE/running/intermediate/` — scalar-only EFT implementation

## `IntermediateMatcheteParsing.py`

**Problem solved:** intermediate Matchete sparse arrays, representations and scalar expressions need robust parsing.

**What it does:** bracket matching, top-level splitting, barred-expression handling, index parsing, representation parsing and sparse-array parsing.

## `ScalarOnlyRenormalisableRunning.py`

**Problem solved:** after removing $F$, the remaining renormalisable SM+scalar couplings still run.

**What it does:** executes the scalar-only RGBeta renormalisable stage and stores its result.

## `ScalarOnlyTensorAdapters.py`

**Problem solved:** exported Matchete Wilson Clebsch-Gordan structures need to become sparse Python tensors in the scalar-only basis.

**What it does:** defines Wilson tensor/leg/CG data structures and converts the exported CG registry into the basis used by the tensor RGE.

## `ScalarOnlyWilsonTensorRGE.py`

**Problem solved:** dimension-five LLSS Wilson coefficients need one-loop tensor beta functions and closure diagnostics.

**What it does:**

- builds the scalar-only RGE context;
- finds seed and closure components;
- evaluates component beta functions;
- checks symmetry;
- diagnoses the Weinberg subspace.

## `ScalarOnlyWilsonFlow.py`

**Problem solved:** tensor beta functions must be transported between two physical scales and written in a stable component representation.

**What it does:** loads the RGE payload, canonicalises component keys and evolves/transports Wilson components.

## `ScalarOnlyWilsonRunning.py`

**Problem solved:** production orchestration needs one named stage for dimension-five scalar-only Wilson running.

**What it does:** calls the tensor-RGE/flow machinery and records the stage result.

## `FlavorWilsonBoundary.py`

**Problem solved:** tree-level matching at the fermion threshold must retain the full lepton-flavor structure of the LLSS coefficient.

**What it does:** exports the full-flavor Wilson boundary used by direct Weinberg running and final matching.

## `DirectWeinbergRunning.py`

**Problem solved:** the scalar-only intermediate EFT directly mixes LLSS into the Weinberg operator at one loop.

**What it does:**

- combines the flavor-blind component prefactor with the full $C_{12}^{pq}$ flavor tensor;
- constructs the direct $O_5$ contribution between the fermion and scalar scales;
- serialises the expression used by the authoritative final-C5 stage.

Historical `EFT1` wording that survives in external data is compatibility naming only.

## `ScalarThresholdMatching.py`

**Problem solved:** after intermediate running, Matchete must resume the stored EFT state and integrate out the remaining scalar threshold.

**What it does:**

- discovers the correct continuation artifact;
- resolves the threshold-2 model/result path;
- builds and executes the resumed Wolfram command;
- validates continuation markers.

## `__init__.py`

Package marker only.

---

# 10. `RGE/running/rgbeta/` — external RGBeta model/export layer

## `RGBetaT3Running.py`

**Problem solved:** Python needs a stable boundary to external RGBeta runs for both UV and scalar-only theories.

**What it does:** launches the relevant Wolfram runner, loads its JSON and returns structured UV/intermediate RGBeta results.

## `RunT3RGBeta.wl`

**Problem solved:** execute RGBeta for the complete UV T3 model.

**What it does:** defines the requested UV model in RGBeta, computes beta functions and exports them in the repository reporting convention.

## `RunT3ScalarOnlyRGBeta.wl`

**Problem solved:** after removing $F$, the reduced scalar-only renormalisable theory needs a separate RGBeta calculation.

**What it does:** constructs that reduced theory and exports its beta functions.

## `T3RGBetaModel.wl`

**Problem solved:** both RGBeta runners need one common implementation of T3 field/representation construction.

**What it does:** validates supported dimensions, converts T3 representations to RGBeta representations and adds heavy fields/couplings.

## `T3RGBetaRunnerCommon.wl`

**Problem solved:** both runners need identical CLI parsing, JSON serialisation and beta-convention conversion.

**What it does:** provides integer parsing, JSON writing, expression serialisation, conventional report betas and LaTeX export.

---

# 11. `RGE/running/weinberg/` — final EFT analytic model and symbolic stages

## `WeinbergRGE.py`

**Problem solved:** below the last heavy threshold, the project needs the one-loop SM+Weinberg beta-function model.

**What it does:** implements the analytic $C_5$ beta matrix and associated SM beta functions needed by the numerical final-EFT runner.

For $C_5$:

$$
16\pi^2\frac{dC_5}{d\ln\mu}
=
(2\lambda_H-3g_2^2+2T)C_5
-\frac32\left[
Y_eY_e^\dagger C_5+C_5(Y_eY_e^\dagger)^T
\right].
$$

This file contains beta-function algebra only.

## `FullFlavorWeinbergStage.py`

**Problem solved:** the symbolic three-generation final-EFT beta matrix needs a pipeline output stage.

**What it does:** constructs full-flavor $C_5$, evaluates the symbolic beta matrix and writes the established stage result.

It does not convert $C_5$ to $m_\nu$.

## `OneGenerationWeinbergBenchmark.py`

**Problem solved:** the older one-generation tensor calculation remains useful as an analytic cross-check but should not be confused with the physical full-flavor stage.

**What it does:** computes the one-generation Weinberg RGE benchmark without file/output orchestration.

## `OneGenerationWeinbergBenchmarkStage.py`

**Problem solved:** the analytic benchmark needs a stable pipeline output contract.

**What it does:** loads inputs, calls the one-generation benchmark and writes benchmark output.

---

# 12. Overall RGE subsystem overview

```text
RGE/general/
    generic bases, generators and Wilson tensor algebra

RGE/group_factors/
    derive and independently validate representation-dependent coefficients

RGE/running/rgbeta/
    obtain renormalisable UV/intermediate beta functions from RGBeta

RGE/running/UVRunning.py
    production UV RGE stage

RGE/running/IntermediateEFTRunning.py
    choose the intermediate backend by physical field content

RGE/running/intermediate/
    scalar-only renormalisable and dimension-five evolution

RGE/matching/FinalWeinbergCoefficient.py
    hard threshold + direct running -> authoritative final C5

RGE/running/weinberg/WeinbergRGE.py
    final SM+Weinberg analytic beta model

Numerical/running/
    solve the corresponding ODEs numerically

physics/
    convert C5 into neutrino masses and observables
```

The main architectural distinction is:

> `RGE/` owns equations, symbolic tensor transport and matching bookkeeping; `Numerical/` owns benchmark substitution, state vectors and ODE integration.

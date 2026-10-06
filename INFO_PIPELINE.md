# T3 Neutrino-Mass Pipeline — Repository Map

This document is the repository-wide map of the implemented T3 radiative neutrino-mass calculation.

The purpose of this file is to answer:

1. **What problem is the pipeline solving?**
2. **Which file owns each part of that problem?**
3. **What is the execution order?**
4. **Where should a future change be made?**

Detailed file-by-file documentation is split into:

```text
INFO_LAGRANGIAN.md   UV model construction and Matchete matching
INFO_RGE.md          RGEs, threshold transport, group factors and final C5
INFO_NUMERICAL.md    numerical states, ODE running, fitting and plots
INFO_TESTS.md        Python/Wolfram tests and what each protects
```

---

## 1. Physics problem solved by the repository

The project implements a multi-threshold T3 radiative neutrino-mass calculation.

```text
T3 representation
    -> UV Lagrangian
    -> UV matching/RGEs
    -> heavy-fermion threshold
    -> scalar-only intermediate EFT
    -> scalar threshold
    -> final Weinberg coefficient C5
    -> SM + Weinberg running
    -> neutrino mass matrix
    -> masses / PMNS / oscillation observables
    -> numerical fit and figures
```

The project uses

$$
Q=T_3+Y,
$$

with

$$
Y(S_1)=\frac{\alpha}{2},\qquad
Y(S_2)=\frac{\alpha+2}{2},\qquad
Y(F)=\frac{\alpha+1}{2},
$$

and therefore

$$
Y_{\rm RZY}=2Y.
$$

The neutrino-mass convention is

$$
m_\nu=-\frac{v_{246}^2}{2}C_5=-v_{174}^2C_5.
$$

Current production supports the common threshold or the verified fermion-first hierarchy

```text
ordinary:      F -> (S1,S2)
shared scalar: F -> S
```

with EFT operators retained through dimension five.

Scalar-first execution is deliberately outside production scope because the first scalar threshold can generate a leading dimension-six operator schematically

$$
\frac{y\lambda_{T3}}{M_S^2}LFHHS.
$$

A scalar-first calculation truncated at $d\le5$ would therefore omit a leading EFT path.

---

## 2. Authoritative execution backbone

### `pipeline.py`

**Problem solved:** the project needs one visible place that states the physical calculation order without also containing all of the detailed algebra.

**What it does:**

- parses one run through `common/PipelineCLI.py`;
- obtains one `PipelinePlan`;
- selects requested T3 model points;
- calls the Lagrangian/matching boundary;
- attaches threshold/EFT metadata;
- runs UV RG evolution;
- dispatches intermediate EFT running by physical field content;
- runs independent intermediate-EFT validation;
- organises the matched Weinberg coefficient;
- runs low-energy Weinberg/neutrino stages;
- optionally invokes the integrated numerical benchmark pipeline;
- finishes with report generation and aggregate summaries.

**Boundary:** model-building algebra, beta functions, numerical ODE implementations, detailed matching formulae and report formatting live elsewhere.

---

# 3. Common pipeline data model

## `common/EFT.py`

**Problem solved:** threshold sequences need one generic description independent of historical labels such as `EFT1`.

**What it does:** defines `EFTContent`, `ThresholdStep`, `EFTTransition`, `EFTRunningInterval` and `EFTTruncation`.

The rest of the pipeline can ask "which fields are active?" rather than "which numbered EFT is this?".

## `common/PipelineCLI.py`

**Problem solved:** command-line parsing should not obscure the physical calculation in `pipeline.py`.

**What it does:** owns CLI options, ordinary/shared-scalar mode resolution, threshold-plan construction and study-directory naming.

## `common/PipelinePlan.py`

**Problem solved:** one run needs a single authoritative description of threshold order, scales and EFT truncation.

**What it does:** combines threshold steps/scales into the physical execution plan and constructs running intervals and stage metadata.

## `common/RunRecords.py`

**Problem solved:** stages need to pass paths, status and metadata without depending on one another's implementation.

**What it does:** defines `RunRecord` and `EFTStageRecord`. These are metadata containers, not physics calculators.

## `common/T3Fields.py`

**Problem solved:** the T3 topology uses formal `S1,S2` roles, but the shared-scalar/scotogenic branch contains one physical scalar.

**What it does:** provides the authoritative translation between

```text
ordinary physical fields: F, S1, S2
shared physical fields:   F, S
formal matching roles:     F, S1, S2
```

Threshold planning uses physical fields. Formal duplication is introduced only at the matching boundary.

## `common/T3Model.py`

**Problem solved:** representation-level checks must be possible before launching Matchete.

**What it does:** implements topology/dimension validity, neutral-component tests, class identification, shared-scalar formal dimensions and the neutral model scan.

## `common/Thresholds.py`

**Problem solved:** threshold ordering and scales must be consistent before any matching/running starts.

**What it does:** normalises threshold names, validates threshold groups, supplies defaults, resolves scales, converts physical plans to formal Wolfram roles and builds stage records.

---

# 4. Model/study selection

## `model/T3Study.py`

**Problem solved:** selecting which models to study should be separate from building their Lagrangians.

**What it does:** defines `T3ModelRequest`, ordinary scan definitions, neutrality filtering and explicit `--dims` selection.

It answers:

> Which UV model points should this study calculate?

## `studies/FullT3Study.py`

**Problem solved:** the honours project needs a reproducible 16-model numerical study plus controlled common-parameter comparisons.

**What it does:**

- enumerates the 16 ordinary neutral-compatible models;
- constructs per-model numerical configs;
- runs optimal benchmark mode;
- runs four common-parameter comparison scenarios;
- snapshots/restores generated model configs;
- builds aggregate interactive dashboards;
- writes `full_numerical_summary.json`.

This is study orchestration, not model physics.

---

# 5. Lagrangian and matching subsystem

Detailed file map: `INFO_LAGRANGIAN.md`.

The main Python route is

```text
model/T3Study.py
    -> Lagrangian/T3ModelMatching.py
    -> Lagrangian/Runner.py
```

with responsibilities split between

```text
Lagrangian/ModelValidation.py
    representation/request validation

Lagrangian/WolframRunner.py
    external wolframscript execution

Lagrangian/MatchingResults.py
    matching-output ingestion and RunRecord construction
```

Wolfram model building lives under

```text
Lagrangian/model/
Lagrangian/interactions/
```

with sequential execution controlled by

```text
Lagrangian/RunModel.wl
Lagrangian/RunMatching.wl
Lagrangian/RunThresholdStage.wl
```

---

# 6. RGE, intermediate EFT and final C5 subsystem

Detailed file map: `INFO_RGE.md`.

The main production path is

```text
RGE/running/UVRunning.py
    -> RGE/running/IntermediateEFTRunning.py
    -> RGE/running/backends/ScalarOnlyAfterFermion.py
    -> RGE/matching/FinalWeinbergCoefficient.py
    -> RGE/running/weinberg/WeinbergRGE.py
```

The verified intermediate theory after removing the fermion is

```text
ordinary:      SM + S1 + S2 + dimension-five Wilson operators
shared scalar: SM + S       + dimension-five Wilson operators
```

The authoritative final coefficient is matching/bookkeeping, not another RGE integrator:

```text
RGE/matching/FinalWeinbergCoefficient.py
```

and is adapted into a symbolic physical Majorana matrix by

```text
RGE/matching/FinalWeinbergAdapter.py
```

---

# 7. Low-energy neutrino physics

## `physics/LowEnergyNeutrino.py`

**Problem solved:** post-matching neutrino stages need one high-level facade while the actual physics remains in dedicated modules.

**What it does:** orchestrates the final-C5 location, one-generation benchmark, full-flavor Weinberg stage, symbolic mass construction and numerical neutrino-observable output.

## `physics/NeutrinoMass.py`

**Problem solved:** the $C_5\to m_\nu$ convention must have one owner.

**What it does:** implements

$$
m_\nu=-\frac{v^2}{2}C_5
$$

for symbolic and numerical use.

## `physics/NeutrinoObservables.py`

**Problem solved:** a complex symmetric Majorana mass matrix must be converted into physical masses, mass splittings and PMNS information.

**What it does:** performs Takagi factorisation and constructs physical neutrino observables.

## `physics/NeutrinoTrajectory.py`

**Problem solved:** a running $C_5(\mu)$ trajectory needs physical interpretation at every scale.

**What it does:**

```text
C5(mu)
    -> m_nu(mu)
    -> charged-lepton mass basis
    -> scale-dependent observables
```

without performing ODE integration itself.

---

# 8. Numerical subsystem

Detailed file map: `INFO_NUMERICAL.md`.

```text
Numerical/
├── core/
├── running/
├── fitting/
├── diagnostics/
├── plotting/
└── orchestration/
```

Ownership:

```text
core/           state representation and beta evaluation
running/        numerical EFT evolution and threshold boundaries
fitting/        scans, chi-square, optimisation and benchmark search
diagnostics/    retained physics information/checks
plotting/       presentation only
orchestration/  config persistence and integrated pipeline execution
```

---

# 9. Independent validation

## `validation/IntermediateEFTValidation.py`

**Problem solved:** production running and scientific validation should not be the same code path.

**What it does:** dispatches independent checks for a completed intermediate-EFT interval by physical field content.

## `validation/backends/ScalarOnlyAfterFermionValidation.py`

**Problem solved:** the scalar-only production backend needs independent checks of component transport and pole/RGE consistency.

**What it does:** owns those checks for the verified fermion-first scalar-only interval.

## `validation/__init__.py`, `validation/backends/__init__.py`

Package markers only.

---

# 10. Reports

## `Reports/PipelineReports.py`

**Problem solved:** report construction and aggregate output should not clutter `pipeline.py`.

**What it does:** prints run summaries, orchestrates human-readable reports and preserves aggregate output schemas.

## `Reports/ReportGeneration.py`

**Problem solved:** all reports need shared LaTeX/paper notation and path rules.

**What it does:** owns report preamble, notation conversion and output-path helpers.

## `Reports/GroupFactorReports.py`

**Problem solved:** representation/group-factor information needs compact cross-model reporting.

**What it does:** constructs group-factor/RGE comparison tables from successful run records.

## `Reports/RGEComparison.py`

**Problem solved:** stored RGBeta/intermediate expressions are machine-oriented and need consistent report conversion.

**What it does:** loads RGE payloads and rewrites symbolic structures into report-ready expressions.

---

# 11. Helper scripts

## `scripts/run_comparison_study.py`

**Problem solved:** debugging one common-parameter model/scenario should not require running all 16 models.

**What it does:** builds one comparison config and launches the normal `pipeline.py` path.

## `scripts/run_regression.py`

**Problem solved:** the external smoke + Wolfram C5 validation needs a reproducible regression command.

**What it does:** runs smoke matching and then `tests/wolfram/RegressionC5.wl`, preventing stale reports from masking a failed Wolfram run.

---

# 12. Tests

Detailed file map: `INFO_TESTS.md`.

The suite has four roles:

```text
unit/numerical tests
    verify formulas, state packing and deterministic scan behaviour

architecture tests
    verify ownership boundaries and retired historical modules

integration-style tests
    verify multi-stage trajectories and output contracts

Wolfram/reference tests
    verify group theory, normalization, RGBeta and matched C5 conventions
```

The Python suite does not replace a real external Wolfram/Matchete/RGBeta run.

---

# 13. Where to make a change

| Desired change | Primary owner |
|---|---|
| CLI option | `common/PipelineCLI.py` |
| threshold plan | `common/Thresholds.py`, `common/PipelinePlan.py` |
| shared scalar mapping | `common/T3Fields.py` |
| representation validity | `common/T3Model.py` |
| model selection | `model/T3Study.py` |
| UV Lagrangian/invariant | `Lagrangian/model/`, `Lagrangian/interactions/` |
| Matchete process execution | `Lagrangian/WolframRunner.py` |
| UV RGBeta model/export | `RGE/running/rgbeta/` |
| intermediate EFT equations | `RGE/running/intermediate/` |
| final C5 construction | `RGE/matching/FinalWeinbergCoefficient.py` |
| final SM+Weinberg beta | `RGE/running/weinberg/WeinbergRGE.py` |
| $C_5\to m_\nu$ | `physics/NeutrinoMass.py` |
| neutrino observables | `physics/NeutrinoObservables.py` |
| numerical integration | `Numerical/running/` |
| fit/scan | `Numerical/fitting/` |
| plots only | `Numerical/plotting/` |
| 16-model study | `studies/FullT3Study.py` |
| reports | `Reports/` |

---

# 14. Overall pipeline overview

Each layer answers one question:

```text
common/
    What theory and threshold plan are we asking for?

model/
    Which T3 model points should be calculated?

Lagrangian/
    What is the UV theory and what EFT coefficients are generated?

RGE/
    How do coefficients run and what is the authoritative final C5?

physics/
    What physical neutrino masses and observables follow from C5?

Numerical/
    What does that calculation give numerically for a benchmark or scan?

validation/
    Do independent checks support the production result?

Reports/ and Numerical/plotting/
    How are stored results presented without changing the physics?

tests/
    Which physics, numerical and architecture contracts must never regress?
```

`pipeline.py` is deliberately the map connecting those questions, not the place where they are individually solved.

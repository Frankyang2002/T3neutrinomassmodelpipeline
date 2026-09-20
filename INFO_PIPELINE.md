# T3 Neutrino-Mass Model Pipeline: End-to-End Guide

This document is the repository-wide overview for the T3 neutrino-mass model pipeline.

It is intentionally broader and shorter than the two detailed implementation notes:

- [`INFO_LAGRANGIAN.md`](INFO_LAGRANGIAN.md): UV model construction, interaction generation, Matchete matching, Weinberg-operator extraction, and Lagrangian conventions.
- [`INFO_RGE.md`](INFO_RGE.md): UV/intermediate-EFT running, Weinberg-coefficient RGEs, numerical running, neutrino masses, and observables.

The purpose of this file is to answer three questions:

1. What does the complete pipeline do?
2. How does information move from one stage to the next?
3. Where does each source file sit in that pipeline?

For detailed equations and implementation notes, the two specialised INFO files remain the primary references.

---

## 1. Big Picture

The project takes a T3 radiative neutrino-mass model from a representation assignment to a low-energy neutrino prediction.

Schematically,

$$
\text{T3 representations}
\longrightarrow
\mathcal L_{\rm UV}
\longrightarrow
\text{EFT matching}
\longrightarrow
C_5
\longrightarrow
\text{RGE running}
\longrightarrow
m_\nu
\longrightarrow
\text{neutrino observables}.
$$

With sequential thresholds, the middle of the pipeline is more accurately

$$
\mathcal L_{\rm UV}
\longrightarrow
{\rm EFT}_1
\longrightarrow
{\rm EFT}_2
\longrightarrow \cdots \longrightarrow
{\rm SMEFT}
\longrightarrow C_5(\mu).
$$

The project therefore has several physically distinct layers:

```text
model definition
    ↓
UV Lagrangian construction
    ↓
one-loop matching
    ↓
threshold bookkeeping / intermediate EFTs
    ↓
UV and EFT RGE running
    ↓
final Weinberg coefficient
    ↓
SMEFT Weinberg running
    ↓
Majorana neutrino mass matrix
    ↓
neutrino observables
    ↓
reports and regression checks
```

---

## 2. Main Physics Conventions

### 2.1 Hypercharge convention

The Lagrangian side uses

$$
Q=T_3+Y.
$$

For the T3 parameter $\alpha$,

$$
Y(S_1)=\frac{\alpha}{2},
\qquad
Y(S_2)=\frac{\alpha+2}{2},
\qquad
Y(F)=\frac{\alpha+1}{2}.
$$

This differs by a factor of two from papers using $Q=T_3+Y/2$.

### 2.2 Historical T3 classes

The historical classes are

| Class | $d_{S_1}$ | $d_{S_2}$ | $d_F$ |
|---|---:|---:|---:|
| T3-A | 1 | 3 | 2 |
| T3-B | 2 | 2 | 1 |
| T3-C | 2 | 2 | 3 |
| T3-D | 3 | 1 | 2 |
| T3-E | 3 | 3 | 2 |

The dimensions are SU(2) representation dimensions, not isospins.

The generic T3 representation conditions are

$$
d_{S_1}=d_F\pm1,
\qquad
d_{S_2}=d_F\pm1,
$$

together with the requirement that the scalar product contain the triplet channel needed by the Higgs pair.

### 2.3 Shared-scalar mode

The ordinary T3 bookkeeping contains the formal roles

$$
S_1,\quad S_2,\quad F.
$$

The shared-scalar branch instead contains one physical scalar $S$ and the fermion $F$. Internally the formal topology may still use $S_1$ and $S_2$, but both roles belong to the same physical scalar.

This distinction is important for thresholds and reporting:

```text
ordinary T3 physical fields:
F, S1, S2

shared-scalar physical fields:
F, S
```

while the Wolfram matching layer can still use the formal T3 roles.

### 2.4 Matching order

The normal matching target is

$$
d_{\rm EFT}=5,
\qquad
L=1,
$$

i.e. dimension-five operators at one loop.

### 2.5 Weinberg and neutrino-mass conventions

The low-energy EFT contains

$$
\mathcal L_{\rm EFT}
\supset
C_5\,\mathcal O_5+\text{h.c.}
$$

and the pipeline uses the convention bridge

$$
m_\nu=-v_{174}^2 C_5
$$

or equivalently

$$
m_\nu=-\frac{v_{246}^2}{2}C_5,
$$

with $v_{246}/\sqrt2=v_{174}$.

The RGE convention uses

$$
t=\ln\mu.
$$

---

## 3. End-to-End Execution

### Stage 1: command-line orchestration

`pipeline.py` reads the requested model, validates the representation assignment, prepares the output directories, chooses the threshold plan, and calls the relevant Lagrangian/RGE stages.

### Stage 2: model definition

The Python-side representation bookkeeping in `common/T3Model.py` and `Lagrangian/Runner.py` is translated into a Wolfram-side model association by `Lagrangian/model/T3ModelCatalog.wl`.

At this point the code knows the SU(2) dimensions, hypercharges, field names, and whether the requested representation assignment is valid.

### Stage 3: UV Lagrangian construction

`Lagrangian/model/LagrangianBuilder.wl` builds the BSM extension of the Standard Model.

The topology-generating interactions are schematically

$$
y_1 LFS_1
+
y_2 LFS_2
+
\lambda_{T3}HHS_1S_2^\dagger
+\text{h.c.}
$$

with the actual conjugations and SU(2) contractions determined from the representations.

Supporting files define:

- fields and couplings;
- SU(2) invariant tensors / Clebsch-Gordan structures;
- scalar-potential terms;
- T3 topology interactions.

The resulting candidate interactions are checked before the full UV Lagrangian is accepted.

### Stage 4: one-loop EFT matching

`Lagrangian/RunMatching.wl` sends the UV theory through the Matchete matching chain,

```text
Match
  ↓
GreensSimplify
  ↓
EOMSimplify
  ↓
EvaluateLoopFunctions
  ↓
ReplaceEffectiveCouplings
```

and also performs a Standard-Model-only match.

The BSM matched EFT is obtained by subtraction,

$$
\mathcal L_{\rm EFT}^{\rm BSM}
=
\mathcal L_{\rm EFT}^{\rm full}
-
\mathcal L_{\rm EFT}^{\rm SM}.
$$

The dimension-five result is then searched for the Weinberg operator and its coefficient $C_5$ is extracted.

### Stage 5: threshold handling

`common/Thresholds.py` defines which physical heavy fields are removed at each threshold.

A threshold plan such as

```text
F -> (S1,S2)
```

means that the heavy fermion is integrated out first, producing an intermediate EFT, and the two scalars are removed at a later common threshold.

For when we have the same scalars, the physical plan can be

```text
F -> S
```

while the Wolfram matcher receives the corresponding formal T3 scalar roles.

`Lagrangian/RunThresholdStage.wl` performs the later threshold matching calculations.

### Stage 6: UV and intermediate-EFT running

The UV T3 theory and supported intermediate EFTs can be run using RGBeta through `RGE/running/rgbeta/` files.

When a fermion threshold leaves a scalar-containing intermediate EFT, the `RGE/running/eft1/` modules construct and run the required Wilson tensors between thresholds.

The threshold-running contribution is kept distinct from the hard matching contribution until the final Weinberg coefficient is assembled.

### Stage 7: final Weinberg coefficient

`RGE/running/weinberg/FinalWeinbergCoefficient.py` combines the pieces that contribute to the final physical Weinberg coefficient after the heavy states have been removed.

This is the point at which the threshold history has been converted into the coefficient that should be evolved in the final SMEFT.

### Stage 8: symbolic Weinberg RGE

`RGE/matching/MatchedWeinbergRGE.py` parses the matched coefficient and evaluates the generic $\psi^2\phi^2$ master RGE.

The one-generation check is

$$
\frac{16\pi^2\beta_{C_5}}{C_5}
=
-3g_2^2
+2\lambda_H
+6|y_u|^2
+6|y_d|^2
-|y_e|^2.
$$

The full three-generation SMEFT equation used later is

$$
16\pi^2\frac{dC_5}{d\ln\mu}
=
(2\lambda_H-3g_2^2+2T)C_5
-\frac32
\left[
Y_eY_e^\dagger C_5
+
C_5(Y_eY_e^\dagger)^T
\right],
$$

with

$$
T=
\operatorname{Tr}
\left(
Y_eY_e^\dagger
+3Y_uY_u^\dagger
+3Y_dY_d^\dagger
\right).
$$

### Stage 9: flavor lift

The one-generation matched expression is promoted to a complex symmetric three-generation matrix.

Schematically,

$$
(C_5)_{pq}
=
\frac12\sum_r F_r
\left[
y_{1,pr}^*y_{2,qr}^*
+
y_{2,pr}^*y_{1,qr}^*
\right].
$$

The flavor construction currently assumes a diagonal heavy-fermion mass basis in the path documented by `INFO_RGE.md`.

### Stage 10: numerical running

When `--numerical CONFIG.json` is supplied, the pipeline numerically evolves:

- $g_Y,g_2,g_3$;
- $\lambda_H$;
- diagonal SM Yukawa couplings;
- the complex symmetric $C_5$ matrix.

The independent variable is $\ln\mu$, and the current implementation uses SciPy `solve_ivp` with the `DOP853` method.

### Stage 11: neutrino observables

At the low scale,

$$
m_\nu=-\frac{v^2}{2}C_5
$$

in the $v\simeq246$ GeV convention.

`RGE/phenomenology/NeutrinoObservables.py` Takagi-factorises the complex symmetric Majorana mass matrix,

$$
U^T m_\nu U
=
\operatorname{diag}(m_1,m_2,m_3),
$$

and extracts masses, mass-squared splittings, and mixing information.

### Stage 12: reports and regression checks

The `Reports/` modules collect readable summaries of the model, matching, RGEs, group factors, and comparisons.

`regression.py` and the files under `tests/` are independent validation layers. They should not be treated as production stages even where they reproduce part of the production calculation.

---

## 4. Repository File Map

The descriptions below are intentionally short. The detailed Lagrangian and RGE documents contain the equations and implementation-specific discussion.

### `common/`

#### `common/T3Model.py`
Shared Python model bookkeeping.

Defines the historical T3 A-E representation assignments, supported dimension checks, benchmark sets, alpha-name encoding, and the formal-versus-physical mapping used by shared-scalar mode.

**Role:** model-definition foundation.  


#### `common/RunRecords.py`
Defines `RunRecord` and `EFTStageRecord`, the metadata objects used to carry model/run and threshold-stage information through the Python pipeline.

**Role:** pipeline bookkeeping.  


#### `common/Thresholds.py`
Defines physical heavy fields, threshold plans, threshold-plan validation, human/JSON labels, the conversion from physical shared-scalar thresholds to formal T3 matching roles, and the construction of EFT-stage records.

**Role:** threshold bookkeeping.  


---

## 5. Lagrangian and Matching Files

### `Lagrangian/model/`

#### `Lagrangian/model/T3ModelCatalog.wl`
Defines the Wolfram-side T3 model data: historical classes, generic dimensions, alpha-dependent hypercharges, and representation-level model construction.

**Role:** model definition.  


#### `Lagrangian/model/T3Fields.wl`
Defines the T3 BSM fields, couplings, and their index structures for use by the Lagrangian builder.

**Role:** field/coupling definition.  


#### `Lagrangian/model/SU2Invariants.wl`
Constructs the SU(2) invariant tensors / Clebsch-Gordan structures required for the T3 interactions, using Matchete invariant tensors.

**Role:** SU(2) representation algebra.  


#### `Lagrangian/model/LagrangianBuilder.wl`
Assembles the complete T3 UV Lagrangian from the model data, field definitions, invariant tensors, and interaction modules. Candidate terms are checked before the UV theory is accepted.

**Role:** UV model construction.  


### `Lagrangian/interactions/`

#### `Lagrangian/interactions/T3Topology.wl`
Defines the interactions specifically required to realise the T3 one-loop neutrino-mass topology.

**Role:** topology-generating interactions.  


#### `Lagrangian/interactions/ScalarPotential.wl`
Builds the gauge-invariant scalar interactions needed by the two-scalar T3 UV theory, including the scalar structures used by matching and RGE stages.

**Role:** scalar potential.  


#### `Lagrangian/interactions/ExportUVCGRegistry.wl`
Exports the Clebsch-Gordan/invariant-tensor information needed by later Python RGE and threshold stages.

**Role:** bridge from Wolfram model construction to downstream RGE tensor code.  


### `Lagrangian/`

#### `Lagrangian/Runner.py`
Python wrapper around the Wolfram model runner. It validates model dimensions, converts the physical threshold plan to the formal plan needed by Wolfram, launches `RunModel.wl`, captures output, and returns a `RunRecord`.

**Role:** Python/Wolfram orchestration boundary.  


#### `Lagrangian/RunModel.wl`
Main Wolfram-side orchestration script. It connects model construction, validation, matching, Weinberg extraction, output conversion, and summary generation.

**Role:** Lagrangian/matching orchestrator.  


#### `Lagrangian/RunMatching.wl`
Runs the one-loop Matchete matching, simplifies the EFT, subtracts the SM baseline, identifies the Weinberg operator, and extracts $C_5$.

**Role:** EFT matching.  


#### `Lagrangian/RunThresholdStage.wl`
Runs matching at later thresholds when the heavy fields are integrated out sequentially rather than at one common scale.

**Role:** sequential threshold matching.  


#### `Lagrangian/PhysicsNotation.wl`
Converts Matchete's internal field/index/coupling representation into readable physics notation and LaTeX.

**Role:** presentation/serialization support.  
**Production:** yes, but it does not define the physics.

---

## 6. General RGE Infrastructure

### `RGE/general/`

#### `RGE/general/ScalarBasis.py`
Defines the scalar representation/basis model used by the generic tensor RGE machinery.

**Role:** generic RGE model/basis definition.  


#### `RGE/general/GaugeGenerators.py`
Builds the real-scalar SU(2) and U(1) gauge generators and gauge-sector representation data.

**Role:** gauge group representation algebra.  


#### `RGE/general/FermionBasis.py`
Constructs the Weyl-fermion basis and the corresponding fermion gauge generators.

**Role:** fermion basis for generic RGEs.  


#### `RGE/general/WilsonTensorRGE.py`
Implements the general tensor terms entering the one-loop $\psi^2\phi^2$ Wilson-coefficient RGE.

**Role:** master Wilson-tensor RGE.  


#### `RGE/general/AnomalousDimensions.py`
Builds scalar and fermion anomalous-dimension contributions and combines them with the Wilson-tensor terms into the complete generic RGE.

**Role:** master anomalous-dimension assembly.  


---

## 7. Group-Factor Infrastructure

### `RGE/group_factors/core/`

#### `RepresentationFactors.py`
Provides basic representation-theory factors needed by the RGE group-factor calculations.

#### `YukawaGroupFactors.py`
Computes group-theory factors associated with Yukawa structures.

#### `ScalarQuarticBasis.py`
Defines/constructs the scalar-quartic tensor basis used by the production RGE machinery.

#### `QuarticTensorAlgebra.py`
Performs tensor algebra for scalar quartic structures.

#### `MixingQuarticTensorAlgebra.py`
Handles quartic tensor structures that mix different scalar sectors/channels.

These are **production group-factor utilities**.

### `RGE/group_factors/recoupling/`

#### `PortalQuarticGroupFactors.py`
Constructs the group factors associated with portal-type quartic structures.

#### `PortalRecouplingDerivations.py`
Carries out the recoupling derivations needed to express portal quartics in the basis used by the RGE pipeline.

#### `MixingQuarticRecoupling.py`
Handles recoupling for scalar mixing quartics.

#### `RGBetaMixingQuarticRecoupling.py`
Connects the internal mixing-quartic basis to RGBeta conventions/results.

These files are part of the **production recoupling layer**, although some functions also support validation.

### `RGE/group_factors/validation/`

#### `ValidateDirectWeinbergGroupFactors.py`
Checks group factors entering direct Weinberg-operator contributions.

#### `ValidateMassGroupFactors.py`
Checks group factors associated with mass-dependent structures.

#### `ValidateMixingQuarticGroupFactors.py`
Checks the scalar mixing-quartic group factors.

#### `ValidateRGBetaRecouplings.py`
Cross-checks the internal recoupling conventions against RGBeta outputs.

#### `YukawaBetaGroupFactors.py`
Provides/validates group factors associated with Yukawa beta functions.

These files are **validation/research cross-check code**, not the primary high-level pipeline entry point.

---

## 8. Matching-to-RGE Bridge

### `RGE/matching/`

#### `RGE/matching/MatchedWeinbergRGE.py`
Parses the matched $C_5$, embeds it in the generic Wilson-tensor machinery, evaluates its one-loop RGE, and checks the one-generation SMEFT benchmark.

**Role:** matched-coefficient RGE validation.  


#### `RGE/matching/FlavorC5Matching.py`
Promotes the one-generation matched result to the symmetric three-generation flavor coefficient matrix.

**Role:** flavor reconstruction.  


#### `RGE/matching/WeinbergTensorAdapter.py`
Maps the Weinberg coefficient into the general $\psi^2\phi^2$ tensor representation used by the master RGE code.

**Role:** adapter between model-specific $C_5$ and generic tensor RGEs.  


---

## 9. UV and Intermediate-EFT Running

### `RGE/running/rgbeta/`

#### `RGBetaT3Running.py`
Python interface that launches and reads the UV/EFT1 RGBeta calculations.

#### `RunT3RGBeta.wl`
Runs RGBeta for the full UV T3 theory.

#### `RunT3EFT1RGBeta.wl`
Runs RGBeta for the supported first intermediate EFT after the first threshold.

#### `T3RGBetaModel.wl`
Builds the RGBeta model definition from the T3 data exported by the Lagrangian side.

#### `T3RGBetaRunnerCommon.wl`
Shared runner utilities used by both UV and EFT1 RGBeta scripts.

These files provide the **UV/intermediate-EFT beta functions**. The documented current RGBeta path is limited by the representations supported by the present model setup.

### `RGE/running/eft1/`

#### `EFT1TensorAdapters.py`
Converts the first intermediate EFT's Matchete exports into the Wilson and quartic tensors expected by the generic RGE machinery.

#### `EFT1WilsonRGE.py`
Evaluates the intermediate-EFT Wilson-coefficient RGE.

#### `EFT1WilsonFlow.py`
Transports the EFT1 Wilson coefficients between thresholds and exports the flavor seed used by the next threshold calculation.

#### `EFT1DirectWeinberg.py`
Builds the direct one-loop $C_{12}\to C_5$ running contribution.

#### `EFT1ThresholdResume.py`
Restarts the later threshold matching calculation using the Wilson coefficients transported through EFT1.

#### `MatcheteParsing.py`
Parses Matchete matching expressions and CG registries needed by the EFT1 tensor construction.

These files are **production intermediate-EFT running code**.

---

## 10. Final Weinberg Running

### `RGE/running/weinberg/`

#### `FinalWeinbergCoefficient.py`
Combines hard threshold matching and threshold-running contributions into the final Weinberg coefficient.

#### `FlavorMatchedWeinbergStage.py`
Runs the symbolic full-flavor Weinberg beta matrix and constructs the symbolic neutrino-mass matrix.

#### `WeinbergRunning.py`
Contains the three-generation symbolic SMEFT Weinberg RGE and the numerical coupled SM+$C_5$ evolution.

#### `NumericalWeinbergStage.py`
Reads/evaluates the matched coefficient for a numerical configuration and drives the numerical low-energy RGE stage.

These files are the **final SMEFT running layer**.

---

## 11. Phenomenology

### `RGE/phenomenology/NeutrinoObservables.py`

Takagi-factorises the low-scale complex symmetric Majorana mass matrix and extracts the physical neutrino masses, mass-squared splittings, and mixing information.

**Role:** phenomenology / observable extraction.  


---

## 12. Reports

### `Reports/ReportGeneration.py`
Generates the main human-readable model/matching/RGE report products and the LaTeX/PDF outputs used by the pipeline.

### `Reports/RGEComparison.py`
Builds comparisons between the different RGE calculations/stages.

### `Reports/GroupFactorReports.py`
Produces reports focused on the group factors and recoupling structures.

The `Reports/` directory is **production presentation code**: it should consume physics outputs rather than define the physics itself.

---

## 13. Root Orchestration

### `pipeline.py`

Main user-facing entry point.

It connects:

```text
model selection
→ UV construction
→ matching
→ threshold handling
→ UV/EFT running
→ final C5
→ symbolic/numerical Weinberg running
→ neutrino observables
→ reports
```

It also implements the multi-model comparison studies and manages the output directory structure.

### `regression.py`

Repository-level regression gate.

It reruns the historical benchmark models and validates important matching outputs, including the scotogenic/T3-B limit.

**Role:** validation, not physics production.

---

## 14. Tests

The tests deliberately duplicate or independently reproduce some calculations. This is useful: a cross-check should not simply call the same production function it is supposed to validate.

### `tests/python/`

#### `check_t3_scalar_quartics.py`
Independent Python check of T3 scalar-quartic structures/group factors.

#### `test_rgbeta_t3_running.py`
Regression tests for the Python RGBeta T3-running interface.

#### `test_shared_scalar_mode.py`
Checks the physical/formal shared-scalar mapping, threshold-plan expansion, and scalar-block counting.

#### `test_t3_scotogenic_normalization.py`
Checks the T3/scotogenic matching and normalisation conventions.

#### `test_weinberg_normalization_cleanup.py`
Checks Weinberg-operator normalisation and conversion conventions.

### `tests/wolfram/`

#### `ProbeT3SU2Factors.wl`
Independent probe of SU(2) factors.

#### `RegressionC5.wl`
Wolfram regression check of the matched Weinberg coefficient.

#### `T3CouplingConventions.wl`
Reference definitions for T3 coupling conventions used by Wolfram-side tests.

#### `TestT3Conventions.wl`
Checks the Wolfram T3 convention definitions.

#### `TestT3RGBeta.wl`
Checks RGBeta-related Wolfram calculations.

#### `TestT3RGEExport.wl`
Checks export of the tensors/data required by the RGE pipeline.

#### `TestWeinberg.wl`
Checks Weinberg-operator extraction and related conventions.

### `tests/reference/`

#### `MaUVRGE.py`
Independent Ma/scotogenic UV-RGE reference implementation.

It is intentionally **not** a production module.

#### `T3RGEComponentExport.wl`
Independent export/check of T3 RGE component structures.

#### `T3RGEPotentialExport.wl`
Independent export/check of scalar-potential RGE data.

#### `T3RGETensorExport.wl`
Independent export/check of full T3 RGE tensor structures.

The `tests/reference/` files are best understood as **reference/provenance and validation calculations**, not code that the main pipeline should depend on.


## 15. Typical Output Flow

A normal model run creates a directory under

```text
output/<model>/
```

The main conceptual products are

```text
UV model / matching summaries
        ↓
matched Weinberg coefficient
data/c5_coefficient.txt
        ↓
symbolic RGE results
data/c5_beta.txt
data/c5_flavor_beta_matrix.txt
        ↓
symbolic neutrino mass
data/neutrino_mass_matrix.txt
        ↓
optional numerical evolution
data/c5_flavor_matrix_low_scale.txt
data/neutrino_mass_matrix_low_scale.txt
        ↓
phenomenology
data/neutrino_observables.json
        ↓
human-readable reports
c5_coefficient.pdf
rge_report.pdf
```

With `--debug-reports`, additional intermediate algebra and Wolfram logs are written under `debug/`.

---

## 16. Useful Commands

Run the five historical smoke models:

```bash
python pipeline.py --smoke
```

Run an ordinary T3 representation assignment:

```bash
python pipeline.py --dims 2 2 1 --alpha -1
```

Run the shared-scalar branch:

```bash
python pipeline.py --dims 2 1
```

Run with sequential thresholds:

```bash
python pipeline.py --smoke --threshold F --threshold S1 S2
```

Run with detailed diagnostics:

```bash
python pipeline.py --smoke --debug-reports
```

Run the repository regression gate:

```bash
python regression.py
```

Run it without redoing the expensive matching:

```bash
python regression.py --no-rematch
```

Run the Python tests:

```bash
python -m pytest -q tests/python
```
---


For understanding the code rather than simply running it, a useful dependency order is:

```text
1. common/
2. Lagrangian/model/
3. Lagrangian/interactions/
4. Lagrangian/Runner.py + RunModel.wl
5. RunMatching.wl + RunThresholdStage.wl
6. RGE/general/
7. RGE/group_factors/core/
8. RGE/group_factors/recoupling/
9. RGE/group_factors/validation/
10. RGE/matching/
11. RGE/running/eft1/
12. RGE/running/rgbeta/
13. RGE/running/weinberg/
14. RGE/phenomenology/
15. Reports/
16. pipeline.py
17. regression.py
18. tests/python/
19. tests/wolfram/
20. tests/reference/
21. INFO_LAGRANGIAN.md + INFO_RGE.md
```

The most important distinction while reading is whether a file is:

- defining the physical model;
- doing matching;
- transporting information across a threshold;
- computing a beta function;
- converting a coefficient into an observable;
- reporting a result;
- or independently validating the calculation.

That distinction is more useful than the programming language of the file.

---

## 17. Short Summary

The complete project can be thought of as three connected calculations:

### A. Build and match the radiative model

$$
(S_1,S_2,F)
\rightarrow
\mathcal L_{\rm UV}
\rightarrow
\mathcal L_{\rm EFT}
\rightarrow
C_5.
$$

### B. Run through the relevant scales

$$
{\rm UV}
\rightarrow
{\rm intermediate\ EFTs}
\rightarrow
C_5(M)
\rightarrow
C_5(\mu).
$$

### C. Convert to neutrino physics

$$
C_5(\mu)
\rightarrow
m_\nu(\mu)
\rightarrow
m_i,\Delta m^2,U_{\rm PMNS}.
$$

`pipeline.py` is the orchestration layer connecting all three, while the files under `tests/` and `RGE/group_factors/validation/` provide independent checks that the calculation is being performed consistently.

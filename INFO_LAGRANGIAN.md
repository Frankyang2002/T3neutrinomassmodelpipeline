# T3 UV Lagrangian and Matching — File-by-File Guide

This document explains the Lagrangian/model-building side of the project file by file.

For every file the emphasis is:

- **Problem solved:** why the file exists.
- **What it does:** the concrete implementation responsibility.
- **Boundary:** what should remain outside that file.

The subsystem starts from a selected T3 representation and ends with matched EFT artifacts plus a Python `RunRecord`.

---

# 1. Physics conventions

The project uses

$$
Q=T_3+Y.
$$

For the Restrepo-Zapata-Yaguna parameter $\alpha$,

$$
Y(S_1)=\frac{\alpha}{2},\qquad
Y(S_2)=\frac{\alpha+2}{2},\qquad
Y(F)=\frac{\alpha+1}{2},
$$

and therefore

$$
Y_{\rm RZY}=2Y.
$$

The schematic T3 interaction content is

$$
\mathcal L_{\rm T3}
\supset
y_1LFS_1+y_2LFS_2
+\lambda_{T3}HHS_1S_2^\dagger+\mathrm{h.c.},
$$

with representation-dependent conjugations and Clebsch-Gordan tensors supplied by the Wolfram builder.

Ordinary models have physical fields `F,S1,S2`.

The shared-scalar/scotogenic branch has physical fields `F,S`, where the same scalar fills the formal topology roles `S1,S2` up to conjugation. `common/T3Fields.py` is the Python source of truth for this translation.

---

# 2. Python boundary files

## `Lagrangian/T3ModelMatching.py`

**Problem solved:** the central pipeline needs a clean Python entry point for "build this selected T3 model and perform the configured matching", without knowing how Matchete is launched.

**What it does:**

- accepts an already-selected `T3ModelRequest`;
- converts the `PipelinePlan` into the matching call expected by the Lagrangian stack;
- calls the stable `Lagrangian.Runner` interface;
- returns the resulting `RunRecord`;
- handles lists of selected models for one study.

**Boundary:** it does not reimplement any matching formula and does not decide which model points belong to the study.

## `Lagrangian/Runner.py`

**Problem solved:** older code expects stable Python entry points such as `run_model`, while the implementation responsibilities are now separated.

**What it does:** provides the stable public facade and delegates to:

```text
ModelValidation.py
WolframRunner.py
MatchingResults.py
```

It also retains compatibility helpers for dimension validation and class-dimension lookup.

**Boundary:** it should remain thin. Detailed subprocess handling and JSON interpretation belong elsewhere.

## `Lagrangian/ModelValidation.py`

**Problem solved:** invalid/unsupported T3 requests must be rejected before starting an expensive external Wolfram process.

**What it does:**

- validates ordinary and shared-scalar dimensions;
- constructs `ModelRunSpecification`;
- derives stable model/output names;
- maps T3 class labels to dimensions.

**Boundary:** pure Python validation/configuration only. It does not launch Wolfram or read output files.

## `Lagrangian/WolframRunner.py`

**Problem solved:** `wolframscript` execution has process-management concerns that should not be mixed with model physics or output interpretation.

**What it does:**

- prepares the model output directory;
- converts physical threshold steps to formal Wolfram roles;
- constructs the external command;
- streams stdout;
- handles interruption/failure;
- writes process/debug logs.

**Boundary:** it does not decide whether a T3 representation is valid and does not interpret the matching result.

## `Lagrangian/MatchingResults.py`

**Problem solved:** Wolfram outputs need to be converted into stable Python runtime metadata.

**What it does:**

- loads `comparison_summary.json`;
- converts formal shared-scalar metadata back to the physical scalar language;
- verifies that sequential stage export is consistent with the requested threshold plan;
- constructs the final `RunRecord`.

**Boundary:** it does not launch Wolfram and does not validate the original requested representation.

---

# 3. Top-level Wolfram execution files

## `Lagrangian/RunModel.wl`

**Problem solved:** one executable Wolfram entry point must turn command-line input into a complete model/matching run.

**What it does:**

- parses signed integer and threshold-plan CLI arguments;
- loads the model-building/matching stack;
- constructs the requested T3 model;
- executes the requested matching sequence;
- exports Weinberg/matching artifacts and summary metadata.

This is the main Wolfram-side orchestration file.

## `Lagrangian/RunMatching.wl`

**Problem solved:** matching itself needs a reusable staged implementation separate from CLI parsing.

**What it does:**

- defines safe stage execution;
- performs the requested Matchete matching stages;
- runs the simplify/evaluate sequence;
- manages the sequence used by `RunModel.wl`.

Conceptually:

```text
UV Lagrangian
    -> Match
    -> GreensSimplify
    -> EOMSimplify
    -> EvaluateLoopFunctions
    -> ReplaceEffectiveCouplings
    -> EFT result
```

## `Lagrangian/RunThresholdStage.wl`

**Problem solved:** sequential threshold calculations need to resume from a previous EFT and integrate out the next physical threshold group.

**What it does:** implements the lower-threshold continuation/resume stage used by sequential matching.

This is particularly important for the verified fermion-first hierarchy, where the scalar threshold is processed after the fermion has already been removed.

## `Lagrangian/PhysicsNotation.wl`

**Problem solved:** dynamically generated Wolfram/Matchete symbols and indices need deterministic names that remain consistent across construction, export and downstream parsing.

**What it does:** provides common field, coupling and index naming utilities such as `PhysicsField`, `PhysicsIndex`, conjugation/transpose helpers and coupling-base names.

---

# 4. Model-definition files

## `Lagrangian/model/T3ModelCatalog.wl`

**Problem solved:** the allowed T3 representations and hypercharges must be generated systematically rather than hard-coded independently throughout the project.

**What it does:**

- computes T3 hypercharges from $\alpha$;
- tests whether both Yukawa structures are allowed;
- tests whether the scalar-mixing structure is allowed;
- constructs the representation definition for the requested T3 class/model.

## `Lagrangian/model/T3Fields.wl`

**Problem solved:** Matchete needs explicit fields, masses and representations for the selected T3 point.

**What it does:**

- defines the heavy fermion and scalar fields;
- assigns gauge representations;
- constructs the mass specification;
- sets the currently active heavy-field set;
- preserves mass couplings needed across sequential threshold stages.

**Boundary:** field declaration, not interaction selection.

## `Lagrangian/model/SU2Invariants.wl`

**Problem solved:** arbitrary supported SU(2) multiplets need representation-correct invariant tensors.

**What it does:**

- maps dimensions to SU(2) Dynkin labels;
- defines BSM index types;
- constructs the Clebsch-Gordan/invariant tensors needed by the T3 interactions.

## `Lagrangian/model/LagrangianBuilder.wl`

**Problem solved:** several representation-dependent contraction candidates may exist, but the project needs one validated full UV Lagrangian.

**What it does:**

- checks candidate T3 invariant structures;
- selects valid group contractions;
- verifies that required T3 ingredients are present;
- combines fields and interactions into the final Matchete Lagrangian.

---

# 5. Interaction-definition files

## `Lagrangian/interactions/T3Topology.wl`

**Problem solved:** the two Yukawa vertices and scalar/Higgs mixing vertex depend on representation and conjugation choices.

**What it does:**

- constructs candidate $y_1$ contractions;
- constructs candidate $y_2$ contractions;
- constructs candidate $\lambda_{T3}$ mixing contractions;
- packages valid T3 interaction candidates for the builder.

## `Lagrangian/interactions/ScalarPotential.wl`

**Problem solved:** the renormalisable scalar sector must include all representation-required self, portal, adjoint and cross invariants needed for RGE closure.

**What it does:**

- builds scalar norms/self quartics;
- builds Higgs portals;
- builds inter-scalar portals;
- constructs and symmetrises tensor structures;
- finds independent invariant tensors;
- defines dynamic real quartic couplings where the representation requires them.

Missing allowed quartics here would make the later RGE system incomplete.

## `Lagrangian/interactions/ExportUVCGRegistry.wl`

**Problem solved:** downstream Python RGE/group-factor code needs the exact invariant tensors used by the Wolfram UV model.

**What it does:** exports the UV Clebsch-Gordan/invariant registry in a machine-readable form.

This avoids reconstructing the same group contractions independently in Python.

---

# 6. Closely related model-selection files

These are outside `Lagrangian/`, but they are part of the model-building boundary.

## `common/T3Model.py`

**Problem solved:** fast representation/topology validity and neutrality checks are needed before Matchete.

**What it does:** implements topology dimension rules, neutral-component detection, class identification and shared-scalar formal-dimension mapping.

## `common/T3Fields.py`

**Problem solved:** one physical scalar in the shared branch must not accidentally be treated as two independent threshold fields.

**What it does:** maps between physical heavy fields and formal topology roles.

## `model/T3Study.py`

**Problem solved:** study selection and model construction are separate concerns.

**What it does:** decides which T3 points should be built and hands those requests to `Lagrangian/T3ModelMatching.py`.

---

# 7. Matching outputs consumed downstream

The matching subsystem produces the data required by the RGE and numerical layers, including:

- matched EFT expressions;
- extracted one-loop Weinberg information;
- threshold/stage summaries;
- physical and formal field metadata;
- representation/invariant information;
- intermediate Wilson seeds for fermion-first sequential matching;
- paths recorded in `RunRecord`/`EFTStageRecord`.

The final hierarchical physical $C_5$ is completed later by

```text
RGE/matching/FinalWeinbergCoefficient.py
```

after direct intermediate running is known.

---

# 8. Production threshold scope

The model/matching layer can represent threshold plans generically, but current production execution supports:

```text
common threshold
F -> (S1,S2)
F -> S       [shared scalar]
```

Scalar-first paths remain outside production scope because the first scalar threshold can generate a leading dimension-six intermediate operator.

That restriction is a physics truncation decision, not a limitation of the threshold data structures.

---

# 9. Overall Lagrangian subsystem overview

The subsystem can be read as a sequence of questions:

```text
common/T3Model.py
    Is this representation physically/topologically allowed?

model/T3Study.py
    Do we want to calculate this point in the current study?

Lagrangian/ModelValidation.py
    Is the requested run supported and how should it be named?

Lagrangian/model/ + interactions/
    What are the actual fields, invariants and UV interactions?

RunModel.wl / RunMatching.wl / RunThresholdStage.wl
    What EFT is generated when the requested heavy fields are integrated out?

Lagrangian/MatchingResults.py
    How do those external results become stable Python pipeline metadata?
```

The architectural rule is:

> Python owns study selection, physical threshold configuration and process/result boundaries; Wolfram/Matchete owns the explicit UV invariant construction and matching algebra.

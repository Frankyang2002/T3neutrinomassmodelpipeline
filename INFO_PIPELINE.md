# T3 Neutrino-Mass Pipeline

This file is the repository-wide map of the implemented T3 radiative neutrino-mass calculation.  The central executable remains `pipeline.py`; detailed algebra, model construction, running, validation, and report generation are delegated to specialised modules.

## 1. Physics scope and conventions

The code uses the Standard Model hypercharge convention

$$
Q=T_3+Y.
$$

The T3 classification parameter $\alpha$ is related to the physical hypercharges by

$$
Y(S_1)=\frac{\alpha}{2},\qquad
Y(S_2)=\frac{\alpha+2}{2},\qquad
Y(F)=\frac{\alpha+1}{2}.
$$

This is one half of the hypercharge number used in references which write $Q=T_3+Y/2$.  In particular, the Restrepo-Zapata-Yaguna classification uses the doubled convention

$$
Y_{\rm RZY}=2Y.
$$

The historical T3 representation classes are

| class | $d_{S_1}$ | $d_{S_2}$ | $d_F$ |
|---|---:|---:|---:|
| A | 1 | 3 | 2 |
| B | 2 | 2 | 1 |
| C | 2 | 2 | 3 |
| D | 3 | 1 | 2 |
| E | 3 | 3 | 2 |

Normal production mode is restricted to the currently supported singlet/doublet/triplet representations.  `--force` can bypass this production-support restriction only when the representation assignment still satisfies the T3 topology conditions.

The generic T3 interactions are schematically

$$
y_1LFS_1+y_2LFS_2+\lambda_{T3}HHS_1S_2^\dagger+\text{h.c.},
$$

with representation-dependent conjugations and Clebsch-Gordan contractions supplied by the Wolfram/Matchete model builder.

The project uses the Weinberg-to-neutrino-mass convention

$$
m_\nu=-v_{174}^2C_5=-\frac{v_{246}^2}{2}C_5.
$$

## 2. Physical heavy fields

Ordinary T3 models contain the physical heavy fields

```text
F, S1, S2
```

The shared-scalar/scotogenic branch contains

```text
F, S
```

where the one physical scalar `S` fills both formal topology roles `S1` and `S2` up to conjugation.  The authoritative physical/formal translation is centralised in `common/T3Fields.py`.  Threshold planning and RGE dispatch use physical fields; expansion back to formal `S1,S2` roles occurs only at the matching boundary.

For the T3-B scotogenic point with $\alpha=-1$, the formal hypercharges are

$$
Y(S_1)=-\frac12,\qquad Y(S_2)=+\frac12,\qquad Y(F)=0,
$$

and the two scalar roles correspond to one physical doublet and its conjugate.

## 3. Central calculation order

`pipeline.py` is intentionally kept as the authoritative map of the calculation.  Its implemented order is

```text
CLI/configuration
    -> PipelinePlan
    -> model-point selection
    -> UV model construction + threshold matching
    -> attach physical EFT-stage metadata
    -> organise matched C5
    -> UV RGE
    -> each implemented intermediate EFT running interval
    -> independent validation of completed intermediate artifacts
    -> final SMEFT C5 RGE
    -> full-flavor Weinberg RGE
    -> symbolic Majorana mass matrix
    -> optional numerical running + neutrino observables
    -> summaries and reports
```

The detailed implementations live outside `pipeline.py`; the backbone owns their order and failure propagation.

## 4. Threshold plans and EFT content

`common/PipelinePlan.py` is the Python source of truth for the threshold sequence.  It contains an ordered set of `ThresholdStep` objects and derives:

- the active heavy fields before and after each threshold;
- the intermediate `EFTRunningInterval` objects;
- threshold labels and scales;
- stage metadata used by reports;
- the explicit EFT operator-dimension truncation.

`common/EFT.py` contains the generic EFT-content objects.  New calculation code should identify an intermediate theory by its active physical fields, not by ordinal labels such as `EFT1` or `EFT2`.

If no `--threshold` options are supplied, all physical heavy fields are integrated out at one common threshold.  A verified hierarchical fermion-first example is

```powershell
python pipeline.py --dims 2 2 1 --alpha -1 `
    --threshold F `
    --threshold S1 S2
```

For a shared scalar,

```powershell
python pipeline.py --dims 2 1 `
    --threshold F `
    --threshold S
```

### Scalar-first scope

Scalar-first threshold orderings are intentionally outside the production scope used for this project. Integrating a scalar first can generate an intermediate dimension-six operator schematically

$$
\frac{y\lambda_{T3}}{M_S^2}\,L F H H S,
$$

which can be the leading EFT representation of the T3 amplitude at that stage. A $d\le5$ scalar-first calculation would therefore be incomplete rather than merely a smaller correction.

The user-facing pipeline now rejects scalar-first and partially split scalar hierarchies before matching. The supported production choices are:

```text
common threshold:          (F,S1,S2)   or (F,S)
verified hierarchy:        F -> (S1,S2) or F -> S
```

`PipelinePlan` remains field-content based internally, so the architecture does not need to be rewritten if a higher-dimensional scalar-first treatment is added in a future project. No scalar-first $d=6$ production code is retained in the current source tree.

## 5. Model selection and UV matching

`model/T3Study.py` decides which model points are requested.  It owns the smoke/hypercharge/dimension study definitions and the direct `--dims` request object.  It does not construct a Lagrangian.

`Lagrangian/T3ModelMatching.py` is the Python compatibility boundary into `Lagrangian/Runner.py`.  `Runner.py` validates the requested T3 topology, converts the physical threshold plan to the formal Wolfram roles, launches the existing Wolfram/Matchete calculation, and returns `RunRecord` metadata.

The Wolfram side constructs the UV theory and performs matching through files under `Lagrangian/`, including:

- `model/T3ModelCatalog.wl`;
- `model/T3Fields.wl`;
- `model/SU2Invariants.wl`;
- `model/LagrangianBuilder.wl`;
- `interactions/T3Topology.wl`;
- `interactions/ScalarPotential.wl`;
- `RunMatching.wl`;
- `RunThresholdStage.wl`.

The production matching target is fixed at dimension five at one loop.

## 6. UV and intermediate-EFT running

`RGE/running/UVRunning.py` is the high-level UV-RGE entry point.

`RGE/running/IntermediateEFTRunning.py` is a generic field-content dispatcher.  It receives each `EFTRunningInterval` from `PipelinePlan` and selects a production backend from the actual active heavy fields.

The currently verified nontrivial hierarchical backend is

```text
RGE/running/backends/ScalarOnlyAfterFermion.py
```

and applies only when:

1. `F` is integrated out on entry;
2. the complete physical scalar sector remains active;
3. all remaining physical T3 scalars are integrated out together on exit.

The production backend calls descriptive interfaces under `RGE/running/intermediate/`.  The scalar-only implementation now lives directly in physically named modules, including

```text
RGE/running/intermediate/ScalarOnlyTensorAdapters.py
RGE/running/intermediate/ScalarOnlyWilsonTensorRGE.py
RGE/running/intermediate/ScalarOnlyWilsonFlow.py
RGE/running/intermediate/DirectWeinbergRunning.py
RGE/running/intermediate/ScalarThresholdMatching.py
```

The former ordinal `EFT1` implementation modules are no longer production dependencies.  Existing JSON keys, filenames, and a small number of Wolfram symbols containing `EFT1` are retained only as serialization/report compatibility names; they are not the dispatch model used by new production code.

## 7. Validation versus production

Production code constructs the authoritative EFT result.  Independent scientific checks are run afterwards under `validation/`.

For the scalar-only fermion-first interval, the validation backend checks items such as:

- the independent one-generation Wilson transport regression;
- Weinberg-subspace diagnostics emitted by the component RGE;
- one-generation reduction of the full-flavor seed;
- equal-scale vanishing of direct running;
- pole/RGE consistency metadata;
- presence of the final authoritative C5 artifact.

Some consistency conditions remain operational prerequisites inside historical production modules because the present construction of the authoritative C5 depends on them.  They should not be removed merely to make the validation split cosmetically complete.

## 8. Final Weinberg coefficient and neutrino physics

`physics/LowEnergyNeutrino.py` owns the high-level post-matching stages while `pipeline.py` keeps their order explicit.

The sequence is

```text
matched/final C5
    -> one-generation SMEFT Weinberg RGE output
    -> full-flavor Weinberg RGE
    -> symbolic Majorana mass matrix
    -> optional numerical SM + C5 running
    -> neutrino observables
```

The three-generation SMEFT equation used by `RGE/running/weinberg/WeinbergRunning.py` is

$$
16\pi^2\frac{dC_5}{d\ln\mu}
=
(2\lambda_H-3g_2^2+2T)C_5
-\frac32\left[
Y_eY_e^\dagger C_5+C_5(Y_eY_e^\dagger)^T
\right],
$$

with

$$
T=\operatorname{Tr}(Y_eY_e^\dagger+3Y_uY_u^\dagger+3Y_dY_d^\dagger).
$$

The numerical stage evolves the SM parameters and the complex symmetric `C5` matrix and then uses `RGE/phenomenology/NeutrinoObservables.py` for the Takagi factorisation and observable extraction.

## 9. Reports and stable compatibility outputs

`Reports/PipelineReports.py` owns report orchestration.  Detailed LaTeX construction remains in the existing `Reports/` modules.

Per-study raw outputs are written below

```text
output/<study>/
```

and human-readable reports below

```text
Reports/output/<study>/
```

The aggregate comparison file is

```text
output/<study>/t3_model_comparison.json
```

and `--full` writes the master study summary below `output/full/`.

Historical report paths, summary keys, and serialized filenames are intentionally retained where downstream tooling already depends on them.  In particular, `EFT1...` compatibility keys can still appear in JSON even though new orchestration identifies the theory by field content.

## 10. Main repository map

| Area | Main responsibility |
|---|---|
| `pipeline.py` | authoritative calculation order and failure propagation |
| `common/PipelineCLI.py` | CLI configuration |
| `common/PipelinePlan.py` | threshold sequence and EFT intervals |
| `common/T3Fields.py` | physical/formal heavy-field identity |
| `common/RunRecords.py` | run/stage metadata |
| `model/T3Study.py` | requested T3 model points |
| `Lagrangian/T3ModelMatching.py` | Python boundary to UV construction/matching |
| `Lagrangian/` | Wolfram/Matchete UV model and matching |
| `RGE/running/UVRunning.py` | UV RGE orchestration |
| `RGE/running/IntermediateEFTRunning.py` | intermediate field-content dispatch |
| `RGE/running/backends/` | concrete intermediate-EFT production backends |
| `RGE/running/intermediate/` | scalar-only intermediate-EFT tensor, Wilson-flow, and threshold implementations |
| `RGE/general/` | generic scalar, fermion, gauge, anomalous-dimension, and dimension-five Wilson-tensor RGE machinery |
| `RGE/running/weinberg/` | final C5 construction and SMEFT running |
| `Numerical/IntermediateScalarState.py` | numerical state for the post-F scalar-only intermediate EFT |
| `physics/LowEnergyNeutrino.py` | post-matching neutrino calculation orchestration |
| `validation/` | independent scientific/regression checks |
| `Reports/` | report generation only |
| `studies/FullT3Study.py` | `--full` study orchestration |
| `tests/` | architecture and physics regressions |

## 11. Recommended checks after a refactor

Fast Python regression:

```powershell
python -m pytest -q
```

Then run a real Wolfram/Matchete smoke study from the actual project checkout:

```powershell
python pipeline.py --smoke
```

For hierarchical fermion-first matching/running, also run an explicit sequential point such as

```powershell
python pipeline.py --dims 2 2 1 --alpha -1 `
    --threshold F `
    --threshold S1 S2
```

The Python architecture tests cannot replace the real Wolfram/Matchete calculation; both levels are needed before treating a structural refactor as numerically verified.

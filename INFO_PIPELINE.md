# T3 Neutrino-Mass Pipeline

This file is the repository-wide map of the implemented T3 radiative neutrino-mass calculation. `pipeline.py` remains the authoritative execution backbone; detailed model construction, matching, running, numerical integration, neutrino physics, validation, and reporting are delegated to specialised modules.

## 1. Physics scope and conventions

The project uses

$$
Q=T_3+Y.
$$

For the Restrepo-Zapata-Yaguna T3 classification parameter $\alpha$,

$$
Y(S_1)=\frac{\alpha}{2},\qquad
Y(S_2)=\frac{\alpha+2}{2},\qquad
Y(F)=\frac{\alpha+1}{2}.
$$

The classification convention is doubled:

$$
Y_{\rm RZY}=2Y.
$$

The project Weinberg-to-mass convention is

$$
m_\nu=-v_{174}^2C_5=-\frac{v_{246}^2}{2}C_5.
$$

Ordinary T3 models have physical heavy fields

```text
F, S1, S2
```

while the shared-scalar/scotogenic branch has

```text
F, S
```

with one physical scalar `S` filling the formal topology roles `S1` and `S2` up to conjugation. `common/T3Fields.py` is the authoritative physical/formal translation. Expansion back to formal roles occurs only at the matching boundary.

## 2. Central calculation order

`pipeline.py` owns the calculation order and failure propagation:

```text
CLI/configuration
    -> common/PipelinePlan.py
    -> model/T3Study.py
    -> Lagrangian/T3ModelMatching.py
    -> UV construction + threshold matching
    -> UV RGE
    -> implemented intermediate EFT interval(s)
    -> validation/
    -> authoritative final C5
    -> final SM+Weinberg RGE
    -> C5 -> m_nu
    -> neutrino observables
    -> reports
```

Detailed stage implementations should not be reintroduced into `pipeline.py`.

## 3. Threshold plans and production scope

`common/PipelinePlan.py` represents the ordered physical thresholds and derives the active-field `EFTRunningInterval` objects. New code identifies an intermediate theory by physical field content rather than ordinal labels such as `EFT1`.

Supported production choices are

```text
common threshold:          (F,S1,S2) or (F,S)
verified hierarchy:        F -> (S1,S2) or F -> S
```

Scalar-first orderings are outside current production scope. Integrating a scalar first can generate a leading intermediate dimension-six operator schematically

$$
\frac{y\lambda_{T3}}{M_S^2}LFHHS,
$$

so a scalar-first $d\le5$ truncation would omit the leading EFT path. The CLI therefore rejects unsupported scalar-first and partially split hierarchies before Wolfram/Matchete is launched.

## 4. Lagrangian and matching architecture

`model/T3Study.py` selects model requests. `Lagrangian/T3ModelMatching.py` is the compatibility-facing Python boundary into the matching stack.

The Python matching responsibilities are split as follows:

```text
Lagrangian/ModelValidation.py
    T3 request validation, model naming, run specification

Lagrangian/WolframRunner.py
    wolframscript process execution, streamed output, logs

Lagrangian/MatchingResults.py
    summary loading/physicalisation, stage validation, RunRecord construction

Lagrangian/Runner.py
    thin orchestration and historical public entrypoints
```

The Wolfram model and matching definitions remain under `Lagrangian/model/`, `Lagrangian/interactions/`, `RunMatching.wl`, and `RunThresholdStage.wl`.

## 5. UV and intermediate running

`RGE/running/UVRunning.py` orchestrates the UV RGBeta stage.

`RGE/running/IntermediateEFTRunning.py` dispatches each nonzero interval by active field content. The verified fermion-first scalar-only backend is

```text
RGE/running/backends/ScalarOnlyAfterFermion.py
```

and calls physically named interfaces under `RGE/running/intermediate/`, including

```text
ScalarOnlyRenormalisableRunning.py
ScalarOnlyWilsonTensorRGE.py
ScalarOnlyWilsonFlow.py
DirectWeinbergRunning.py
ScalarThresholdMatching.py
```

Historical `EFT1...` filenames and JSON keys remain only where required by compatibility contracts.

## 6. Final C5 ownership

The authoritative hierarchical final coefficient is matching/bookkeeping, not an RGE integrator.

```text
RGE/matching/FinalWeinbergCoefficient.py
    hard threshold + direct running combination
    MSbar/pole bookkeeping
    Majorana symmetrisation
    final_weinberg_coefficient.json

RGE/matching/WeinbergFlavorMatching.py
    one-generation Matchete coefficient -> full flavor

RGE/matching/FinalWeinbergAdapter.py
    hierarchical final-C5 JSON -> symbolic physical Majorana C5
```

Historical paths under `RGE/running/weinberg/` and `RGE/matching/FlavorC5Matching.py` are compatibility boundaries only.

## 7. Final SM+Weinberg RGE and neutrino physics

The analytic final-EFT beta model is

```text
RGE/running/weinberg/WeinbergRGE.py
```

and implements

$$
16\pi^2\frac{dC_5}{d\ln\mu}
=
(2\lambda_H-3g_2^2+2T)C_5
-\frac32\left[
Y_eY_e^\dagger C_5+C_5(Y_eY_e^\dagger)^T
\right].
$$

The symbolic full-flavor output stage is

```text
RGE/running/weinberg/FullFlavorWeinbergStage.py
```

and the one-generation tensor cross-check is explicitly named

```text
RGE/running/weinberg/OneGenerationWeinbergBenchmark.py
RGE/running/weinberg/OneGenerationWeinbergBenchmarkStage.py
```

Physical interpretation lives under `physics/`:

```text
physics/NeutrinoMass.py
physics/NeutrinoObservables.py
physics/NeutrinoTrajectory.py
physics/LowEnergyNeutrino.py
```

## 8. Numerical ownership

Numerical ODE/state/trajectory work stays under `Numerical/`.

```text
Numerical/SMWeinbergEvolution.py
    final SM+C5 state packing and solve_ivp evolution

Numerical/SMWeinbergStage.py
    numerical C5 evaluation, integration-stage serialization

Numerical/FermionThresholdBoundary.py
    UV -> scalar-only renormalisable boundary

Numerical/ScalarThresholdBoundary.py
    scalar-only -> final-SM boundary

Numerical/T3Trajectory.py
    multi-segment UV/intermediate/final trajectory orchestration

Numerical/FinalC5TrajectoryAdapter.py
    trajectory -> authoritative final-C5 numerical evaluator inputs
```

Historical modules such as `Numerical/WeinbergRunning.py`, `Numerical/WeinbergStage.py`, `Numerical/ThresholdMatching.py`, `Numerical/FinalSMBoundary.py`, and `Numerical/FinalC5Bridge.py` are retained as compatibility surfaces.

## 9. Validation and reports

`validation/` contains independent scientific/regression checks only. `Reports/PipelineReports.py` owns report orchestration, while detailed report construction stays under `Reports/`.

Stable serialized names and historical summary keys are retained where downstream tooling already depends on them.

## 10. Main repository map

| Area | Main responsibility |
|---|---|
| `pipeline.py` | authoritative calculation order and failure propagation |
| `common/PipelinePlan.py` | thresholds and EFT intervals |
| `common/T3Fields.py` | physical/formal heavy-field identity |
| `model/T3Study.py` | requested model points |
| `Lagrangian/T3ModelMatching.py` | compatibility boundary into matching |
| `Lagrangian/ModelValidation.py` | model validation/specification |
| `Lagrangian/WolframRunner.py` | external Wolfram process execution |
| `Lagrangian/MatchingResults.py` | matching outputs and RunRecord construction |
| `RGE/running/IntermediateEFTRunning.py` | intermediate-EFT dispatch |
| `RGE/running/backends/ScalarOnlyAfterFermion.py` | verified hierarchical backend |
| `RGE/running/intermediate/ScalarOnlyWilsonTensorRGE.py` | scalar-only dimension-five tensor RGE |
| `RGE/matching/FinalWeinbergCoefficient.py` | authoritative final-C5 construction |
| `RGE/running/weinberg/WeinbergRGE.py` | final SM+Weinberg beta model |
| `Numerical/IntermediateScalarState.py` | scalar-only numerical state |
| `Numerical/SMWeinbergEvolution.py` | final SM+Weinberg numerical integration |
| `physics/LowEnergyNeutrino.py` | post-matching neutrino orchestration |
| `validation/` | independent checks only |
| `Reports/PipelineReports.py` | report orchestration |

## 11. Regression procedure

Run the Python suite first:

```powershell
python -m pytest -q
```

Then exercise the actual external calculation:

```powershell
python pipeline.py --smoke
```

For the verified hierarchical benchmark also run

```powershell
python pipeline.py --dims 2 2 1 --alpha -1 `
    --threshold F `
    --threshold S1 S2
```

Python tests do not replace a real Wolfram/Matchete/RGBeta calculation.

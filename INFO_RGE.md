# T3 RGE and Neutrino-Mass Pipeline

This document describes the implemented running, threshold transport, final Weinberg coefficient, and low-energy neutrino stages.

## 1. Conventions

The project uses

$$
Q=T_3+Y,\qquad t=\ln\mu,
$$

and

$$
m_\nu=-\frac{v^2}{2}C_5=-v_{174}^2C_5.
$$

The final SM+Weinberg beta model is implemented in

```text
RGE/running/weinberg/WeinbergRGE.py
```

with

$$
16\pi^2\frac{dC_5}{d\ln\mu}
=
(2\lambda_H-3g_2^2+2T)C_5
-\frac32\left[
Y_eY_e^\dagger C_5+C_5(Y_eY_e^\dagger)^T
\right],
$$

where

$$
T=\operatorname{Tr}
\left(
Y_eY_e^\dagger+3Y_uY_u^\dagger+3Y_dY_d^\dagger
\right).
$$

## 2. High-level RGE flow

```text
UV T3
    -> UV renormalisable running
    -> fermion threshold
    -> scalar-only intermediate EFT
    -> direct LLSS -> Weinberg running
    -> scalar threshold
    -> authoritative final C5
    -> final SM+Weinberg running
    -> m_nu
    -> observables
```

`pipeline.py` owns the global sequence. RGE modules own equations and stage-specific calculations.

## 3. UV running

`RGE/running/UVRunning.py` is the production UV entry point and uses the RGBeta stack under `RGE/running/rgbeta/`.

Numerically, the UV state and `solve_ivp` integration live under `Numerical/State.py`, `Numerical/StateVector.py`, and `Numerical/UVRunner.py`.

## 4. Intermediate scalar-only EFT

The production dispatcher is

```text
RGE/running/IntermediateEFTRunning.py
```

and the verified hierarchical backend is

```text
RGE/running/backends/ScalarOnlyAfterFermion.py
```

for

$$
F\rightarrow(S_1,S_2)
$$

or $F\rightarrow S$ in the shared-scalar branch.

The production sequence is:

1. scalar-only renormalisable running;
2. dimension-five LLSS Wilson RGE;
3. full-flavor tree-level Wilson boundary;
4. direct one-loop LLSS $\to$ Weinberg mixing;
5. resumed scalar-threshold matching;
6. authoritative final-C5 construction.

The fixed-order reason for excluding LLSS self-running from the authoritative one-loop C5 is

$$
C_{\rm LLSS}^{(0)}
\xrightarrow{\text{1-loop self-running}}
C_{\rm LLSS}^{(1)}
\xrightarrow{\text{scalar loop}}
C_5^{(2)},
$$

which first contributes at two-loop order.

## 5. Final C5 is matching/bookkeeping

The canonical implementation is

```text
RGE/matching/FinalWeinbergCoefficient.py
```

not a running module. It combines

$$
C_5^{pq}(M_S)
=
C_{5,\mathrm{hard}}^{pq}(M_F)
+
\hbar\,\Delta C_{5,\mathrm{direct}}^{(1),pq}(M_F\to M_S),
$$

performs the existing MSbar/pole-consistency bookkeeping, and constructs the physical symmetric Majorana coefficient.

Flavor lifting is split into

```text
RGE/matching/WeinbergFlavorMatching.py
RGE/matching/FinalWeinbergAdapter.py
```

Historical files such as `RGE/running/weinberg/FinalWeinbergCoefficient.py` and `RGE/matching/FlavorC5Matching.py` are compatibility surfaces.

## 6. Symbolic final-EFT stages

The symbolic three-generation stage is

```text
RGE/running/weinberg/FullFlavorWeinbergStage.py
```

and calls the analytic model in `WeinbergRGE.py`.

The former ambiguous "matched Weinberg RGE" name has been replaced by an explicit one-generation benchmark:

```text
RGE/running/weinberg/OneGenerationWeinbergBenchmark.py
RGE/running/weinberg/OneGenerationWeinbergBenchmarkStage.py
```

Historical `MatchedWeinbergRGE.py`, `MatchedWeinbergStage.py`, and `FlavorMatchedWeinbergStage.py` files remain compatibility-only.

## 7. Numerical final SM+Weinberg evolution

The canonical numerical modules are

```text
Numerical/SMWeinbergEvolution.py
Numerical/SMWeinbergStage.py
Numerical/WeinbergTrajectory.py
```

`SMWeinbergEvolution.py` owns state packing/unpacking and `solve_ivp`. It does not own the physical $C_5\to m_\nu$ interpretation.

Threshold-state construction is split into

```text
Numerical/FermionThresholdBoundary.py
Numerical/ScalarThresholdBoundary.py
```

and `Numerical/T3Trajectory.py` connects UV, intermediate, and final segments.

The final-C5 numerical adapter is

```text
Numerical/FinalC5TrajectoryAdapter.py
```

which extracts the values at the correct threshold scales and calls the established final-C5 evaluator. `Numerical/FinalC5Bridge.py` remains a compatibility wrapper for historical callers and mocks.

## 8. Neutrino physics ownership

Physical interpretation is under `physics/`:

```text
physics/NeutrinoMass.py
    symbolic and numerical C5 -> m_nu

physics/NeutrinoObservables.py
    Takagi factorisation, masses, mass splittings, |U_PMNS|

physics/NeutrinoTrajectory.py
    charged-lepton basis rotation and scale-dependent observables

physics/LowEnergyNeutrino.py
    orchestration only
```

The oscillation-fit numerical layer imports the physical observables from `physics.NeutrinoObservables` directly.

## 9. Production versus validation

Independent scientific checks live under `validation/`. The production backend may still consume a validated pole/RGE relation where the established final-C5 construction requires it; that is an operational prerequisite rather than a reason to move independent validation logic back into production.

## 10. Scalar-first scope

Scalar-first matching is outside production scope because a first scalar threshold can generate a leading dimension-six operator schematically

$$
\frac{y\lambda_{T3}}{M_S^2}LFHHS.
$$

Dropping it would remove the leading scalar-first EFT path. Current production therefore supports only the common threshold and the verified fermion-first hierarchy.

## 11. Regression procedure

Run:

```powershell
python -m pytest -q
python pipeline.py --smoke
```

and for the hierarchical benchmark:

```powershell
python pipeline.py --dims 2 2 1 --alpha -1 `
    --threshold F `
    --threshold S1 S2
```

The real Wolfram/Matchete/RGBeta run is required in addition to Python architecture tests.

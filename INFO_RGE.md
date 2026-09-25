# T3 RGE and Neutrino-Mass Pipeline

This document describes the implemented RGE path after the T3 UV model has been built and matched.  `INFO_PIPELINE.md` is the repository-wide map; this file focuses on running, threshold transport, the final Weinberg coefficient, and neutrino observables.

## 1. Conventions

The code uses

$$
Q=T_3+Y,
\qquad t=\ln\mu.
$$

The matched Weinberg coefficient is related to the Majorana neutrino mass by

$$
m_\nu=-\frac{v^2}{2}C_5
$$

when $v\simeq246\,\mathrm{GeV}$, equivalently $m_\nu=-v_{174}^2C_5$.

The final SMEFT Weinberg equation implemented in `RGE/running/weinberg/WeinbergRunning.py` is

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
T=\operatorname{Tr}\left(
Y_eY_e^\dagger+3Y_uY_u^\dagger+3Y_dY_d^\dagger
\right).
$$

## 2. RGE architecture

The implemented high-level flow is

```text
UV T3 model
    -> UV renormalisable RGEs
    -> threshold matching
    -> optional intermediate heavy-field EFT
    -> final hard + running C5
    -> SMEFT Weinberg RGE
    -> full-flavor C5
    -> neutrino mass matrix
    -> optional numerical evolution
    -> observables
```

`pipeline.py` owns this order.  RGE modules calculate individual stages but do not decide the global execution sequence.

## 3. UV running

`RGE/running/UVRunning.py` is the high-level Python entry point for the UV theory.  It uses the existing RGBeta machinery under `RGE/running/rgbeta/`.

The UV stage is always associated with the complete physical T3 heavy-field content from `RunRecord.physical_heavy_fields`.

Normal production support remains tied to the representation range supported by the present T3/RGBeta model setup.  Larger topology-compatible representations can be attempted with `--force`, but that does not imply the downstream RGE implementation supports them.

## 4. Generic intermediate-EFT dispatch

The threshold sequence is represented by `PipelinePlan`.  For every nonzero interval between two thresholds it produces an `EFTRunningInterval` containing:

- the physical heavy fields active in that interval;
- the threshold transition which entered the interval;
- the threshold transition which exits it;
- the high and low matching scales.

`RGE/running/IntermediateEFTRunning.py` dispatches from that physical content.  Ordinal names such as `EFT1` are not used to decide which calculation applies.

If no production backend exists for an interval, the outcome is explicitly

```text
NotImplementedForFieldContent
```

and the main pipeline marks the route incomplete.  The low-energy neutrino stages are then blocked for that run.

## 5. Verified fermion-first scalar-only backend

The implemented nontrivial hierarchical backend is

```text
RGE/running/backends/ScalarOnlyAfterFermion.py
```

with backend identifier

```text
scalar_only_after_fermion_d5
```

It applies when

$$
\text{UV}\xrightarrow{\;F\;}
\text{SM + complete T3 scalar sector}
\xrightarrow{\;\text{all remaining scalars}\;}
\text{SMEFT}.
$$

For ordinary T3 this means the intermediate active heavy fields are $\{S_1,S_2\}$.  For the shared-scalar branch the active set is $\{S\}$.

The production backend performs, in order:

1. the renormalisable scalar-EFT RGE;
2. the dimension-five $\psi^2\phi^2$ Wilson RGE;
3. construction of the full-flavor tree-level Wilson boundary;
4. direct one-loop mixing into the Weinberg operator;
5. resumed scalar-threshold matching;
6. construction of the final authoritative Weinberg coefficient.

The descriptive interfaces are under `RGE/running/intermediate/`.

## 6. Fixed-order bookkeeping

After integrating out the heavy fermion, a tree-level LLSS-type dimension-five Wilson coefficient is present.  Its one-loop direct mixing into the Weinberg operator contributes at the one-loop order retained by the project.

The backend deliberately does not feed one-loop self-running of the heavy LLSS coefficient back through a scalar loop into the authoritative one-loop C5, because that combination would be of order $\hbar^2$.

This fixed-order separation is why the code distinguishes:

- hard threshold matching;
- direct LLSS $\to$ Weinberg running;
- diagnostic heavy-operator transport.

## 7. Scalar-only intermediate-EFT implementation

The original implementation called the scalar-only fermion-first theory `EFT1`.  The refactored production API no longer depends on that stage number, and the implementation now lives directly under `RGE/running/intermediate/`.

The main descriptive interfaces are:

```text
run_scalar_only_renormalisable_rge
run_scalar_only_dimension_five_wilson_rge
export_full_flavor_wilson_boundary
build_direct_weinberg_running
resume_scalar_threshold_with_running
```

The component Wilson calculation itself is implemented in `ScalarOnlyWilsonTensorRGE.py`, with tensor construction in `ScalarOnlyTensorAdapters.py` and transport logic in `ScalarOnlyWilsonFlow.py`.

Existing `EFT1...` JSON keys, filenames, and Wolfram insertion symbols are preserved for report and regression compatibility.  They are serialized legacy names, not the conceptual EFT identifier used by production dispatch.

## 8. General tensor RGE machinery

The intermediate Wilson calculation ultimately uses the generic tensor infrastructure under `RGE/general/`:

- `ScalarBasis.py` defines real scalar blocks;
- `GaugeGenerators.py` builds gauge generators;
- `FermionBasis.py` builds the Weyl-fermion basis and gauge sectors;
- `WilsonTensorRGE.py` evaluates the one-loop $\psi^2\phi^2$ tensor terms;
- `Psi2Phi3RGE.py` evaluates the one-loop $\psi^2\phi^3\to\psi^2\phi^3$ self-running terms needed by the scalar-first $d=6$ path;
- `AnomalousDimensions.py` supplies the anomalous-dimension contributions.

The historical scalar-only Wilson module constructs its context from matched Wilson/quartic seeds and RGBeta metadata, projects onto the required coefficient symmetry, evaluates the nonzero tensor beta functions, and records Weinberg-subspace diagnostics.

## 9. Direct Weinberg transport

The current full-flavor direct-running path constructs a symmetric flavor boundary coefficient and derives the direct Weinberg beta in the form

$$
16\pi^2\,\beta_{\kappa,pq}
= r\,C_{12,pq},
$$

where $r$ is the representation/coupling prefactor extracted from the component calculation and $C_{12,pq}$ carries the full flavor structure.

Between two threshold scales $\mu_h$ and $\mu_l$, the fixed-order leading-log correction is proportional to

$$
\ln\!\frac{\mu_l}{\mu_h}.
$$

The equal-scale limit therefore vanishes.  The corresponding equal-scale and one-generation reductions are scientific diagnostics recorded by the validation layer.

## 10. Final C5 construction

`RGE/running/weinberg/FinalWeinbergCoefficient.py` assembles the final Weinberg coefficient after the last heavy threshold.  It combines the hard threshold contribution with the separately carried direct-running contribution and preserves the project's current MSbar/pole-consistency handling.

The pole/RGE relation is currently not merely a cosmetic regression: parts of the historical final-C5 construction use the validated relation to authorize the existing subtraction/normalization step.  That condition therefore remains a production prerequisite until the construction API is redesigned.

## 11. Production versus validation

Independent post-production checks for the verified scalar-only interval are under

```text
validation/backends/ScalarOnlyAfterFermionValidation.py
```

They include the one-generation Wilson-transport regression and interpretation of diagnostic metadata already emitted by production artifacts.

A validation failure still produces a nonzero pipeline status, but the production status is stored separately so it is possible to distinguish:

```text
calculation failed
```

from

```text
calculation completed but an independent cross-check failed
```

## 12. Scalar-first threshold ordering

A scalar-first threshold can generate an intermediate operator of dimension six, schematically

$$
\frac{y\lambda_{T3}}{M_S^2}LFHHS.
$$

After later integration of $F$ and the remaining scalar, this operator can contribute to the leading Weinberg coefficient.  Consequently, a dimension-five-only scalar-first intermediate theory is a truncation of the leading path, not a complete alternative ordering.

The implemented $d=6$ foundation now contains three concrete pieces:

- `--eft-max-dimension 6` propagates the requested operator dimension to the Wolfram/Matchete matching call;
- `RGE/running/intermediate/ScalarFirstDimensionSixSeed.py` extracts exact candidate $LFHHS$ terms from the Matchete tree EFT and preserves them verbatim as a boundary JSON artifact;
- `RGE/general/Psi2Phi3RGE.py` implements the explicit one-loop $\psi^2\phi^3\leftarrow\psi^2\phi^3$ kernel of Eq. (4.25) with the same gauge, Yukawa, quartic, and collinear-anomalous-dimension conventions as the existing general tensor machinery.

The general-EFT source states that this class must additionally be supplemented by the EOM/redundant-$D^2\phi^4$ contribution of Eq. (2.36); that piece is not yet implemented.  This is therefore not yet an authoritative scalar-first production backend.  The exact Matchete boundary still has to be mapped into the symmetric tensor $C_{ijabc}=C_{(ij)(abc)}$, and any additional dimension-six operator mixing/EOM contributions required by that matched boundary must be verified before the later thresholds are allowed to construct the physical $C_5$.  Until then `IntermediateEFTRunning.py` deliberately leaves the interval as `NotImplementedForFieldContent`, so downstream neutrino physics remains blocked.

The old $d\le5$ scalar-first diagnostic is still available only with `--allow-truncated-scalar-first`.

## 13. Low-energy full-flavor and numerical running

After an authoritative final C5 is available, `physics/LowEnergyNeutrino.py` calls the final EFT stages in this order:

1. `RGE/matching/MatchedWeinbergRGE.py` for the existing one-generation SMEFT RGE output/check;
2. `RGE/running/weinberg/FlavorMatchedWeinbergStage.py` for the symbolic three-generation coefficient and beta matrix;
3. the symbolic Majorana mass construction;
4. optionally `NumericalWeinbergStage.py` and `WeinbergRunning.py` for numerical evolution;
5. `RGE/phenomenology/NeutrinoObservables.py` for masses, splittings, and mixing quantities.

The flavor coefficient is complex symmetric, and the neutrino mass matrix is Takagi-factorised rather than diagonalised as a generic Hermitian matrix.

## 14. Regression procedure

After changing RGE architecture, run the full Python suite in the real checkout:

```powershell
python -m pytest -q
```

Then exercise the actual external calculation:

```powershell
python pipeline.py --smoke
```

and one explicit hierarchical fermion-first point:

```powershell
python pipeline.py --dims 2 2 1 --alpha -1 `
    --threshold F `
    --threshold S1 S2
```

The external Wolfram/Matchete/RGBeta smoke run is required to verify the generated physics artifacts; static architecture tests alone cannot establish numerical equivalence.

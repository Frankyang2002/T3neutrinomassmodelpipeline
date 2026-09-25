# T3 UV Lagrangian and Matching

This file documents the UV-model and EFT-matching side of the T3 neutrino-mass pipeline.

## 1. T3 representations and hypercharge

The code uses

$$
Q=T_3+Y.
$$

For the Restrepo-Zapata-Yaguna parameter $\alpha$,

$$
Y(S_1)=\frac{\alpha}{2},\qquad
Y(S_2)=\frac{\alpha+2}{2},\qquad
Y(F)=\frac{\alpha+1}{2}.
$$

The classification convention is

$$
Y_{\rm RZY}=2Y.
$$

The T3 topology requires $d_{S_i}=d_F\pm1$ together with the required scalar/Higgs triplet contraction.

## 2. Physical versus formal scalar roles

Ordinary models use formal and physical fields `F, S1, S2`. The shared-scalar/scotogenic branch contains one physical scalar `S`, which fills the formal topology roles `S1` and `S2` up to conjugation.

`common/T3Fields.py` is the Python source of truth. Threshold planning uses physical fields; expansion to formal `S1,S2` roles occurs only at the matching boundary.

For the T3-B point $\alpha=-1$,

$$
(d_{S_1},d_{S_2},d_F)=(2,2,1),
$$

with

$$
Y(S_1)=-\frac12,\qquad
Y(S_2)=+\frac12,\qquad
Y(F)=0.
$$

## 3. UV interaction content

Schematically,

$$
\mathcal L_{\rm T3}
\supset
y_1LFS_1+y_2LFS_2+\lambda_{T3}HHS_1S_2^\dagger+\mathrm{h.c.}
$$

with representation-dependent conjugations and Clebsch-Gordan contractions supplied by the Wolfram model builder.

## 4. Python matching architecture

`model/T3Study.py` selects requested model points. `Lagrangian/T3ModelMatching.py` preserves the historical public boundary into the matching stack.

Responsibilities are now split explicitly:

```text
Lagrangian/ModelValidation.py
    validates T3 topology/support
    derives names and ModelRunSpecification

Lagrangian/WolframRunner.py
    builds/runs wolframscript commands
    streams process output
    handles logs and interruption

Lagrangian/MatchingResults.py
    loads and physicalises summaries
    validates threshold-stage counts
    constructs RunRecord objects

Lagrangian/Runner.py
    thin orchestration and compatibility entrypoints
```

`Runner.py` no longer owns all validation, subprocess execution, and result construction itself.

## 5. Wolfram model construction

The principal model files remain:

```text
Lagrangian/model/T3ModelCatalog.wl
Lagrangian/model/T3Fields.wl
Lagrangian/model/SU2Invariants.wl
Lagrangian/model/LagrangianBuilder.wl
Lagrangian/interactions/T3Topology.wl
Lagrangian/interactions/ScalarPotential.wl
```

The production target remains

$$
d_{\rm EFT}=5,\qquad L=1.
$$

The Matchete flow is conceptually

```text
UV Lagrangian
    -> Match
    -> GreensSimplify
    -> EOMSimplify
    -> EvaluateLoopFunctions
    -> ReplaceEffectiveCouplings
    -> matched EFT
```

with the Standard-Model-only contribution subtracted from the full result.

## 6. Weinberg extraction and convention

The matching stage extracts one holomorphic Weinberg orientation without double-counting its Hermitian conjugate.

Downstream the project uses

$$
m_\nu=-\frac{v_{246}^2}{2}C_5.
$$

The scotogenic normalization is protected by regression tests; structural cleanup does not modify the one-loop normalization.

## 7. Sequential thresholds

`common/PipelinePlan.py` represents physical threshold orderings and `Lagrangian/RunThresholdStage.wl` performs the requested sequential matching.

Production supports:

```text
common threshold
ordinary hierarchy: F -> (S1,S2)
shared hierarchy:   F -> S
```

Scalar-first orderings are intentionally outside production scope because the first scalar threshold can generate a leading dimension-six operator

$$
\frac{y\lambda_{T3}}{M_S^2}LFHHS.
$$

A $d\le5$ scalar-first calculation would therefore be incomplete.

## 8. Matching data passed downstream

The matching stage exports:

- matched EFT expressions;
- extracted Weinberg coefficient;
- threshold summaries;
- field/representation metadata;
- invariant-tensor information used by RGE adapters;
- tree-level intermediate Wilson seeds for fermion-first sequential matching.

`RunRecord` and `EFTStageRecord` contain metadata and paths; they do not perform the algebra themselves.

## 9. Practical execution

Ordinary benchmark:

```powershell
python pipeline.py --dims 2 2 1 --alpha -1
```

Shared-scalar benchmark:

```powershell
python pipeline.py --dims 2 1
```

Verified hierarchy:

```powershell
python pipeline.py --dims 2 2 1 --alpha -1 `
    --threshold F `
    --threshold S1 S2
```

After a matching-layer change:

```powershell
python -m pytest -q
python pipeline.py --smoke
```

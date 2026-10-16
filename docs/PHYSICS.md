# T3 model physics and symbolic implementation

This document consolidates the technical overview from `INFO_PIPELINE.md`, `INFO_LAGRANGIAN.md`, `INFO_RGE.md`, `INFO_T3_MULTIPLICITY.md` and relevant sections of `INFO_IMPLEMENTATION.md`. For an unabridged file-by-file specification and specialist caveats, consult those original notes until the consolidation has been reviewed.

## 1. Fields, charges and representation scope

The project uses the electric-charge convention

\[
Q=T_3+Y,\qquad
Y(S_1)=\frac{\alpha}{2},\quad
Y(S_2)=\frac{\alpha+2}{2},\quad
Y(F)=\frac{\alpha+1}{2}.
\]

The Restrepo–Zapata–Yaguna convention in some literature is instead \(Y_{\rm RZY}=2Y\). The production model catalogue covers the ordinary T3 classes:

| Class | \((d_{S_1},d_{S_2},d_F)\) |
|---|---|
| A | (1,3,2) |
| B | (2,2,1) |
| C | (2,2,3) |
| D | (3,1,2) |
| E | (3,3,2) |

Production code is restricted to the supported singlet/doublet/triplet cases; topology and neutral-component selection constrain which representation/hypercharge combinations are admitted. The shared-scalar/scotogenic branch and Majorana versus vectorlike fermion structures require their respective conventions and should not be conflated with generic ordinary-T3 output.

Schematically the relevant interactions are

\[
\mathcal L_{\rm T3}\supset y_1 LFS_1+y_2LFS_2
+\lambda_{T3}HH S_1S_2^\dagger+\mathrm{h.c.},
\]

with the conjugations and representation-dependent invariant tensors supplied by the Wolfram model builder. This schematic line is not a substitute for the generated invariant Lagrangian.

## 2. Source ownership

| Subsystem | Responsibility |
|---|---|
| `common/T3Model.py`, `common/T3Fields.py`, `model/T3Study.py` | Physical model identity and study selection |
| `Lagrangian/T3ModelMatching.py`, `Lagrangian/Runner.py` | Orchestration and matching boundaries |
| `Lagrangian/model/`, `Lagrangian/interactions/` | Explicit fields, SU(2) invariants and interaction vertices |
| `Lagrangian/RunModel.wl`, `RunMatching.wl`, `RunThresholdStage.wl` | External Wolfram/Matchete stages |
| `RGE/group_factors/` | Representation/group-factor algebra, recoupling, validation |
| `RGE/running/rgbeta/` | RGBeta model, extraction and beta-function exports |
| `RGE/running/intermediate/` | Scalar-only EFT evolution, LLSS and direct Weinberg mixing |
| `RGE/matching/` | Authoritative matched Weinberg coefficient and flavour/tensor adapters |
| `RGE/running/weinberg/` | Final SM + Weinberg running model |
| `physics/` | Low-energy physical masses, mixing and observables |

Python chooses physical inputs and coordinates processes; Wolfram/Matchete provides generated invariant matching expressions; RGBeta exports the appropriate renormalisable beta functions. **Numerical integrators should evaluate these expressions, not introduce competing matching or RGE formulae.**

## 3. Fermion-first threshold flow and EFTs

The production sequence is:

```text
UV T3 model and RGEs
  -> integrate out F at mu_F
  -> intermediate scalar-only EFT (EFT1)
  -> evolve renormalisable couplings and applicable Wilson structures
  -> scalar threshold at mu_S
  -> form final C5 = hard + direct contribution
  -> run the SM + Weinberg EFT to mu_low
  -> neutrino masses and flavour observables
```

The ordinary documented benchmark hierarchy is \(\mu_{\rm UV}=10^7\,\mathrm{GeV}>\mu_F=10^5\,\mathrm{GeV}>\mu_S=10^3\,\mathrm{GeV}>\mu_{\rm low}=100\,\mathrm{GeV}\). Scalar-first threshold orders and arbitrary split scalar thresholds are **not production-supported** according to the implementation notes. Threshold matching is retained even when a study freezes beta evolution in one stage.

For one-loop running the convention is

\[
t=\ln\mu,\qquad 16\pi^2\,\frac{dX}{dt}=\beta_X^{(1)}.
\]

The final mass convention is

\[
m_\nu=-\frac{v^2}{2}C_5,
\]

with the appropriately symmetrised, flavour-resolved Wilson coefficient and charged-lepton basis prescription.

## 4. LLSS and one-loop power counting

The intermediate LLSS Wilson structure can be shown diagnostically using fixed-order leading-log transport:

\[
C_{\rm LLSS}(\mu)=C_{\rm LLSS}^{(0)}+
\frac{\ln(\mu/M_F)}{16\pi^2}\,\beta^{(1)}[C_{\rm LLSS}^{(0)}]
+\mathcal O(\hbar^2).
\]

This display is **not** the final Weinberg coefficient. Feeding the one-loop self-running correction into an additional scalar loop would be a higher-order contribution under the project's fixed one-loop truncation. The documented authoritative expression is

\[
C_5(M_S)=C_{5,\rm hard}+\Delta C_{5,\rm direct}.
\]

The *direct* \(LLSS\to O_5\) mixing contribution is distinct from treating the fixed-order LLSS self-running diagnostic as a physical Weinberg coefficient. This distinction affects labels on reports and interpretation of ablation studies.

## 5. Charge-channel multiplicity and Class E

The recorded Matchete CG-tensor audit verified 16 models/48 topology tensors under its covariance checks. The compatible internal-charge count in the unbroken-SU(2), multiplet-mass benchmark is one channel for A–D and two equal-weight channels for E.

| Class | Compatible channels | Documented contracted relative channel weights |
|---|---:|---|
| A | 1 | \(-2\sqrt3/3\) |
| B | 1 | \(-2\sqrt3/3\) |
| C | 1 | \(2/3\) |
| D | 1 | \(-2\sqrt3/3\) |
| E | 2 | \(\sqrt6/3\) per channel |

The *benchmark-only* normalisation in `studies/FullT3Study.py` is

```python
LAMBDA_T3_NORMALISATION = {"A": 1.0, "B": 1.0, "C": 1.0, "D": 1.0, "E": 0.5}
```

That halves Class-E **initial** \(\lambda_{T3}\) relative to the reference to compensate the two equal channel contributions in the simplified mass model. It is **not** an extra multiplicity multiplier to apply to the physical neutrino-mass amplitude, nor a replacement for physical Clebsch–Gordan tensors. With charge-dependent mass splittings a channel sum \(\sum_q w_qI_q\) would be needed rather than one common loop function.

## 6. Physics and implementation limits

- The analytic implementation uses a specified one-loop truncation, threshold hierarchy and matching convention, not every possible T3 EFT ordering.
- Underlying heavy flavour matrices and threshold scales are benchmark-dependent. SU(2) multiplet masses in the noted charge audit are not fully electroweak-broken component mass matrices.
- Numerical comparisons should preserve the hypercharge, Yukawa orientation, scalar quartic, Majorana symmetry and mass-basis conventions in `INFO_IMPLEMENTATION.md` until they are fully reproduced in documentation/tests.
- Group-factor and tensor validations constrain algebra, but do not alone establish an experimentally viable neutrino spectrum.

## Source notes

`INFO_PIPELINE.md`, `INFO_LAGRANGIAN.md`, `INFO_RGE.md`, `INFO_T3_MULTIPLICITY.md`, `INFO_IMPLEMENTATION.md`, `docs/archive/INFO_LLSS_DIAGNOSTIC.md`, `docs/archive/INFO_T3_BENCHMARK_VALIDATION.md`.

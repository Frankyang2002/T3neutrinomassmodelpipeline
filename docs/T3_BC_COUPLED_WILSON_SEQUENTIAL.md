# B(-1)/C(-1) coupled EFT1 Wilson coefficients — conditional numerical test

## Scope

This stage replaces the **frozen tree-level LLSS coefficient** in the previous
`RunT3BCExtraWilsonSequential` selected-source approximation with a time-evolved
Wilson tensor. It is a **separate candidate study**, not a replacement for any
published or previous output.

The numerical tensor carries all independent `L1,L2,eC` pairs and all real
scalar pairs, giving **6 × 78 = 468 independent gauge-index components** per
physical lepton-flavour pair. Six physical flavour pairs are evolved. The
numerical RHS implements all **seven terms** of
`RGE/general/WilsonTensorRGE.py`, including terms proportional to `C*`.

Background: previous continuous full-general-scalar EFT1 running, reintegrated
from the saved F-threshold state with the same scalar/RGBeta beta functions.
No UV running is repeated.

F matching: existing one-generation Matchete PL LLSS Wilson tree seed lifted
using positive heavy mass eigenvalues and UV-evolved mass-basis Yukawa tensors;
plus the saved MSbar `C5` hard matching and required grouped-threshold log
relocation. The group-index tensor is normalized using the existing Weinberg
operator template (30 independent real-field components).

Scalar quartics: Matchete scalar quartic seed for the ordinary scalar sector;
**Stage17-certified** B/C two-Higgs extra quartic potential for
`lambdaHHdagS2S2`, `lambdaHHdagS1barS1bar`;
other five additional contractions constructed from Stage22 Matchete tensors
and conditional Stage29 coupling blocks. Stage29 is not a physically certified
Matchete ↔ RGBeta map. It is held conditional in the output.

## Mandatory consistency tests

Before numerical integration, the script:

1. Verifies the saved 16-model scalar threshold report and EFT1 endpoint.
2. Requires the original B/C Stage17/Stage19 tensor and normalization
   certificates, and validates saved direct and selected-extra results.
3. Evolved *frozen-seed* direct Weinberg source, projected from the numerical
   seven-term master, must reproduce the independently saved hierarchical
   direct-only `C5` source to relative tolerance `3e-4`.
4. The projected beta from the two active additional quartics must reproduce
   the independently saved Stage17/19 extra source to relative tolerance
   `3e-4`.
5. The scalar-threshold `LLHH` tensor must remain on the Weinberg subspace,
   with projection residual < `2e-5`.

If one of these checks fails, the script **does not produce neutrino masses**.
Do not relax the checks to make a run pass; review the tensor conventions.

## Install

Extract the ZIP in the local repository root
`C:\thesistemp\T3neutrinomassmodelpipeline`, preserving paths.

- `Numerical/experimental/T3BCCoupledWilsonTensor.py` — numerical master,
  real-basis quartic tensor compilation and ODE
- `scripts/RunT3BCCoupledWilsonSequential.py` — B/C runner
- `tests/test_t3_bc_coupled_wilson.py` — regression tests
- `docs/T3_BC_COUPLED_WILSON_SEQUENTIAL.md` — this file

No existing project code or previous results are overwritten.

## PowerShell commands

```powershell
python -X utf8 -m pytest -q tests/test_t3_bc_coupled_wilson.py

python -m scripts.RunT3BCCoupledWilsonSequential --model B_-1 --run-id continuous_v1

python -m scripts.RunT3BCCoupledWilsonSequential --all --run-id continuous_v1 --resume
```

If the first run is blocked, its error identifies the parity or source check.
The second command should be used only after that first run is successful.

For an independent run with a fresh output directory:

```powershell
python -m scripts.RunT3BCCoupledWilsonSequential --all --run-id continuous_v1 --output-id bc_coupled_wilson_v2
```

## Outputs

Under:

`output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/bc_coupled_wilson_v1/`

- `B_m1_coupled_wilson.json`
- `C_m1_coupled_wilson.json`
- `summary.csv`, `summary.json`

The actual filenames are `<model_tag>_coupled_wilson.json`, with model tags
as returned by `continuous.model_tag`, e.g. `B_m1_coupled_wilson.json`.

## Equations and interpretation

`C` is the one-generation group-index Wilson tensor with its physical lepton
flavour matrix retained in each group component:

\[
16\pi^2\frac{dC_{ijab}^{pq}}{d\ln\mu}
=\beta_{ijab}^{pq}[C,C^*,g,Y,\lambda].
\]

The seven-term master equation used is **exactly** the project's existing
`RGE.general.WilsonTensorRGE` at the chosen one-generation gauge-index closure.
For the three SM Yukawa matrices, their Frobenius norms project onto the
one-generation master; this matches the Yukawa trace but is not a complete
three-generation flavour mixing calculation.

The old selected-source result was approximately

\[
C_5^{\rm previous}=C_5^{\rm direct}+
\int_{\ln\mu_F}^{\ln\mu_S} \!dt\;
 \beta_{5,\rm extra}[C_{LLSS}(\mu_F),\lambda(t)]/(16\pi^2).
\]

The new calculation transports the **whole L/eC gauge-index coefficient
vector** and projects its Higgs-Higgs sector after running, using the same
physical Weinberg tensor convention.

Uncertified or absent: full three-family charged-lepton Wilson tensor mixing,
finite one-loop scalar-threshold matching, diagonalization/decoupling at each
individual mixed-scalar mass threshold, and a physical certification of all
Stage29 coupling transformations. The output **must not** be called complete
sequential general-scalar matching or a newly fitted neutrino prediction.

The new objective is quantitative: measure whether the earlier 26.4%/33.1%
selected-source corrections remain stable when the intermediate `LLSS` tensor
is evolved, while exposing additional interactions from the declared
conditional scalar basis.

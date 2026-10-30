# A(-2): anomalous one-scale C5 from the general scalar quartic

## Scope

This **new verification stage** is independent of production matching, RGBeta generation, and numerical fitting. It depends on the successful prior A(-2) bilinear normalization certificate:

`output/studies/general_scalar_extended/a2_bilinear_normalization/a2_norm_v1/certificate.json`

It checks the **existing** Matchete SU(2) vertex tensors, verifies the exact absolute one-scale loop coefficient against a separately integrated Euclidean loop, and constructs a conditional three-generation anomalous Weinberg coefficient. It **does not** establish an explicit physical Matchete–RGBeta spherical/Cartesian field map or full one-loop sequential Wilson matching. Do not use this as a production neutrino prediction.

## Physics

For A(-2), the general scalar basis includes the additional independent complex quartic

\[
Q= (H^\dagger H^\dagger S_1^\dagger S_2^\dagger)_{\boldsymbol1}.
\]

Repository identifiers:

- Matchete `lambdaT3ExtraQ020101I1`
- RGBeta `lambdaT3GeneralQ020101I1`

In the **same component convention**, Matchete constructs `+lambda_M Q+h.c.` in the Lagrangian, while RGBeta uses `-lambda_R Q+h.c.`; hence `lambda_M=-lambda_R`. The HH channel uses the **Hermitian conjugate** of Q and therefore `conjugate(lambda_M)`.

With scalar potential

\[
V_2=m_{S_2}^2S_2^\dagger S_2 +\big[\mu_R^2(S_2^\dagger)^TC(S_2^\dagger)+\mathrm{h.c.}\big],
\]

its two real-triplet mass eigenvalues are \(m_\pm^2=m_{S_2}^2\pm 2|\mu_R^2|\). The anomalous scalar propagator contains the phase \(\mu_R^2/|\mu_R^2|\), and the normal propagator is their average.

Let \(K_\pm\) denote the **three-generation** C5 from the original Matchete one-scale kernel with \(\lambda_{T3}=1\), evaluated at the mass \(m_\pm\). The native Matchete tensor test checks that the extra quartic's HH-channel tensor, contracted with the self-bilinear metric, equals the original topology tensor with ratio exactly +1. The resulting *conditional* one-scale expression is

\[
C_5^{N}=\frac{\lambda_{T3}}2(K_++K_-),\qquad
\Delta C_5^{A}=-\frac{\lambda_R^*}2\frac{\mu_R^2}{|\mu_R^2|}(K_+-K_-).
\]

**Important restriction:** The second equality assumes the remaining RGBeta–Matchete triplet component-field identification does not introduce an additional factor. The current stage explicitly leaves this unverified.

The stripped Matchete A(-2) one-scale loop kernel has the absolute form

\[
K(M,m_1,m_2)=-\frac2{\sqrt3}\,\frac{1}{16\pi^2}\, M
\int_0^\infty\!\frac{t\,dt}{(t+M^2)(t+m_1^2)(t+m_2^2)}.
\]

No normalization is fit from the integral: the algebraic \(-2/\sqrt3\) factor is read from the known Matchete expression, and the integral is evaluated numerically, independently, for all three benchmark heavy masses and both split scalar masses. The usual symmetric flavour lift is \((Y_1^*\operatorname{diag}(K_r)Y_2^{*T}+Y_2^*\operatorname{diag}(K_r)Y_1^{*T})\), with this source's coefficient convention.

## Install

Extract the supplied ZIP at the **repository root**. It contains only new files:

- `scripts/VerifyA2ExtraQuarticC5.py`
- `tests/test_a2_extra_quartic_c5.py`
- `docs/T3_A_MINUS2_EXTRA_QUARTIC_C5.md`

No existing sources/benchmarks change.

## Execute in PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Aminus2_extra_quartic_C5.zip' -DestinationPath . -Force
python -X utf8 -m pytest -q tests/test_a2_extra_quartic_c5.py
python -u -X utf8 -m scripts.VerifyA2ExtraQuarticC5 --normalization-run-id a2_norm_v1 --run-id a2_extra_c5_v1
```

Default illustrative extra RGBeta coupling is **0.02+0.01i**, *not a fitted value*. To choose a different one:

```powershell
python -u -X utf8 -m scripts.VerifyA2ExtraQuarticC5 --normalization-run-id a2_norm_v1 --run-id a2_extra_c5_v2 --extra-rgbeta-real 0.01 --extra-rgbeta-imag 0.00
```

Report: `output/studies/general_scalar_extended/a2_extra_quartic_c5/a2_extra_c5_v1/certificate.json`.

## Interpreting the result

`A2NativeTensorAndConditionalOneScaleAnomalousC5Checked` means:

1. Native Matchete tensor proportionality with the topology is checked from **source tensors**, not inferred from operator names.
2. The relative Matchete/RGBeta potential sign and HH-channel Hermitian conjugation are sourced and applied **conditionally** on a common triplet component convention.
3. All three heavy masses receive the **absolute** numerical loop-integral check for both split masses.
4. The numerical three-generation C5 normal/anomalous matrices agree with independent integral evaluation.

It does **not** certify the underlying RGBeta–Matchete field-component intertwiner, UV/EFT matching thresholds, higher-order diagrams, or a full general-scalar neutrino-mass calculation.

Next needed for physical certification: explicitly compare native RGBeta and Matchete **component tensors under a single unitary SU(2) triplet intertwiner**, including Yukawa contractions, followed by mass-eigenstate Matchete matching and scheme-consistent EFT matching.

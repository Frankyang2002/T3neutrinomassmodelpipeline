# T3 E_-2: unified general-model physics verification

S1 = (3,-1), S2 = (3,0), F = (2,-1/2).

**Source-matched structural inventory:** VerifiedStructural.  **Gauge/group coefficients verified within RGBeta:** 100/100.

## Complete general BSM Lagrangian and UV/EFT1 RGE expressions

[Full model report](../../GeneralSymbolic/models/E_m2.md) displays the native Matchete BSM expression and all exported RGBeta UV/EFT1 renormalisable beta functions.

These are candidate RGBeta RGEs in the extended interaction basis; the physical cross-package coupling dictionary is not yet certified.

## Gauge group factors

The tables independently project the coefficients of `g^3` (gauge), `g^2*lambda` (scalar), and `g^2*Matrix[y]` (Yukawa). These checks **do not** test all non-gauge scalar recoupling or matrix-product terms.

| Stage | Beta | Gauge | Expected | Extracted | Status |
|---|---|---|---:|---:|---|
| UV | `gY` | gY | 59/6 | 59/6 | VerifiedRGBetaGaugeCoefficient |
| UV | `g2` | g2 | 1/6 | 1/6 | VerifiedRGBetaGaugeCoefficient |
| UV | `g3` | g3 | -7 | -7 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaH` | g2 | -9 | -9 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaH` | gY | -3 | -3 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaS1` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaS1` | gY | -12 | -12 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaS2` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaS2` | gY | 0 | 0 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaH1` | g2 | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaH1` | gY | -15/2 | -15/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaH2` | g2 | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaH2` | gY | -3/2 | -3/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambda12` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambda12` | gY | -6 | -6 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaT3` | g2 | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaT3` | gY | -9/2 | -9/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `mS1Sq` | g2 | -12 | -12 | VerifiedRGBetaGaugeCoefficient |
| UV | `mS1Sq` | gY | -6 | -6 | VerifiedRGBetaGaugeCoefficient |
| UV | `mS2Sq` | g2 | -12 | -12 | VerifiedRGBetaGaugeCoefficient |
| UV | `mS2Sq` | gY | 0 | 0 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaT3GeneralQ000004I1` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaT3GeneralQ000004I1` | gY | 0 | 0 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaT3GeneralQ000013I1` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaT3GeneralQ000013I1` | gY | 0 | 0 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaT3GeneralQ001102I1` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaT3GeneralQ001102I1` | gY | -6 | -6 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaT3GeneralQ001102I2` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaT3GeneralQ001102I2` | gY | -6 | -6 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaT3GeneralQ020101I1` | g2 | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaT3GeneralQ020101I1` | gY | -9/2 | -9/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaT3GeneralQ110002I1` | g2 | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaT3GeneralQ110002I1` | gY | -3/2 | -3/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaH1Adj` | g2 | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaH1Adj` | gY | -15/2 | -15/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaH2Adj` | g2 | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaH2Adj` | gY | -3/2 | -3/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaS1Adj` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaS1Adj` | gY | -12 | -12 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaS2Adj` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaS2Adj` | gY | 0 | 0 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambda12Adj` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambda12Adj` | gY | -6 | -6 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambda12Cross` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambda12Cross` | gY | -6 | -6 | VerifiedRGBetaGaugeCoefficient |
| UV | `muT3GeneralS2Sq` | g2 | -12 | -12 | VerifiedRGBetaGaugeCoefficient |
| UV | `muT3GeneralS2Sq` | gY | 0 | 0 | VerifiedRGBetaGaugeCoefficient |
| UV | `y1` | g2 | -9/2 | -9/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `y1` | gY | -3/2 | -3/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `y2` | g2 | -9/2 | -9/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `y2` | gY | -3/2 | -3/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `y2Conj` | g2 | -9/2 | -9/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `y2Conj` | gY | -3/2 | -3/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `gY` | gY | 47/6 | 47/6 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `g2` | g2 | -11/6 | -11/6 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `g3` | g3 | -7 | -7 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaH` | g2 | -9 | -9 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaH` | gY | -3 | -3 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaS1` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaS1` | gY | -12 | -12 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaS2` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaS2` | gY | 0 | 0 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaH1` | g2 | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaH1` | gY | -15/2 | -15/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaH2` | g2 | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaH2` | gY | -3/2 | -3/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambda12` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambda12` | gY | -6 | -6 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaT3` | g2 | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaT3` | gY | -9/2 | -9/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `mS1Sq` | g2 | -12 | -12 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `mS1Sq` | gY | -6 | -6 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `mS2Sq` | g2 | -12 | -12 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `mS2Sq` | gY | 0 | 0 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaT3GeneralQ000004I1` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaT3GeneralQ000004I1` | gY | 0 | 0 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaT3GeneralQ000013I1` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaT3GeneralQ000013I1` | gY | 0 | 0 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaT3GeneralQ001102I1` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaT3GeneralQ001102I1` | gY | -6 | -6 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaT3GeneralQ001102I2` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaT3GeneralQ001102I2` | gY | -6 | -6 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaT3GeneralQ020101I1` | g2 | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaT3GeneralQ020101I1` | gY | -9/2 | -9/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaT3GeneralQ110002I1` | g2 | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaT3GeneralQ110002I1` | gY | -3/2 | -3/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaH1Adj` | g2 | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaH1Adj` | gY | -15/2 | -15/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaH2Adj` | g2 | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaH2Adj` | gY | -3/2 | -3/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaS1Adj` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaS1Adj` | gY | -12 | -12 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaS2Adj` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaS2Adj` | gY | 0 | 0 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambda12Adj` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambda12Adj` | gY | -6 | -6 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambda12Cross` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambda12Cross` | gY | -6 | -6 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `muT3GeneralS2Sq` | g2 | -12 | -12 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `muT3GeneralS2Sq` | gY | 0 | 0 | VerifiedRGBetaGaugeCoefficient |

## Physical coupling correspondence

The coupling ledger records the scalar-invariant basis block and additional Yukawa vertices. The relative Matchete/RGBeta action sign, conjugation, phase and normalisation remain **NotCertified**.

| Sector | Native Matchete coupling | RGBeta candidate(s) | Tensor rank |
|---|---|---|---:|
| extra_quartic | `lambdaT3ExtraQ000004I1` | `lambdaT3GeneralQ000004I1` | 1 |
| extra_quartic | `lambdaT3ExtraQ000013I1` | `lambdaT3GeneralQ000013I1` | 1 |
| extra_quartic | `lambdaT3ExtraQ001102I1` | `lambdaT3GeneralQ001102I1;lambdaT3GeneralQ001102I2` | 2 |
| extra_quartic | `lambdaT3ExtraQ001102I2` | `lambdaT3GeneralQ001102I1;lambdaT3GeneralQ001102I2` | 2 |
| extra_quartic | `lambdaT3ExtraQ020101I1` | `lambdaT3GeneralQ020101I1` | 1 |
| extra_quartic | `lambdaT3ExtraQ110002I1` | `lambdaT3GeneralQ110002I1` | 1 |
| extra_bilinear | `lambdaT3ExtraM000002I1` | `muT3GeneralS2Sq` | 1 |
| yukawa | `y1` | `y1` | 1 |
| yukawa | `y2` | `y2` | 1 |
| yukawa | `y2Conj` | `y2Conj` | 1 |

## EFT scope and open physics

- UV and EFT1 **renormalisable** RGEs are exported; EFT1 dimension-five Wilson-operator evolution is not rederived here.
- Full extended matching and final C5 are not certified. B(-1)/C(-1) require scalar mass diagonalisation.
- The separate legacy restricted-model matching must not be labelled as matching the general BSM Lagrangian.

# T3 A_-4: unified general-model physics verification

S1 = (1,-2), S2 = (3,-1), F = (2,-3/2).

**Source-matched structural inventory:** VerifiedStructural.  **Gauge/group coefficients verified within RGBeta:** 54/54.

## Complete general BSM Lagrangian and UV/EFT1 RGE expressions

[Full model report](../../GeneralSymbolic/models/A_m4.md) displays the native Matchete BSM expression and all exported RGBeta UV/EFT1 renormalisable beta functions.

These are candidate RGBeta RGEs in the extended interaction basis; the physical cross-package coupling dictionary is not yet certified.

## Gauge group factors

The tables independently project the coefficients of `g^3` (gauge), `g^2*lambda` (scalar), and `g^2*Matrix[y]` (Yukawa). These checks **do not** test all non-gauge scalar recoupling or matrix-product terms.

| Stage | Beta | Gauge | Expected | Extracted | Status |
|---|---|---|---:|---:|---|
| UV | `gY` | gY | 163/6 | 163/6 | VerifiedRGBetaGaugeCoefficient |
| UV | `g2` | g2 | -1/2 | -1/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `g3` | g3 | -7 | -7 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaH` | g2 | -9 | -9 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaH` | gY | -3 | -3 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaS1` | g2 | 0 | 0 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaS1` | gY | -48 | -48 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaS2` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaS2` | gY | -12 | -12 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaH1` | g2 | -9/2 | -9/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaH1` | gY | -51/2 | -51/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaH2` | g2 | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaH2` | gY | -15/2 | -15/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambda12` | g2 | -12 | -12 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambda12` | gY | -30 | -30 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaT3` | g2 | -21/2 | -21/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaT3` | gY | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `mS1Sq` | g2 | 0 | 0 | VerifiedRGBetaGaugeCoefficient |
| UV | `mS1Sq` | gY | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| UV | `mS2Sq` | g2 | -12 | -12 | VerifiedRGBetaGaugeCoefficient |
| UV | `mS2Sq` | gY | -6 | -6 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaH2Adj` | g2 | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaH2Adj` | gY | -15/2 | -15/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaS2Adj` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| UV | `lambdaS2Adj` | gY | -12 | -12 | VerifiedRGBetaGaugeCoefficient |
| UV | `y1` | g2 | -9/2 | -9/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `y1` | gY | -15/2 | -15/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `y2` | g2 | -9/2 | -9/2 | VerifiedRGBetaGaugeCoefficient |
| UV | `y2` | gY | -15/2 | -15/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `gY` | gY | 55/6 | 55/6 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `g2` | g2 | -5/2 | -5/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `g3` | g3 | -7 | -7 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaH` | g2 | -9 | -9 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaH` | gY | -3 | -3 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaS1` | g2 | 0 | 0 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaS1` | gY | -48 | -48 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaS2` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaS2` | gY | -12 | -12 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaH1` | g2 | -9/2 | -9/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaH1` | gY | -51/2 | -51/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaH2` | g2 | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaH2` | gY | -15/2 | -15/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambda12` | g2 | -12 | -12 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambda12` | gY | -30 | -30 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaT3` | g2 | -21/2 | -21/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaT3` | gY | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `mS1Sq` | g2 | 0 | 0 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `mS1Sq` | gY | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `mS2Sq` | g2 | -12 | -12 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `mS2Sq` | gY | -6 | -6 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaH2Adj` | g2 | -33/2 | -33/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaH2Adj` | gY | -15/2 | -15/2 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaS2Adj` | g2 | -24 | -24 | VerifiedRGBetaGaugeCoefficient |
| EFT1 | `lambdaS2Adj` | gY | -12 | -12 | VerifiedRGBetaGaugeCoefficient |

## Physical coupling correspondence

The coupling ledger records the scalar-invariant basis block and additional Yukawa vertices. The relative Matchete/RGBeta action sign, conjugation, phase and normalisation remain **NotCertified**.

| Sector | Native Matchete coupling | RGBeta candidate(s) | Tensor rank |
|---|---|---|---:|
| yukawa | `y1` | `y1` | 1 |
| yukawa | `y2` | `y2` | 1 |

## EFT scope and open physics

- UV and EFT1 **renormalisable** RGEs are exported; EFT1 dimension-five Wilson-operator evolution is not rederived here.
- Full extended matching and final C5 are not certified. B(-1)/C(-1) require scalar mass diagonalisation.
- The separate legacy restricted-model matching must not be labelled as matching the general BSM Lagrangian.

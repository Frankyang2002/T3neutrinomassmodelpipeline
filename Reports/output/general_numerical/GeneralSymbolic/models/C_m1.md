# T3 C_-1: general-interaction symbolic results

Fields: S1=(2,-1/2), S2=(2,1/2), F=(3,0) [SU(2) dimension, hypercharge].

Additional independent quartics: **7**; additional complex bilinears: **1**; additional Yukawa: **none**.

Heavy scalar mass basis: **Blocked**. General sequential matching: **not performed**. Matchete/RGBeta physical coupling dictionary: **not certified**.

## Independent interaction inventory

- `Yukawa1_GenericCG`
- `Yukawa2_GenericCG`
- `Scalar1SelfInv1`
- `Scalar2SelfInv1`
- `HiggsPortal1Inv1`
- `HiggsPortal1Inv2`
- `HiggsPortal2Inv1`
- `HiggsPortal2Inv2`
- `ScalarCrossPortalInv1`
- `ScalarCrossPortalInv2`
- `T3ScalarMix_GenericCG`
- `T3ExtraQ_0_0_0_1_1_2_1`
- `T3ExtraQ_0_0_0_2_0_2_1`
- `T3ExtraQ_0_0_1_2_0_1_1`
- `T3ExtraQ_0_2_0_0_2_0_1`
- `T3ExtraQ_0_2_0_2_0_0_1`
- `T3ExtraQ_1_1_0_1_0_1_1`
- `T3ExtraQ_1_1_0_1_0_1_2`
- `T3ExtraM_0_0_0_1_0_1_1`

## Complete BSM UV Lagrangian (native Matchete InputForm)

The expression below includes BSM free fields and interactions; it excludes the pure Standard Model Lagrangian.

```text
Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd2]}, {}]]*CG[T3ExtraCGM000101I1, {Index[d$$1, T3BSMd2], Index[d$$2, T3BSMd2]}]*Coupling[lambdaT3ExtraM000101I1, {}, 0] + Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*Bar[Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd2]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$3, T3BSMd2]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd2]}, {}]]*CG[T3ExtraCGQ000202I1, {Index[d$$1, T3BSMd2], Index[d$$2, T3BSMd2], Index[d$$3, T3BSMd2], Index[d$$4, T3BSMd2]}]*Coupling[lambdaT3ExtraQ000202I1, {}, 0] + Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*Bar[Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd2]}, {}]]*CG[T3ExtraCGQ020200I1, {Index[d$$1, SU2L[fund]], Index[d$$2, SU2L[fund]], Index[d$$1, T3BSMd2], Index[d$$2, T3BSMd2]}]*Coupling[lambdaT3ExtraQ020200I1, {}, 0] + Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd2]}, {}]]*CG[T3ExtraCGQ110101I1, {Bar[Index[d$$2, SU2L[fund]]], Index[d$$1, SU2L[fund]], Index[d$$1, T3BSMd2], Index[d$$2, T3BSMd2]}]*Coupling[lambdaT3ExtraQ110101I1, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}] + Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd2]}, {}]]*CG[T3ExtraCGQ110101I2, {Bar[Index[d$$2, SU2L[fund]]], Index[d$$1, SU2L[fund]], Index[d$$1, T3BSMd2], Index[d$$2, T3BSMd2]}]*Coupling[lambdaT3ExtraQ110101I2, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}] + Bar[Coupling[lambdaT3ExtraQ020020I1, {}, 0]]*Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd2]}, {}]]*CG[Bar[T3ExtraCGQ020020I1], {Bar[Index[d$$1, SU2L[fund]]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$1, T3BSMd2], Index[d$$2, T3BSMd2]}]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}] - Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*Coupling[MS1, {}, 0]^2*Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}] + Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {Index[d$$1, Lorentz]}]]*Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {Index[d$$1, Lorentz]}] + Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*CG[T3HiggsPortal1CGInv1, {Index[d$$1, SU2L[fund]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$1, T3BSMd2], Bar[Index[d$$2, T3BSMd2]]}]*Coupling[lambdaH1Inv1, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd2]}, {}] + Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*CG[T3HiggsPortal1CGInv2, {Index[d$$1, SU2L[fund]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$1, T3BSMd2], Bar[Index[d$$2, T3BSMd2]]}]*Coupling[lambdaH1Inv2, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd2]}, {}] + Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*CG[T3MixCG, {Bar[Index[d$$1, SU2L[fund]]], Bar[Index[d$$2, SU2L[fund]]], Bar[Index[d$$2, T3BSMd2]], Index[d$$1, T3BSMd2]}]*Coupling[lambdaT3, {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd2]}, {}] + Bar[Coupling[lambdaT3ExtraQ020200I1, {}, 0]]*CG[Bar[T3ExtraCGQ020200I1], {Bar[Index[d$$1, SU2L[fund]]], Bar[Index[d$$2, SU2L[fund]]], Bar[Index[d$$1, T3BSMd2]], Bar[Index[d$$2, T3BSMd2]]}]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]*Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd2]}, {}] + Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*Bar[Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd2]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$3, T3BSMd2]}, {}]]*CG[T3ExtraCGQ001201I1, {Bar[Index[d$$4, T3BSMd2]], Index[d$$1, T3BSMd2], Index[d$$2, T3BSMd2], Index[d$$3, T3BSMd2]}]*Coupling[lambdaT3ExtraQ001201I1, {}, 0]*Field[NewScalar1, Scalar, {Index[d$$4, T3BSMd2]}, {}] + (Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*Bar[Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd2]}, {}]]*CG[T3Scalar1SelfCGInv1, {Index[d$$1, T3BSMd2], Bar[Index[d$$3, T3BSMd2]], Index[d$$2, T3BSMd2], Bar[Index[d$$4, T3BSMd2]]}]*Coupling[lambdaS1Inv1, {}, 0]*Field[NewScalar1, Scalar, {Index[d$$3, T3BSMd2]}, {}]*Field[NewScalar1, Scalar, {Index[d$$4, T3BSMd2]}, {}])/2 - Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*Coupling[MS2, {}, 0]^2*Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd2]}, {}] + Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd2]}, {Index[d$$1, Lorentz]}]]*Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd2]}, {Index[d$$1, Lorentz]}] + Bar[Coupling[lambdaT3, {}, 0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*CG[Bar[T3MixCG], {Index[d$$1, SU2L[fund]], Index[d$$2, SU2L[fund]], Index[d$$1, T3BSMd2], Bar[Index[d$$2, T3BSMd2]]}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd2]}, {}] + Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*CG[T3HiggsPortal2CGInv1, {Index[d$$1, SU2L[fund]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$1, T3BSMd2], Bar[Index[d$$2, T3BSMd2]]}]*Coupling[lambdaH2Inv1, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd2]}, {}] + Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*CG[T3HiggsPortal2CGInv2, {Index[d$$1, SU2L[fund]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$1, T3BSMd2], Bar[Index[d$$2, T3BSMd2]]}]*Coupling[lambdaH2Inv2, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd2]}, {}] + Bar[Coupling[lambdaT3ExtraM000101I1, {}, 0]]*CG[Bar[T3ExtraCGM000101I1], {Bar[Index[d$$1, T3BSMd2]], Bar[Index[d$$2, T3BSMd2]]}]*Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd2]}, {}] + Bar[Coupling[lambdaT3ExtraQ110101I1, {}, 0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*CG[Bar[T3ExtraCGQ110101I1], {Index[d$$1, SU2L[fund]], Bar[Index[d$$2, SU2L[fund]]], Bar[Index[d$$1, T3BSMd2]], Bar[Index[d$$2, T3BSMd2]]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd2]}, {}] + Bar[Coupling[lambdaT3ExtraQ110101I2, {}, 0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*CG[Bar[T3ExtraCGQ110101I2], {Index[d$$1, SU2L[fund]], Bar[Index[d$$2, SU2L[fund]]], Bar[Index[d$$1, T3BSMd2]], Bar[Index[d$$2, T3BSMd2]]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd2]}, {}] + Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*CG[T3ExtraCGQ020020I1, {Index[d$$1, SU2L[fund]], Index[d$$2, SU2L[fund]], Bar[Index[d$$1, T3BSMd2]], Bar[Index[d$$2, T3BSMd2]]}]*Coupling[lambdaT3ExtraQ020020I1, {}, 0]*Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd2]}, {}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd2]}, {}] + Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd2]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$3, T3BSMd2]}, {}]]*CG[T3ExtraCGQ000112I1, {Index[d$$1, T3BSMd2], Bar[Index[d$$4, T3BSMd2]], Index[d$$2, T3BSMd2], Index[d$$3, T3BSMd2]}]*Coupling[lambdaT3ExtraQ000112I1, {}, 0]*Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd2]}, {}] + Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd2]}, {}]]*CG[T3CrossScalarCGInv1, {Index[d$$1, T3BSMd2], Bar[Index[d$$3, T3BSMd2]], Index[d$$2, T3BSMd2], Bar[Index[d$$4, T3BSMd2]]}]*Coupling[lambda12Inv1, {}, 0]*Field[NewScalar1, Scalar, {Index[d$$3, T3BSMd2]}, {}]*Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd2]}, {}] + Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd2]}, {}]]*CG[T3CrossScalarCGInv2, {Index[d$$1, T3BSMd2], Bar[Index[d$$3, T3BSMd2]], Index[d$$2, T3BSMd2], Bar[Index[d$$4, T3BSMd2]]}]*Coupling[lambda12Inv2, {}, 0]*Field[NewScalar1, Scalar, {Index[d$$3, T3BSMd2]}, {}]*Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd2]}, {}] + Bar[Coupling[lambdaT3ExtraQ001201I1, {}, 0]]*Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*CG[Bar[T3ExtraCGQ001201I1], {Index[d$$1, T3BSMd2], Bar[Index[d$$2, T3BSMd2]], Bar[Index[d$$3, T3BSMd2]], Bar[Index[d$$4, T3BSMd2]]}]*Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd2]}, {}]*Field[NewScalar1, Scalar, {Index[d$$3, T3BSMd2]}, {}]*Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd2]}, {}] + (Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd2]}, {}]]*CG[T3Scalar2SelfCGInv1, {Index[d$$1, T3BSMd2], Bar[Index[d$$3, T3BSMd2]], Index[d$$2, T3BSMd2], Bar[Index[d$$4, T3BSMd2]]}]*Coupling[lambdaS2Inv1, {}, 0]*Field[NewScalar2, Scalar, {Index[d$$3, T3BSMd2]}, {}]*Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd2]}, {}])/2 + Bar[Coupling[lambdaT3ExtraQ000112I1, {}, 0]]*Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*CG[Bar[T3ExtraCGQ000112I1], {Bar[Index[d$$2, T3BSMd2]], Index[d$$1, T3BSMd2], Bar[Index[d$$3, T3BSMd2]], Bar[Index[d$$4, T3BSMd2]]}]*Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd2]}, {}]*Field[NewScalar2, Scalar, {Index[d$$3, T3BSMd2]}, {}]*Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd2]}, {}] + Bar[Coupling[lambdaT3ExtraQ000202I1, {}, 0]]*CG[Bar[T3ExtraCGQ000202I1], {Bar[Index[d$$1, T3BSMd2]], Bar[Index[d$$2, T3BSMd2]], Bar[Index[d$$3, T3BSMd2]], Bar[Index[d$$4, T3BSMd2]]}]*Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]*Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd2]}, {}]*Field[NewScalar2, Scalar, {Index[d$$3, T3BSMd2]}, {}]*Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd2]}, {}] + Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*CG[T3Y2CG, {Index[d$$1, SU2L[fund]], Index[d$$1, T3BSMd3], Index[d$$1, T3BSMd2]}]*Coupling[y2, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]*NCM[Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]], DiracProduct[Proj[1]], Field[NewFermion, Fermion, {Index[d$$1, T3BSMd3], Index[d$$1, NFlavor]}, {}]] + CG[T3Y1CG, {Index[d$$1, SU2L[fund]], Index[d$$1, T3BSMd3], Bar[Index[d$$1, T3BSMd2]]}]*Coupling[y1, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]*Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]*NCM[Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]], DiracProduct[Proj[1]], Field[NewFermion, Fermion, {Index[d$$1, T3BSMd3], Index[d$$1, NFlavor]}, {}]] - (Coupling[MF, {}, 0]*NCM[Transp[Field[NewFermion, Fermion, {Index[d$$1, T3BSMd3], Index[d$$1, NFlavor]}, {}]], DiracProduct[GammaCC], Field[NewFermion, Fermion, {Index[d$$1, T3BSMd3], Index[d$$1, NFlavor]}, {}]])/2 + (I/2)*NCM[Transp[Field[NewFermion, Fermion, {Index[d$$1, T3BSMd3], Index[d$$1, NFlavor]}, {}]], DiracProduct[GammaCC, GammaM[Index[d$$1, Lorentz]]], Field[NewFermion, Fermion, {Index[d$$1, T3BSMd3], Index[d$$1, NFlavor]}, {Index[d$$1, Lorentz]}]] + Bar[Coupling[y1, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]]*Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd2]}, {}]]*CG[Bar[T3Y1CG], {Bar[Index[d$$1, SU2L[fund]]], Index[d$$1, T3BSMd3], Index[d$$1, T3BSMd2]}]*NCM[Transp[Field[NewFermion, Fermion, {Index[d$$1, T3BSMd3], Index[d$$1, NFlavor]}, {}]], DiracProduct[GammaCC, Proj[-1]], Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]] + Bar[Coupling[y2, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]]*CG[Bar[T3Y2CG], {Bar[Index[d$$1, SU2L[fund]]], Index[d$$1, T3BSMd3], Bar[Index[d$$1, T3BSMd2]]}]*Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd2]}, {}]*NCM[Transp[Field[NewFermion, Fermion, {Index[d$$1, T3BSMd3], Index[d$$1, NFlavor]}, {}]], DiracProduct[GammaCC, Proj[-1]], Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]]
```

## One-loop renormalisable beta functions

RGBeta's report convention is `16*pi^2*dX/dln(mu)`. These are candidate native RGBeta expressions, not independently mapped/certified Matchete coefficients.

### UV: 29 beta functions

**beta_gY**

```text
(43*gY^3)/6
```

**beta_g2**

```text
(7*g2^3)/6
```

**beta_g3**

```text
-7*g3^3
```

**beta_yu**

```text
(-9*g2^2*Matrix[yu][gen[$i], gen[$j]])/4 - 8*g3^2*Matrix[yu][gen[$i], gen[$j]] - (17*gY^2*Matrix[yu][gen[$i], gen[$j]])/12 + 3*Tr[yd . Trans[Bar[yd]]]*Matrix[yu][gen[$i], gen[$j]] + Tr[ye . Trans[Bar[ye]]]*Matrix[yu][gen[$i], gen[$j]] + 3*Tr[yu . Trans[Bar[yu]]]*Matrix[yu][gen[$i], gen[$j]] - (3*(Matrix[yd, Trans[Bar[yd]], yu][gen[$i], gen[$j]] - Matrix[yu, Trans[Bar[yu]], yu][gen[$i], gen[$j]]))/2
```

**beta_yd**

```text
(-9*g2^2*Matrix[yd][gen[$i], gen[$j]])/4 - 8*g3^2*Matrix[yd][gen[$i], gen[$j]] - (5*gY^2*Matrix[yd][gen[$i], gen[$j]])/12 + 3*Tr[yd . Trans[Bar[yd]]]*Matrix[yd][gen[$i], gen[$j]] + Tr[ye . Trans[Bar[ye]]]*Matrix[yd][gen[$i], gen[$j]] + 3*Tr[yu . Trans[Bar[yu]]]*Matrix[yd][gen[$i], gen[$j]] + (3*(Matrix[yd, Trans[Bar[yd]], yd][gen[$i], gen[$j]] - Matrix[yu, Trans[Bar[yu]], yd][gen[$i], gen[$j]]))/2
```

**beta_ye**

```text
(-9*g2^2*Matrix[ye][gen[$i], gen[$j]])/4 - (15*gY^2*Matrix[ye][gen[$i], gen[$j]])/4 + 3*Tr[yd . Trans[Bar[yd]]]*Matrix[ye][gen[$i], gen[$j]] + Tr[ye . Trans[Bar[ye]]]*Matrix[ye][gen[$i], gen[$j]] + 3*Tr[yu . Trans[Bar[yu]]]*Matrix[ye][gen[$i], gen[$j]] + (3*(Matrix[y1, Trans[Bar[y1]], ye][gen[$i], gen[$j]] + Matrix[y2, Trans[Bar[y2]], ye][gen[$i], gen[$j]] + 2*Matrix[ye, Trans[Bar[ye]], ye][gen[$i], gen[$j]]))/4
```

**beta_y1**

```text
(-33*g2^2*Matrix[y1][gen[$i], heavy[$j]])/4 - (3*gY^2*Matrix[y1][gen[$i], heavy[$j]])/4 + (3*Tr[y1 . Trans[Bar[y1]]]*Matrix[y1][gen[$i], heavy[$j]])/2 + (3*Tr[y1 . Trans[Bar[y2]]]*Matrix[y2][gen[$i], heavy[$j]])/2 + (5*Matrix[y1, Trans[Bar[y1]], y1][gen[$i], heavy[$j]] + 2*Matrix[y1, Trans[Bar[y2]], y2][gen[$i], heavy[$j]] + 3*Matrix[y2, Trans[Bar[y2]], y1][gen[$i], heavy[$j]] + 2*Matrix[ye, Trans[Bar[ye]], y1][gen[$i], heavy[$j]])/4
```

**beta_y2**

```text
(3*Tr[y2 . Trans[Bar[y1]]]*Matrix[y1][gen[$i], heavy[$j]])/2 - (33*g2^2*Matrix[y2][gen[$i], heavy[$j]])/4 - (3*gY^2*Matrix[y2][gen[$i], heavy[$j]])/4 + (3*Tr[y2 . Trans[Bar[y2]]]*Matrix[y2][gen[$i], heavy[$j]])/2 + (3*Matrix[y1, Trans[Bar[y1]], y2][gen[$i], heavy[$j]] + 2*Matrix[y2, Trans[Bar[y1]], y1][gen[$i], heavy[$j]] + 5*Matrix[y2, Trans[Bar[y2]], y2][gen[$i], heavy[$j]] + 2*Matrix[ye, Trans[Bar[ye]], y2][gen[$i], heavy[$j]])/4
```

**beta_MF**

```text
-6*g2^2*(Matrix[MF][heavy[$i], heavy[$j]] + Matrix[Trans[MF]][heavy[$i], heavy[$j]]) + (Matrix[Trans[MF], Trans[Bar[y1]], y1][heavy[$i], heavy[$j]] + Matrix[Trans[MF], Trans[Bar[y2]], y2][heavy[$i], heavy[$j]] + Matrix[Trans[y1], Bar[y1], MF][heavy[$i], heavy[$j]] + Matrix[Trans[y2], Bar[y2], MF][heavy[$i], heavy[$j]])/2
```

**beta_mS1Sq**

```text
(-9*g2^2*mS1Sq)/2 - (3*gY^2*mS1Sq)/2 + 4*lambda12*mS2Sq + 6*muT3GeneralS1S2*Bar[lambdaS1S1bar2S2bar] + 6*lambdaS1S1bar2S2bar*Bar[muT3GeneralS1S2] + mS1Sq*(6*lambdaS1 + 3*Tr[y1 . Trans[Bar[y1]]]) - 12*Tr[MF . Bar[MF] . Trans[y1] . Bar[y1]] - 12*Tr[MF . Trans[Bar[y1]] . y1 . Bar[MF]]
```

**beta_mS2Sq**

```text
4*lambda12*mS1Sq - (9*g2^2*mS2Sq)/2 - (3*gY^2*mS2Sq)/2 + 6*muT3GeneralS1S2*Bar[lambdaS1barS2S2bar2] + 6*lambdaS1barS2S2bar2*Bar[muT3GeneralS1S2] + mS2Sq*(6*lambdaS2 + 3*Tr[y2 . Trans[Bar[y2]]]) - 12*Tr[MF . Bar[MF] . Trans[y2] . Bar[y2]] - 12*Tr[MF . Trans[Bar[y2]] . y2 . Bar[MF]]
```

**beta_lambdaH**

```text
(9*g2^4)/4 + (3*gY^4)/4 + gY^2*((3*g2^2)/2 - 3*lambdaH) - 9*g2^2*lambdaH + 12*lambdaH^2 + 4*lambdaH1^2 + lambdaH1Adj^2/4 + 4*lambdaH2^2 + lambdaH2Adj^2/4 + 8*lambdaHHdagS1barS1bar*Bar[lambdaHHdagS1barS1bar] + 4*lambdaHHdagS1barS2barCross*Bar[lambdaHHdagS1barS2barCross] + 8*lambdaHHdagS2S2*Bar[lambdaHHdagS2S2] + 4*lambdaT3*Bar[lambdaT3] + 4*lambdaT3GeneralQ110101I1*Bar[lambdaT3GeneralQ110101I1] + lambdaH*(12*Tr[yd . Trans[Bar[yd]]] + 4*Tr[ye . Trans[Bar[ye]]] + 12*Tr[yu . Trans[Bar[yu]]]) - 12*Tr[yd . Trans[Bar[yd]] . yd . Trans[Bar[yd]]] - 4*Tr[ye . Trans[Bar[ye]] . ye . Trans[Bar[ye]]] - 12*Tr[yu . Trans[Bar[yu]] . yu . Trans[Bar[yu]]]
```

**beta_lambdaS1**

```text
(9*g2^4)/4 + (3*gY^4)/4 + 4*lambda12^2 + lambda12Adj^2/4 + 4*lambdaH1^2 + lambdaH1Adj^2/4 + gY^2*((3*g2^2)/2 - 3*lambdaS1) - 9*g2^2*lambdaS1 + 12*lambdaS1^2 + 8*lambdaHHdagS1barS1bar*Bar[lambdaHHdagS1barS1bar] + 8*lambdaS1bar2S2bar2*Bar[lambdaS1bar2S2bar2] + 24*lambdaS1S1bar2S2bar*Bar[lambdaS1S1bar2S2bar] + 6*lambdaS1*Tr[y1 . Trans[Bar[y1]]] - 5*Tr[y1 . Trans[Bar[y1]] . y1 . Trans[Bar[y1]]]
```

**beta_lambdaS2**

```text
(9*g2^4)/4 + (3*gY^4)/4 + 4*lambda12^2 + lambda12Adj^2/4 + 4*lambdaH2^2 + lambdaH2Adj^2/4 + gY^2*((3*g2^2)/2 - 3*lambdaS2) - 9*g2^2*lambdaS2 + 12*lambdaS2^2 + 8*lambdaHHdagS2S2*Bar[lambdaHHdagS2S2] + 8*lambdaS1bar2S2bar2*Bar[lambdaS1bar2S2bar2] + 24*lambdaS1barS2S2bar2*Bar[lambdaS1barS2S2bar2] + 6*lambdaS2*Tr[y2 . Trans[Bar[y2]]] - 5*Tr[y2 . Trans[Bar[y2]] . y2 . Trans[Bar[y2]]]
```

**beta_lambdaH1**

```text
(9*g2^4)/4 + (3*gY^4)/4 - 9*g2^2*lambdaH1 - 3*gY^2*lambdaH1 + 6*lambdaH*lambdaH1 + 4*lambdaH1^2 + (3*lambdaH1Adj^2)/4 + 4*lambda12*lambdaH2 + 6*lambdaH1*lambdaS1 + 24*lambdaHHdagS1barS1bar*Bar[lambdaHHdagS1barS1bar] + 6*lambdaT3*Bar[lambdaT3] + 3*lambdaS1S1bar2S2bar*(Bar[lambdaHHdagS1barS2barCross] - Bar[lambdaT3GeneralQ110101I1]) + lambdaHHdagS1barS2barCross*(4*Bar[lambdaHHdagS1barS2barCross] + 3*Bar[lambdaS1S1bar2S2bar] + 2*Bar[lambdaT3GeneralQ110101I1]) + lambdaT3GeneralQ110101I1*(2*Bar[lambdaHHdagS1barS2barCross] - 3*Bar[lambdaS1S1bar2S2bar] + 4*Bar[lambdaT3GeneralQ110101I1]) + lambdaH1*(3*Tr[y1 . Trans[Bar[y1]]] + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]]) - 3*Tr[y1 . Trans[Bar[y1]] . ye . Trans[Bar[ye]]]
```

**beta_lambdaH2**

```text
(9*g2^4)/4 + (3*gY^4)/4 + 4*lambda12*lambdaH1 - 9*g2^2*lambdaH2 - 3*gY^2*lambdaH2 + 6*lambdaH*lambdaH2 + 4*lambdaH2^2 + (3*lambdaH2Adj^2)/4 + 6*lambdaH2*lambdaS2 + 24*lambdaHHdagS2S2*Bar[lambdaHHdagS2S2] + 6*lambdaT3*Bar[lambdaT3] + 3*lambdaS1barS2S2bar2*(Bar[lambdaHHdagS1barS2barCross] - Bar[lambdaT3GeneralQ110101I1]) + lambdaHHdagS1barS2barCross*(4*Bar[lambdaHHdagS1barS2barCross] + 3*Bar[lambdaS1barS2S2bar2] + 2*Bar[lambdaT3GeneralQ110101I1]) + lambdaT3GeneralQ110101I1*(2*Bar[lambdaHHdagS1barS2barCross] - 3*Bar[lambdaS1barS2S2bar2] + 4*Bar[lambdaT3GeneralQ110101I1]) + lambdaH2*(3*Tr[y2 . Trans[Bar[y2]]] + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]]) - 3*Tr[y2 . Trans[Bar[y2]] . ye . Trans[Bar[ye]]]
```

**beta_lambda12**

```text
(9*g2^4)/4 + (3*gY^4)/4 - 9*g2^2*lambda12 - 3*gY^2*lambda12 + 4*lambda12^2 + (3*lambda12Adj^2)/4 + 4*lambdaH1*lambdaH2 + 6*lambda12*lambdaS1 + 6*lambda12*lambdaS2 + 24*lambdaS1bar2S2bar2*Bar[lambdaS1bar2S2bar2] + 9*lambdaS1barS2S2bar2*(Bar[lambdaS1barS2S2bar2] + Bar[lambdaS1S1bar2S2bar]) + 9*lambdaS1S1bar2S2bar*(Bar[lambdaS1barS2S2bar2] + Bar[lambdaS1S1bar2S2bar]) + 3*lambdaT3*Bar[lambdaT3] + lambdaHHdagS1barS2barCross*(2*Bar[lambdaHHdagS1barS2barCross] + Bar[lambdaT3GeneralQ110101I1]) + lambdaT3GeneralQ110101I1*(Bar[lambdaHHdagS1barS2barCross] + 2*Bar[lambdaT3GeneralQ110101I1]) + lambda12*(3*Tr[y1 . Trans[Bar[y1]]] + 3*Tr[y2 . Trans[Bar[y2]]]) - (9*Tr[y1 . Trans[Bar[y1]] . y2 . Trans[Bar[y2]]])/2 - 3*Tr[y1 . Trans[Bar[y2]] . y2 . Trans[Bar[y1]]]
```

**beta_lambdaT3**

```text
-9*g2^2*lambdaT3 - 3*gY^2*lambdaT3 + 2*lambda12*lambdaT3 + 2*lambdaH*lambdaT3 + 4*lambdaH1*lambdaT3 + 4*lambdaH2*lambdaT3 + 3*Bar[lambdaHHdagS2S2]*Tr[y1 . Trans[Bar[y2]]] - Bar[lambdaHHdagS1barS1bar]*Tr[y2 . Trans[Bar[y1]]] + lambdaT3*(lambda12Adj/6 + (2*lambdaH1Adj)/3 + 2*lambdaH2Adj + (3*Tr[y1 . Trans[Bar[y1]]])/2 + (3*Tr[y2 . Trans[Bar[y2]]])/2 + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaH1Adj**

```text
gY^2*(-6*g2^2 - 3*lambdaH1Adj) - 9*g2^2*lambdaH1Adj + 2*lambdaH*lambdaH1Adj + 8*lambdaH1*lambdaH1Adj + lambda12Adj*lambdaH2Adj + 2*lambdaH1Adj*lambdaS1 - 64*lambdaHHdagS1barS1bar*Bar[lambdaHHdagS1barS1bar] + 4*lambdaHHdagS1barS2barCross*(2*Bar[lambdaHHdagS1barS2barCross] + Bar[lambdaS1S1bar2S2bar]) + (16*lambdaT3*Bar[lambdaT3])/3 + 4*lambdaT3GeneralQ110101I1*(Bar[lambdaS1S1bar2S2bar] - 2*Bar[lambdaT3GeneralQ110101I1]) + 4*lambdaS1S1bar2S2bar*(Bar[lambdaHHdagS1barS2barCross] + Bar[lambdaT3GeneralQ110101I1]) + lambdaH1Adj*(3*Tr[y1 . Trans[Bar[y1]]] + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]]) - (4*Tr[y1 . Trans[Bar[y1]] . ye . Trans[Bar[ye]]])/3
```

**beta_lambdaH2Adj**

```text
lambda12Adj*lambdaH1Adj + gY^2*(6*g2^2 - 3*lambdaH2Adj) - 9*g2^2*lambdaH2Adj + 2*lambdaH*lambdaH2Adj + 8*lambdaH2*lambdaH2Adj + 2*lambdaH2Adj*lambdaS2 + 64*lambdaHHdagS2S2*Bar[lambdaHHdagS2S2] - 4*lambdaHHdagS1barS2barCross*(2*Bar[lambdaHHdagS1barS2barCross] + Bar[lambdaS1barS2S2bar2]) + 16*lambdaT3*Bar[lambdaT3] - 4*lambdaS1barS2S2bar2*(Bar[lambdaHHdagS1barS2barCross] + Bar[lambdaT3GeneralQ110101I1]) + lambdaT3GeneralQ110101I1*(-4*Bar[lambdaS1barS2S2bar2] + 8*Bar[lambdaT3GeneralQ110101I1]) + lambdaH2Adj*(3*Tr[y2 . Trans[Bar[y2]]] + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]]) - 4*Tr[y2 . Trans[Bar[y2]] . ye . Trans[Bar[ye]]]
```

**beta_lambda12Adj**

```text
gY^2*(-6*g2^2 - 3*lambda12Adj) - 9*g2^2*lambda12Adj + 8*lambda12*lambda12Adj + lambdaH1Adj*lambdaH2Adj + 2*lambda12Adj*lambdaS1 + 2*lambda12Adj*lambdaS2 + 4*lambdaT3GeneralQ110101I1*Bar[lambdaHHdagS1barS2barCross] - 64*lambdaS1bar2S2bar2*Bar[lambdaS1bar2S2bar2] - 4*lambdaS1barS2S2bar2*(5*Bar[lambdaS1barS2S2bar2] + Bar[lambdaS1S1bar2S2bar]) - 4*lambdaS1S1bar2S2bar*(Bar[lambdaS1barS2S2bar2] + 5*Bar[lambdaS1S1bar2S2bar]) + (4*lambdaT3*Bar[lambdaT3])/3 + 4*lambdaHHdagS1barS2barCross*Bar[lambdaT3GeneralQ110101I1] + lambda12Adj*(3*Tr[y1 . Trans[Bar[y1]]] + 3*Tr[y2 . Trans[Bar[y2]]]) - (2*Tr[y1 . Trans[Bar[y1]] . y2 . Trans[Bar[y2]]])/3 - (8*Tr[y1 . Trans[Bar[y2]] . y2 . Trans[Bar[y1]]])/3
```

**beta_lambdaHHdagS2S2**

```text
-9*g2^2*lambdaHHdagS2S2 - 3*gY^2*lambdaHHdagS2S2 + 2*lambdaH*lambdaHHdagS2S2 + 8*lambdaH2*lambdaHHdagS2S2 + 4*lambdaH2Adj*lambdaHHdagS2S2 + 2*lambdaHHdagS2S2*lambdaS2 + 4*lambdaHHdagS1barS1bar*Bar[lambdaS1bar2S2bar2] + (3*Bar[lambdaT3]*Tr[y1 . Trans[Bar[y2]]])/2 + lambdaHHdagS2S2*(3*Tr[y2 . Trans[Bar[y2]]] + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaHHdagS1barS1bar**

```text
-9*g2^2*lambdaHHdagS1barS1bar - 3*gY^2*lambdaHHdagS1barS1bar + 2*lambdaH*lambdaHHdagS1barS1bar + 8*lambdaH1*lambdaHHdagS1barS1bar - 4*lambdaH1Adj*lambdaHHdagS1barS1bar + 2*lambdaHHdagS1barS1bar*lambdaS1 + 4*lambdaHHdagS2S2*lambdaS1bar2S2bar2 - (Bar[lambdaT3]*Tr[y2 . Trans[Bar[y1]]])/2 + lambdaHHdagS1barS1bar*(3*Tr[y1 . Trans[Bar[y1]]] + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaS1bar2S2bar2**

```text
lambdaHHdagS1barS2barCross^2 - 9*g2^2*lambdaS1bar2S2bar2 - 3*gY^2*lambdaS1bar2S2bar2 + 8*lambda12*lambdaS1bar2S2bar2 - 4*lambda12Adj*lambdaS1bar2S2bar2 + 2*lambdaS1*lambdaS1bar2S2bar2 + 5*lambdaS1barS2S2bar2^2 + 2*lambdaS1barS2S2bar2*lambdaS1S1bar2S2bar + 5*lambdaS1S1bar2S2bar^2 + 2*lambdaS1bar2S2bar2*lambdaS2 + lambdaT3GeneralQ110101I1^2 + 4*lambdaHHdagS1barS1bar*Bar[lambdaHHdagS2S2] + lambdaS1bar2S2bar2*(3*Tr[y1 . Trans[Bar[y1]]] + 3*Tr[y2 . Trans[Bar[y2]]]) + (5*Tr[y2 . Trans[Bar[y1]] . y2 . Trans[Bar[y1]]])/6
```

**beta_lambdaS1barS2S2bar2**

```text
-9*g2^2*lambdaS1barS2S2bar2 - 3*gY^2*lambdaS1barS2S2bar2 + lambda12Adj*((-5*lambdaS1barS2S2bar2)/2 - lambdaS1S1bar2S2bar/2) + lambda12*(6*lambdaS1barS2S2bar2 + 6*lambdaS1S1bar2S2bar) + 12*lambdaS1barS2S2bar2*lambdaS2 + lambdaH2*(2*lambdaHHdagS1barS2barCross - 2*lambdaT3GeneralQ110101I1) + lambdaH2Adj*(-1/2*lambdaHHdagS1barS2barCross - lambdaT3GeneralQ110101I1/2) + 4*lambdaS1bar2S2bar2*(5*Bar[lambdaS1barS2S2bar2] + Bar[lambdaS1S1bar2S2bar]) + lambdaS1barS2S2bar2*((3*Tr[y1 . Trans[Bar[y1]]])/2 + (9*Tr[y2 . Trans[Bar[y2]]])/2)
```

**beta_lambdaS1S1bar2S2bar**

```text
lambda12Adj*(-1/2*lambdaS1barS2S2bar2 - (5*lambdaS1S1bar2S2bar)/2) - 9*g2^2*lambdaS1S1bar2S2bar - 3*gY^2*lambdaS1S1bar2S2bar + 12*lambdaS1*lambdaS1S1bar2S2bar + lambda12*(6*lambdaS1barS2S2bar2 + 6*lambdaS1S1bar2S2bar) + lambdaH1*(2*lambdaHHdagS1barS2barCross - 2*lambdaT3GeneralQ110101I1) + lambdaH1Adj*(lambdaHHdagS1barS2barCross/2 + lambdaT3GeneralQ110101I1/2) + 4*lambdaS1bar2S2bar2*(Bar[lambdaS1barS2S2bar2] + 5*Bar[lambdaS1S1bar2S2bar]) + lambdaS1S1bar2S2bar*((9*Tr[y1 . Trans[Bar[y1]]])/2 + (3*Tr[y2 . Trans[Bar[y2]]])/2)
```

**beta_lambdaHHdagS1barS2barCross**

```text
-9*g2^2*lambdaHHdagS1barS2barCross - 3*gY^2*lambdaHHdagS1barS2barCross + 2*lambda12*lambdaHHdagS1barS2barCross + lambdaH2*(4*lambdaHHdagS1barS2barCross + 6*lambdaS1barS2S2bar2) + lambdaH1*(4*lambdaHHdagS1barS2barCross + 6*lambdaS1S1bar2S2bar) + lambdaH*(4*lambdaHHdagS1barS2barCross - 2*lambdaT3GeneralQ110101I1) + lambdaH2Adj*(-2*lambdaHHdagS1barS2barCross - lambdaS1barS2S2bar2/2 - lambdaT3GeneralQ110101I1) + lambda12Adj*(-1/2*lambdaHHdagS1barS2barCross + lambdaT3GeneralQ110101I1) + lambdaH1Adj*(2*lambdaHHdagS1barS2barCross + lambdaS1S1bar2S2bar/2 + lambdaT3GeneralQ110101I1) + lambdaS1bar2S2bar2*(8*Bar[lambdaHHdagS1barS2barCross] - 4*Bar[lambdaT3GeneralQ110101I1]) + lambdaHHdagS1barS2barCross*((3*Tr[y1 . Trans[Bar[y1]]])/2 + (3*Tr[y2 . Trans[Bar[y2]]])/2 + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaT3GeneralQ110101I1**

```text
lambdaH1Adj*(-lambdaHHdagS1barS2barCross + lambdaS1S1bar2S2bar/2 - 2*lambdaT3GeneralQ110101I1) + lambda12Adj*(lambdaHHdagS1barS2barCross - lambdaT3GeneralQ110101I1/2) - 9*g2^2*lambdaT3GeneralQ110101I1 - 3*gY^2*lambdaT3GeneralQ110101I1 + 2*lambda12*lambdaT3GeneralQ110101I1 + lambdaH2Adj*(lambdaHHdagS1barS2barCross - lambdaS1barS2S2bar2/2 + 2*lambdaT3GeneralQ110101I1) + lambdaH*(-2*lambdaHHdagS1barS2barCross + 4*lambdaT3GeneralQ110101I1) + lambdaH2*(-6*lambdaS1barS2S2bar2 + 4*lambdaT3GeneralQ110101I1) + lambdaH1*(-6*lambdaS1S1bar2S2bar + 4*lambdaT3GeneralQ110101I1) + lambdaS1bar2S2bar2*(-4*Bar[lambdaHHdagS1barS2barCross] + 8*Bar[lambdaT3GeneralQ110101I1]) + lambdaT3GeneralQ110101I1*((3*Tr[y1 . Trans[Bar[y1]]])/2 + (3*Tr[y2 . Trans[Bar[y2]]])/2 + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_muT3GeneralS1S2**

```text
6*lambdaS1S1bar2S2bar*mS1Sq + 6*lambdaS1barS2S2bar2*mS2Sq - (9*g2^2*muT3GeneralS1S2)/2 - (3*gY^2*muT3GeneralS1S2)/2 + 2*lambda12*muT3GeneralS1S2 - (3*lambda12Adj*muT3GeneralS1S2)/2 + 12*lambdaS1bar2S2bar2*Bar[muT3GeneralS1S2] + muT3GeneralS1S2*((3*Tr[y1 . Trans[Bar[y1]]])/2 + (3*Tr[y2 . Trans[Bar[y2]]])/2)
```

### EFT1: 26 beta functions

**beta_gY**

```text
(43*gY^3)/6
```

**beta_g2**

```text
(-17*g2^3)/6
```

**beta_g3**

```text
-7*g3^3
```

**beta_yu**

```text
(-9*g2^2*Matrix[yu][gen[$i], gen[$j]])/4 - 8*g3^2*Matrix[yu][gen[$i], gen[$j]] - (17*gY^2*Matrix[yu][gen[$i], gen[$j]])/12 + 3*Tr[yd . Trans[Bar[yd]]]*Matrix[yu][gen[$i], gen[$j]] + Tr[ye . Trans[Bar[ye]]]*Matrix[yu][gen[$i], gen[$j]] + 3*Tr[yu . Trans[Bar[yu]]]*Matrix[yu][gen[$i], gen[$j]] - (3*(Matrix[yd, Trans[Bar[yd]], yu][gen[$i], gen[$j]] - Matrix[yu, Trans[Bar[yu]], yu][gen[$i], gen[$j]]))/2
```

**beta_yd**

```text
(-9*g2^2*Matrix[yd][gen[$i], gen[$j]])/4 - 8*g3^2*Matrix[yd][gen[$i], gen[$j]] - (5*gY^2*Matrix[yd][gen[$i], gen[$j]])/12 + 3*Tr[yd . Trans[Bar[yd]]]*Matrix[yd][gen[$i], gen[$j]] + Tr[ye . Trans[Bar[ye]]]*Matrix[yd][gen[$i], gen[$j]] + 3*Tr[yu . Trans[Bar[yu]]]*Matrix[yd][gen[$i], gen[$j]] + (3*(Matrix[yd, Trans[Bar[yd]], yd][gen[$i], gen[$j]] - Matrix[yu, Trans[Bar[yu]], yd][gen[$i], gen[$j]]))/2
```

**beta_ye**

```text
(-9*g2^2*Matrix[ye][gen[$i], gen[$j]])/4 - (15*gY^2*Matrix[ye][gen[$i], gen[$j]])/4 + 3*Tr[yd . Trans[Bar[yd]]]*Matrix[ye][gen[$i], gen[$j]] + Tr[ye . Trans[Bar[ye]]]*Matrix[ye][gen[$i], gen[$j]] + 3*Tr[yu . Trans[Bar[yu]]]*Matrix[ye][gen[$i], gen[$j]] + (3*Matrix[ye, Trans[Bar[ye]], ye][gen[$i], gen[$j]])/2
```

**beta_mS1Sq**

```text
(-9*g2^2*mS1Sq)/2 - (3*gY^2*mS1Sq)/2 + 6*lambdaS1*mS1Sq + 4*lambda12*mS2Sq + 6*muT3GeneralS1S2*Bar[lambdaS1S1bar2S2bar] + 6*lambdaS1S1bar2S2bar*Bar[muT3GeneralS1S2]
```

**beta_mS2Sq**

```text
4*lambda12*mS1Sq - (9*g2^2*mS2Sq)/2 - (3*gY^2*mS2Sq)/2 + 6*lambdaS2*mS2Sq + 6*muT3GeneralS1S2*Bar[lambdaS1barS2S2bar2] + 6*lambdaS1barS2S2bar2*Bar[muT3GeneralS1S2]
```

**beta_lambdaH**

```text
(9*g2^4)/4 + (3*gY^4)/4 + gY^2*((3*g2^2)/2 - 3*lambdaH) - 9*g2^2*lambdaH + 12*lambdaH^2 + 4*lambdaH1^2 + lambdaH1Adj^2/4 + 4*lambdaH2^2 + lambdaH2Adj^2/4 + 8*lambdaHHdagS1barS1bar*Bar[lambdaHHdagS1barS1bar] + 4*lambdaHHdagS1barS2barCross*Bar[lambdaHHdagS1barS2barCross] + 8*lambdaHHdagS2S2*Bar[lambdaHHdagS2S2] + 4*lambdaT3*Bar[lambdaT3] + 4*lambdaT3GeneralQ110101I1*Bar[lambdaT3GeneralQ110101I1] + lambdaH*(12*Tr[yd . Trans[Bar[yd]]] + 4*Tr[ye . Trans[Bar[ye]]] + 12*Tr[yu . Trans[Bar[yu]]]) - 12*Tr[yd . Trans[Bar[yd]] . yd . Trans[Bar[yd]]] - 4*Tr[ye . Trans[Bar[ye]] . ye . Trans[Bar[ye]]] - 12*Tr[yu . Trans[Bar[yu]] . yu . Trans[Bar[yu]]]
```

**beta_lambdaS1**

```text
(9*g2^4)/4 + (3*gY^4)/4 + 4*lambda12^2 + lambda12Adj^2/4 + 4*lambdaH1^2 + lambdaH1Adj^2/4 + gY^2*((3*g2^2)/2 - 3*lambdaS1) - 9*g2^2*lambdaS1 + 12*lambdaS1^2 + 8*lambdaHHdagS1barS1bar*Bar[lambdaHHdagS1barS1bar] + 8*lambdaS1bar2S2bar2*Bar[lambdaS1bar2S2bar2] + 24*lambdaS1S1bar2S2bar*Bar[lambdaS1S1bar2S2bar]
```

**beta_lambdaS2**

```text
(9*g2^4)/4 + (3*gY^4)/4 + 4*lambda12^2 + lambda12Adj^2/4 + 4*lambdaH2^2 + lambdaH2Adj^2/4 + gY^2*((3*g2^2)/2 - 3*lambdaS2) - 9*g2^2*lambdaS2 + 12*lambdaS2^2 + 8*lambdaHHdagS2S2*Bar[lambdaHHdagS2S2] + 8*lambdaS1bar2S2bar2*Bar[lambdaS1bar2S2bar2] + 24*lambdaS1barS2S2bar2*Bar[lambdaS1barS2S2bar2]
```

**beta_lambdaH1**

```text
(9*g2^4)/4 + (3*gY^4)/4 - 9*g2^2*lambdaH1 - 3*gY^2*lambdaH1 + 6*lambdaH*lambdaH1 + 4*lambdaH1^2 + (3*lambdaH1Adj^2)/4 + 4*lambda12*lambdaH2 + 6*lambdaH1*lambdaS1 + 24*lambdaHHdagS1barS1bar*Bar[lambdaHHdagS1barS1bar] + 6*lambdaT3*Bar[lambdaT3] + 3*lambdaS1S1bar2S2bar*(Bar[lambdaHHdagS1barS2barCross] - Bar[lambdaT3GeneralQ110101I1]) + lambdaHHdagS1barS2barCross*(4*Bar[lambdaHHdagS1barS2barCross] + 3*Bar[lambdaS1S1bar2S2bar] + 2*Bar[lambdaT3GeneralQ110101I1]) + lambdaT3GeneralQ110101I1*(2*Bar[lambdaHHdagS1barS2barCross] - 3*Bar[lambdaS1S1bar2S2bar] + 4*Bar[lambdaT3GeneralQ110101I1]) + lambdaH1*(6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaH2**

```text
(9*g2^4)/4 + (3*gY^4)/4 + 4*lambda12*lambdaH1 - 9*g2^2*lambdaH2 - 3*gY^2*lambdaH2 + 6*lambdaH*lambdaH2 + 4*lambdaH2^2 + (3*lambdaH2Adj^2)/4 + 6*lambdaH2*lambdaS2 + 24*lambdaHHdagS2S2*Bar[lambdaHHdagS2S2] + 6*lambdaT3*Bar[lambdaT3] + 3*lambdaS1barS2S2bar2*(Bar[lambdaHHdagS1barS2barCross] - Bar[lambdaT3GeneralQ110101I1]) + lambdaHHdagS1barS2barCross*(4*Bar[lambdaHHdagS1barS2barCross] + 3*Bar[lambdaS1barS2S2bar2] + 2*Bar[lambdaT3GeneralQ110101I1]) + lambdaT3GeneralQ110101I1*(2*Bar[lambdaHHdagS1barS2barCross] - 3*Bar[lambdaS1barS2S2bar2] + 4*Bar[lambdaT3GeneralQ110101I1]) + lambdaH2*(6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambda12**

```text
(9*g2^4)/4 + (3*gY^4)/4 - 9*g2^2*lambda12 - 3*gY^2*lambda12 + 4*lambda12^2 + (3*lambda12Adj^2)/4 + 4*lambdaH1*lambdaH2 + 6*lambda12*lambdaS1 + 6*lambda12*lambdaS2 + 24*lambdaS1bar2S2bar2*Bar[lambdaS1bar2S2bar2] + 9*lambdaS1barS2S2bar2*(Bar[lambdaS1barS2S2bar2] + Bar[lambdaS1S1bar2S2bar]) + 9*lambdaS1S1bar2S2bar*(Bar[lambdaS1barS2S2bar2] + Bar[lambdaS1S1bar2S2bar]) + 3*lambdaT3*Bar[lambdaT3] + lambdaHHdagS1barS2barCross*(2*Bar[lambdaHHdagS1barS2barCross] + Bar[lambdaT3GeneralQ110101I1]) + lambdaT3GeneralQ110101I1*(Bar[lambdaHHdagS1barS2barCross] + 2*Bar[lambdaT3GeneralQ110101I1])
```

**beta_lambdaT3**

```text
-9*g2^2*lambdaT3 - 3*gY^2*lambdaT3 + 2*lambda12*lambdaT3 + 2*lambdaH*lambdaT3 + 4*lambdaH1*lambdaT3 + 4*lambdaH2*lambdaT3 + lambdaT3*(lambda12Adj/6 + (2*lambdaH1Adj)/3 + 2*lambdaH2Adj + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaH1Adj**

```text
gY^2*(-6*g2^2 - 3*lambdaH1Adj) - 9*g2^2*lambdaH1Adj + 2*lambdaH*lambdaH1Adj + 8*lambdaH1*lambdaH1Adj + lambda12Adj*lambdaH2Adj + 2*lambdaH1Adj*lambdaS1 - 64*lambdaHHdagS1barS1bar*Bar[lambdaHHdagS1barS1bar] + 4*lambdaHHdagS1barS2barCross*(2*Bar[lambdaHHdagS1barS2barCross] + Bar[lambdaS1S1bar2S2bar]) + (16*lambdaT3*Bar[lambdaT3])/3 + 4*lambdaT3GeneralQ110101I1*(Bar[lambdaS1S1bar2S2bar] - 2*Bar[lambdaT3GeneralQ110101I1]) + 4*lambdaS1S1bar2S2bar*(Bar[lambdaHHdagS1barS2barCross] + Bar[lambdaT3GeneralQ110101I1]) + lambdaH1Adj*(6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaH2Adj**

```text
lambda12Adj*lambdaH1Adj + gY^2*(6*g2^2 - 3*lambdaH2Adj) - 9*g2^2*lambdaH2Adj + 2*lambdaH*lambdaH2Adj + 8*lambdaH2*lambdaH2Adj + 2*lambdaH2Adj*lambdaS2 + 64*lambdaHHdagS2S2*Bar[lambdaHHdagS2S2] - 4*lambdaHHdagS1barS2barCross*(2*Bar[lambdaHHdagS1barS2barCross] + Bar[lambdaS1barS2S2bar2]) + 16*lambdaT3*Bar[lambdaT3] - 4*lambdaS1barS2S2bar2*(Bar[lambdaHHdagS1barS2barCross] + Bar[lambdaT3GeneralQ110101I1]) + lambdaT3GeneralQ110101I1*(-4*Bar[lambdaS1barS2S2bar2] + 8*Bar[lambdaT3GeneralQ110101I1]) + lambdaH2Adj*(6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambda12Adj**

```text
gY^2*(-6*g2^2 - 3*lambda12Adj) - 9*g2^2*lambda12Adj + 8*lambda12*lambda12Adj + lambdaH1Adj*lambdaH2Adj + 2*lambda12Adj*lambdaS1 + 2*lambda12Adj*lambdaS2 + 4*lambdaT3GeneralQ110101I1*Bar[lambdaHHdagS1barS2barCross] - 64*lambdaS1bar2S2bar2*Bar[lambdaS1bar2S2bar2] - 4*lambdaS1barS2S2bar2*(5*Bar[lambdaS1barS2S2bar2] + Bar[lambdaS1S1bar2S2bar]) - 4*lambdaS1S1bar2S2bar*(Bar[lambdaS1barS2S2bar2] + 5*Bar[lambdaS1S1bar2S2bar]) + (4*lambdaT3*Bar[lambdaT3])/3 + 4*lambdaHHdagS1barS2barCross*Bar[lambdaT3GeneralQ110101I1]
```

**beta_lambdaHHdagS2S2**

```text
-9*g2^2*lambdaHHdagS2S2 - 3*gY^2*lambdaHHdagS2S2 + 2*lambdaH*lambdaHHdagS2S2 + 8*lambdaH2*lambdaHHdagS2S2 + 4*lambdaH2Adj*lambdaHHdagS2S2 + 2*lambdaHHdagS2S2*lambdaS2 + 4*lambdaHHdagS1barS1bar*Bar[lambdaS1bar2S2bar2] + lambdaHHdagS2S2*(6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaHHdagS1barS1bar**

```text
-9*g2^2*lambdaHHdagS1barS1bar - 3*gY^2*lambdaHHdagS1barS1bar + 2*lambdaH*lambdaHHdagS1barS1bar + 8*lambdaH1*lambdaHHdagS1barS1bar - 4*lambdaH1Adj*lambdaHHdagS1barS1bar + 2*lambdaHHdagS1barS1bar*lambdaS1 + 4*lambdaHHdagS2S2*lambdaS1bar2S2bar2 + lambdaHHdagS1barS1bar*(6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaS1bar2S2bar2**

```text
lambdaHHdagS1barS2barCross^2 - 9*g2^2*lambdaS1bar2S2bar2 - 3*gY^2*lambdaS1bar2S2bar2 + 8*lambda12*lambdaS1bar2S2bar2 - 4*lambda12Adj*lambdaS1bar2S2bar2 + 2*lambdaS1*lambdaS1bar2S2bar2 + 5*lambdaS1barS2S2bar2^2 + 2*lambdaS1barS2S2bar2*lambdaS1S1bar2S2bar + 5*lambdaS1S1bar2S2bar^2 + 2*lambdaS1bar2S2bar2*lambdaS2 + lambdaT3GeneralQ110101I1^2 + 4*lambdaHHdagS1barS1bar*Bar[lambdaHHdagS2S2]
```

**beta_lambdaS1barS2S2bar2**

```text
-9*g2^2*lambdaS1barS2S2bar2 - 3*gY^2*lambdaS1barS2S2bar2 + lambda12Adj*((-5*lambdaS1barS2S2bar2)/2 - lambdaS1S1bar2S2bar/2) + lambda12*(6*lambdaS1barS2S2bar2 + 6*lambdaS1S1bar2S2bar) + 12*lambdaS1barS2S2bar2*lambdaS2 + lambdaH2*(2*lambdaHHdagS1barS2barCross - 2*lambdaT3GeneralQ110101I1) + lambdaH2Adj*(-1/2*lambdaHHdagS1barS2barCross - lambdaT3GeneralQ110101I1/2) + 4*lambdaS1bar2S2bar2*(5*Bar[lambdaS1barS2S2bar2] + Bar[lambdaS1S1bar2S2bar])
```

**beta_lambdaS1S1bar2S2bar**

```text
lambda12Adj*(-1/2*lambdaS1barS2S2bar2 - (5*lambdaS1S1bar2S2bar)/2) - 9*g2^2*lambdaS1S1bar2S2bar - 3*gY^2*lambdaS1S1bar2S2bar + 12*lambdaS1*lambdaS1S1bar2S2bar + lambda12*(6*lambdaS1barS2S2bar2 + 6*lambdaS1S1bar2S2bar) + lambdaH1*(2*lambdaHHdagS1barS2barCross - 2*lambdaT3GeneralQ110101I1) + lambdaH1Adj*(lambdaHHdagS1barS2barCross/2 + lambdaT3GeneralQ110101I1/2) + 4*lambdaS1bar2S2bar2*(Bar[lambdaS1barS2S2bar2] + 5*Bar[lambdaS1S1bar2S2bar])
```

**beta_lambdaHHdagS1barS2barCross**

```text
-9*g2^2*lambdaHHdagS1barS2barCross - 3*gY^2*lambdaHHdagS1barS2barCross + 2*lambda12*lambdaHHdagS1barS2barCross + lambdaH2*(4*lambdaHHdagS1barS2barCross + 6*lambdaS1barS2S2bar2) + lambdaH1*(4*lambdaHHdagS1barS2barCross + 6*lambdaS1S1bar2S2bar) + lambdaH*(4*lambdaHHdagS1barS2barCross - 2*lambdaT3GeneralQ110101I1) + lambdaH2Adj*(-2*lambdaHHdagS1barS2barCross - lambdaS1barS2S2bar2/2 - lambdaT3GeneralQ110101I1) + lambda12Adj*(-1/2*lambdaHHdagS1barS2barCross + lambdaT3GeneralQ110101I1) + lambdaH1Adj*(2*lambdaHHdagS1barS2barCross + lambdaS1S1bar2S2bar/2 + lambdaT3GeneralQ110101I1) + lambdaS1bar2S2bar2*(8*Bar[lambdaHHdagS1barS2barCross] - 4*Bar[lambdaT3GeneralQ110101I1]) + lambdaHHdagS1barS2barCross*(6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaT3GeneralQ110101I1**

```text
lambdaH1Adj*(-lambdaHHdagS1barS2barCross + lambdaS1S1bar2S2bar/2 - 2*lambdaT3GeneralQ110101I1) + lambda12Adj*(lambdaHHdagS1barS2barCross - lambdaT3GeneralQ110101I1/2) - 9*g2^2*lambdaT3GeneralQ110101I1 - 3*gY^2*lambdaT3GeneralQ110101I1 + 2*lambda12*lambdaT3GeneralQ110101I1 + lambdaH2Adj*(lambdaHHdagS1barS2barCross - lambdaS1barS2S2bar2/2 + 2*lambdaT3GeneralQ110101I1) + lambdaH*(-2*lambdaHHdagS1barS2barCross + 4*lambdaT3GeneralQ110101I1) + lambdaH2*(-6*lambdaS1barS2S2bar2 + 4*lambdaT3GeneralQ110101I1) + lambdaH1*(-6*lambdaS1S1bar2S2bar + 4*lambdaT3GeneralQ110101I1) + lambdaS1bar2S2bar2*(-4*Bar[lambdaHHdagS1barS2barCross] + 8*Bar[lambdaT3GeneralQ110101I1]) + lambdaT3GeneralQ110101I1*(6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_muT3GeneralS1S2**

```text
6*lambdaS1S1bar2S2bar*mS1Sq + 6*lambdaS1barS2S2bar2*mS2Sq - (9*g2^2*muT3GeneralS1S2)/2 - (3*gY^2*muT3GeneralS1S2)/2 + 2*lambda12*muT3GeneralS1S2 - (3*lambda12Adj*muT3GeneralS1S2)/2 + 12*lambdaS1bar2S2bar2*Bar[muT3GeneralS1S2]
```

## Interpretation and limitations

- UV RGEs include symmetry-allowed scalar quartics, bilinears, and the additional Yukawa where gauge invariance permits it.
- EFT1 RGEs describe the renormalisable SM+S1+S2 theory; the heavy fermion and its Yukawas have been removed.
- EFT1 higher-dimensional Wilson RGEs and final C5 have **not** been rederived with the full general interactions.
- Mixed scalar mass terms require mass-basis diagonalisation for B(-1) and C(-1) before general matching.
- The existing heavy-fermion flavour-mass parameterisation is not certified as the most general non-diagonal flavour mass basis.

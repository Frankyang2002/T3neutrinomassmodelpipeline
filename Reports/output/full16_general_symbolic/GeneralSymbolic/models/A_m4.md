# T3 A_-4: general-interaction symbolic results

Fields: S1=(1,-2), S2=(3,-1), F=(2,-3/2) [SU(2) dimension, hypercharge].

Additional independent quartics: **0**; additional complex bilinears: **0**; additional Yukawa: **none**.

Heavy scalar mass basis: **CheckedNotMatched**. General sequential matching: **not performed**. Matchete/RGBeta physical coupling dictionary: **not certified**.

## Independent interaction inventory

- `Yukawa1_GenericCG`
- `Yukawa2_GenericCG`
- `Scalar1SelfInv1`
- `Scalar2SelfInv1`
- `Scalar2SelfInv2`
- `HiggsPortal1Inv1`
- `HiggsPortal2Inv1`
- `HiggsPortal2Inv2`
- `ScalarCrossPortalInv1`
- `T3ScalarMix_GenericCG`

## Complete BSM UV Lagrangian (native Matchete InputForm)

The expression below includes BSM free fields and interactions; it excludes the pure Standard Model Lagrangian.

```text
-(Bar[Field[NewScalar1, Scalar, {}, {}]]*Coupling[MS1, {}, 0]^2*Field[NewScalar1, Scalar, {}, {}]) + Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar1, Scalar, {}, {}]]*CG[T3HiggsPortal1CGInv1, {Index[d$$1, SU2L[fund]], Bar[Index[d$$2, SU2L[fund]]]}]*Coupling[lambdaH1Inv1, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar1, Scalar, {}, {}] + Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3MixCG, {Bar[Index[d$$1, SU2L[fund]]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$1, T3BSMd3]}]*Coupling[lambdaT3, {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar1, Scalar, {}, {}] + (Bar[Field[NewScalar1, Scalar, {}, {}]]^2*Coupling[lambdaS1Inv1, {}, 0]*Field[NewScalar1, Scalar, {}, {}]^2)/2 + Bar[Field[NewScalar1, Scalar, {}, {Index[d$$1, Lorentz]}]]*Field[NewScalar1, Scalar, {}, {Index[d$$1, Lorentz]}] + Bar[Coupling[lambdaT3, {}, 0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*Bar[Field[NewScalar1, Scalar, {}, {}]]*CG[Bar[T3MixCG], {Index[d$$1, SU2L[fund]], Index[d$$2, SU2L[fund]], Index[d$$1, T3BSMd3]}]*Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}] - Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*Coupling[MS2, {}, 0]^2*Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}] + Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {Index[d$$1, Lorentz]}]]*Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {Index[d$$1, Lorentz]}] + (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3HiggsPortal2CGInv1, {Index[d$$1, SU2L[fund]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3]}]*Coupling[lambdaH2Inv1, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}])/2 + (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3HiggsPortal2CGInv1, {Index[d$$1, SU2L[fund]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$2, T3BSMd3], Index[d$$1, T3BSMd3]}]*Coupling[lambdaH2Inv1, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}])/2 + (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3HiggsPortal2CGInv2, {Index[d$$1, SU2L[fund]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3]}]*Coupling[lambdaH2Inv2, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}])/2 + (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3HiggsPortal2CGInv2, {Index[d$$1, SU2L[fund]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$2, T3BSMd3], Index[d$$1, T3BSMd3]}]*Coupling[lambdaH2Inv2, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}])/2 + Bar[Field[NewScalar1, Scalar, {}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3CrossScalarCGInv1, {Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3]}]*Coupling[lambda12Inv1, {}, 0]*Field[NewScalar1, Scalar, {}, {}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}] + (Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}]]*CG[T3Scalar2SelfCGInv1, {Index[d$$1, T3BSMd3], Index[d$$3, T3BSMd3], Index[d$$2, T3BSMd3], Index[d$$4, T3BSMd3]}]*Coupling[lambdaS2Inv1, {}, 0]*Field[NewScalar2, Scalar, {Index[d$$3, T3BSMd3]}, {}]*Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd3]}, {}])/2 + (Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}]]*CG[T3Scalar2SelfCGInv2, {Index[d$$1, T3BSMd3], Index[d$$3, T3BSMd3], Index[d$$2, T3BSMd3], Index[d$$4, T3BSMd3]}]*Coupling[lambdaS2Inv2, {}, 0]*Field[NewScalar2, Scalar, {Index[d$$3, T3BSMd3]}, {}]*Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd3]}, {}])/2 - Coupling[MF, {}, 0]*NCM[Bar[Field[NewFermion, Fermion, {Index[d$$1, T3BSMd2], Index[d$$1, NFlavor]}, {}]], Field[NewFermion, Fermion, {Index[d$$1, T3BSMd2], Index[d$$1, NFlavor]}, {}]] + Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3Y2CG, {Index[d$$1, SU2L[fund]], Bar[Index[d$$1, T3BSMd2]], Index[d$$1, T3BSMd3]}]*Coupling[y2, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]*NCM[Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]], DiracProduct[Proj[1]], Field[NewFermion, Fermion, {Index[d$$1, T3BSMd2], Index[d$$1, NFlavor]}, {}]] + CG[T3Y1CG, {Index[d$$1, SU2L[fund]], Index[d$$1, T3BSMd2]}]*Coupling[y1, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]*Field[NewScalar1, Scalar, {}, {}]*NCM[Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]], DiracProduct[GammaCC, Proj[1]], Transp[Bar[Field[NewFermion, Fermion, {Index[d$$1, T3BSMd2], Index[d$$1, NFlavor]}, {}]]]] + I*NCM[Bar[Field[NewFermion, Fermion, {Index[d$$1, T3BSMd2], Index[d$$1, NFlavor]}, {}]], DiracProduct[GammaM[Index[d$$1, Lorentz]]], Field[NewFermion, Fermion, {Index[d$$1, T3BSMd2], Index[d$$1, NFlavor]}, {Index[d$$1, Lorentz]}]] + Bar[Coupling[y2, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]]*CG[Bar[T3Y2CG], {Bar[Index[d$$1, SU2L[fund]]], Index[d$$1, T3BSMd2], Index[d$$1, T3BSMd3]}]*Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]*NCM[Bar[Field[NewFermion, Fermion, {Index[d$$1, T3BSMd2], Index[d$$1, NFlavor]}, {}]], DiracProduct[Proj[-1]], Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]] + Bar[Coupling[y1, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]]*Bar[Field[NewScalar1, Scalar, {}, {}]]*CG[Bar[T3Y1CG], {Bar[Index[d$$1, SU2L[fund]]], Bar[Index[d$$1, T3BSMd2]]}]*NCM[Transp[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]], DiracProduct[GammaCC, Proj[-1]], Field[NewFermion, Fermion, {Index[d$$1, T3BSMd2], Index[d$$1, NFlavor]}, {}]]
```

## One-loop renormalisable beta functions

RGBeta's report convention is `16*pi^2*dX/dln(mu)`. These are candidate native RGBeta expressions, not independently mapped/certified Matchete coefficients.

### UV: 20 beta functions

**beta_gY**

```text
(163*gY^3)/6
```

**beta_g2**

```text
-1/2*g2^3
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
(-9*g2^2*Matrix[ye][gen[$i], gen[$j]])/4 - (15*gY^2*Matrix[ye][gen[$i], gen[$j]])/4 + 3*Tr[yd . Trans[Bar[yd]]]*Matrix[ye][gen[$i], gen[$j]] + Tr[ye . Trans[Bar[ye]]]*Matrix[ye][gen[$i], gen[$j]] + 3*Tr[yu . Trans[Bar[yu]]]*Matrix[ye][gen[$i], gen[$j]] + (2*Matrix[y1, Trans[Bar[y1]], ye][gen[$i], gen[$j]] + 3*Matrix[y2, Trans[Bar[y2]], ye][gen[$i], gen[$j]] + 6*Matrix[ye, Trans[Bar[ye]], ye][gen[$i], gen[$j]])/4
```

**beta_y1**

```text
(-9*g2^2*Matrix[y1][gen[$i], heavy[$j]])/2 - (15*gY^2*Matrix[y1][gen[$i], heavy[$j]])/2 + 2*Tr[y1 . Trans[Bar[y1]]]*Matrix[y1][gen[$i], heavy[$j]] + Matrix[y1, Trans[Bar[y1]], y1][gen[$i], heavy[$j]] + (3*Matrix[y2, Trans[Bar[y2]], y1][gen[$i], heavy[$j]])/4 + Matrix[ye, Trans[Bar[ye]], y1][gen[$i], heavy[$j]]/2
```

**beta_y2**

```text
(-9*g2^2*Matrix[y2][gen[$i], heavy[$j]])/2 - (15*gY^2*Matrix[y2][gen[$i], heavy[$j]])/2 + Tr[y2 . Trans[Bar[y2]]]*Matrix[y2][gen[$i], heavy[$j]] + (Matrix[y1, Trans[Bar[y1]], y2][gen[$i], heavy[$j]] + 3*Matrix[y2, Trans[Bar[y2]], y2][gen[$i], heavy[$j]] + Matrix[ye, Trans[Bar[ye]], y2][gen[$i], heavy[$j]])/2
```

**beta_MF**

```text
(-9*g2^2*Matrix[MF][heavy[$i], heavy[$j]])/2 - (27*gY^2*Matrix[MF][heavy[$i], heavy[$j]])/2 + (3*Matrix[MF, Trans[Bar[y2]], y2][heavy[$i], heavy[$j]] + 2*Matrix[Trans[y1], Bar[y1], MF][heavy[$i], heavy[$j]])/4
```

**beta_mS1Sq**

```text
-24*gY^2*mS1Sq + 6*lambda12*mS2Sq + mS1Sq*(4*lambdaS1 + 4*Tr[y1 . Trans[Bar[y1]]]) - 8*Tr[MF . Trans[Bar[MF]] . Trans[y1] . Bar[y1]]
```

**beta_mS2Sq**

```text
2*lambda12*mS1Sq - 12*g2^2*mS2Sq - 6*gY^2*mS2Sq + mS2Sq*(8*lambdaS2 + 4*lambdaS2Adj + 2*Tr[y2 . Trans[Bar[y2]]]) - 4*Tr[MF . Trans[Bar[y2]] . y2 . Trans[Bar[MF]]]
```

**beta_lambdaH**

```text
(9*g2^4)/4 + (3*gY^4)/4 + gY^2*((3*g2^2)/2 - 3*lambdaH) - 9*g2^2*lambdaH + 12*lambdaH^2 + 2*lambdaH1^2 + 6*lambdaH2^2 + 3*lambdaH2Adj^2 + 4*lambdaT3*Bar[lambdaT3] + lambdaH*(12*Tr[yd . Trans[Bar[yd]]] + 4*Tr[ye . Trans[Bar[ye]]] + 12*Tr[yu . Trans[Bar[yu]]]) - 12*Tr[yd . Trans[Bar[yd]] . yd . Trans[Bar[yd]]] - 4*Tr[ye . Trans[Bar[ye]] . ye . Trans[Bar[ye]]] - 12*Tr[yu . Trans[Bar[yu]] . yu . Trans[Bar[yu]]]
```

**beta_lambdaS1**

```text
192*gY^4 + 6*lambda12^2 + 4*lambdaH1^2 - 48*gY^2*lambdaS1 + 10*lambdaS1^2 + 8*lambdaS1*Tr[y1 . Trans[Bar[y1]]] - 8*Tr[y1 . Trans[Bar[y1]] . y1 . Trans[Bar[y1]]]
```

**beta_lambdaS2**

```text
22*g2^4 + 12*gY^4 + 2*lambda12^2 + 4*lambdaH2^2 + gY^2*(8*g2^2 - 12*lambdaS2) - 24*g2^2*lambdaS2 + 14*lambdaS2^2 + 16*lambdaS2Adj^2 + lambdaS2*(8*lambdaS2Adj + 4*Tr[y2 . Trans[Bar[y2]]]) - (8*Tr[y2 . Trans[Bar[y2]] . y2 . Trans[Bar[y2]]])/3
```

**beta_lambdaH1**

```text
12*gY^4 - (9*g2^2*lambdaH1)/2 - (51*gY^2*lambdaH1)/2 + 6*lambdaH*lambdaH1 + 4*lambdaH1^2 + 6*lambda12*lambdaH2 + 4*lambdaH1*lambdaS1 + 12*lambdaT3*Bar[lambdaT3] + lambdaH1*(4*Tr[y1 . Trans[Bar[y1]]] + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]]) - 4*Tr[y1 . Trans[Bar[y1]] . ye . Trans[Bar[ye]]]
```

**beta_lambdaH2**

```text
6*g2^4 + 3*gY^4 + 2*lambda12*lambdaH1 - (33*g2^2*lambdaH2)/2 - (15*gY^2*lambdaH2)/2 + 6*lambdaH*lambdaH2 + 4*lambdaH2^2 + 6*lambdaH2Adj^2 + 8*lambdaH2*lambdaS2 + 4*lambdaT3*Bar[lambdaT3] + lambdaH2*(4*lambdaS2Adj + 2*Tr[y2 . Trans[Bar[y2]]] + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]]) - 2*Tr[y2 . Trans[Bar[y2]] . ye . Trans[Bar[ye]]]
```

**beta_lambda12**

```text
48*gY^4 - 12*g2^2*lambda12 - 30*gY^2*lambda12 + 4*lambda12^2 + 4*lambdaH1*lambdaH2 + 4*lambda12*lambdaS1 + 8*lambda12*lambdaS2 + 4*lambdaT3*Bar[lambdaT3] + lambda12*(4*lambdaS2Adj + 4*Tr[y1 . Trans[Bar[y1]]] + 2*Tr[y2 . Trans[Bar[y2]]]) - 4*Tr[y1 . Trans[Bar[y1]] . y2 . Trans[Bar[y2]]]
```

**beta_lambdaT3**

```text
(-21*g2^2*lambdaT3)/2 - (33*gY^2*lambdaT3)/2 + 2*lambda12*lambdaT3 + 2*lambdaH*lambdaT3 + 4*lambdaH1*lambdaT3 + 4*lambdaH2*lambdaT3 + lambdaT3*(4*lambdaH2Adj + 2*Tr[y1 . Trans[Bar[y1]]] + Tr[y2 . Trans[Bar[y2]]] + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaH2Adj**

```text
gY^2*(-4*g2^2 - (15*lambdaH2Adj)/2) - (33*g2^2*lambdaH2Adj)/2 + 2*lambdaH*lambdaH2Adj + 8*lambdaH2*lambdaH2Adj + 2*lambdaH2Adj*lambdaS2 + (8*lambdaT3*Bar[lambdaT3])/3 + lambdaH2Adj*(16*lambdaS2Adj + 2*Tr[y2 . Trans[Bar[y2]]] + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]]) + (4*Tr[y2 . Trans[Bar[y2]] . ye . Trans[Bar[ye]]])/3
```

**beta_lambdaS2Adj**

```text
-2*g2^4 + 3*lambdaH2Adj^2 + gY^2*(8*g2^2 - 12*lambdaS2Adj) - 24*g2^2*lambdaS2Adj + 12*lambdaS2*lambdaS2Adj + 20*lambdaS2Adj^2 + 4*lambdaS2Adj*Tr[y2 . Trans[Bar[y2]]] - (2*Tr[y2 . Trans[Bar[y2]] . y2 . Trans[Bar[y2]]])/3
```

### EFT1: 17 beta functions

**beta_gY**

```text
(55*gY^3)/6
```

**beta_g2**

```text
(-5*g2^3)/2
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
-24*gY^2*mS1Sq + 4*lambdaS1*mS1Sq + 6*lambda12*mS2Sq
```

**beta_mS2Sq**

```text
2*lambda12*mS1Sq - 12*g2^2*mS2Sq - 6*gY^2*mS2Sq + (8*lambdaS2 + 4*lambdaS2Adj)*mS2Sq
```

**beta_lambdaH**

```text
(9*g2^4)/4 + (3*gY^4)/4 + gY^2*((3*g2^2)/2 - 3*lambdaH) - 9*g2^2*lambdaH + 12*lambdaH^2 + 2*lambdaH1^2 + 6*lambdaH2^2 + 3*lambdaH2Adj^2 + 4*lambdaT3*Bar[lambdaT3] + lambdaH*(12*Tr[yd . Trans[Bar[yd]]] + 4*Tr[ye . Trans[Bar[ye]]] + 12*Tr[yu . Trans[Bar[yu]]]) - 12*Tr[yd . Trans[Bar[yd]] . yd . Trans[Bar[yd]]] - 4*Tr[ye . Trans[Bar[ye]] . ye . Trans[Bar[ye]]] - 12*Tr[yu . Trans[Bar[yu]] . yu . Trans[Bar[yu]]]
```

**beta_lambdaS1**

```text
192*gY^4 + 6*lambda12^2 + 4*lambdaH1^2 - 48*gY^2*lambdaS1 + 10*lambdaS1^2
```

**beta_lambdaS2**

```text
22*g2^4 + 12*gY^4 + 2*lambda12^2 + 4*lambdaH2^2 + gY^2*(8*g2^2 - 12*lambdaS2) - 24*g2^2*lambdaS2 + 14*lambdaS2^2 + 8*lambdaS2*lambdaS2Adj + 16*lambdaS2Adj^2
```

**beta_lambdaH1**

```text
12*gY^4 - (9*g2^2*lambdaH1)/2 - (51*gY^2*lambdaH1)/2 + 6*lambdaH*lambdaH1 + 4*lambdaH1^2 + 6*lambda12*lambdaH2 + 4*lambdaH1*lambdaS1 + 12*lambdaT3*Bar[lambdaT3] + lambdaH1*(6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaH2**

```text
6*g2^4 + 3*gY^4 + 2*lambda12*lambdaH1 - (33*g2^2*lambdaH2)/2 - (15*gY^2*lambdaH2)/2 + 6*lambdaH*lambdaH2 + 4*lambdaH2^2 + 6*lambdaH2Adj^2 + 8*lambdaH2*lambdaS2 + 4*lambdaT3*Bar[lambdaT3] + lambdaH2*(4*lambdaS2Adj + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambda12**

```text
48*gY^4 - 12*g2^2*lambda12 - 30*gY^2*lambda12 + 4*lambda12^2 + 4*lambdaH1*lambdaH2 + 4*lambda12*lambdaS1 + 8*lambda12*lambdaS2 + 4*lambda12*lambdaS2Adj + 4*lambdaT3*Bar[lambdaT3]
```

**beta_lambdaT3**

```text
(-21*g2^2*lambdaT3)/2 - (33*gY^2*lambdaT3)/2 + 2*lambda12*lambdaT3 + 2*lambdaH*lambdaT3 + 4*lambdaH1*lambdaT3 + 4*lambdaH2*lambdaT3 + lambdaT3*(4*lambdaH2Adj + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaH2Adj**

```text
gY^2*(-4*g2^2 - (15*lambdaH2Adj)/2) - (33*g2^2*lambdaH2Adj)/2 + 2*lambdaH*lambdaH2Adj + 8*lambdaH2*lambdaH2Adj + 2*lambdaH2Adj*lambdaS2 + (8*lambdaT3*Bar[lambdaT3])/3 + lambdaH2Adj*(16*lambdaS2Adj + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaS2Adj**

```text
-2*g2^4 + 3*lambdaH2Adj^2 + gY^2*(8*g2^2 - 12*lambdaS2Adj) - 24*g2^2*lambdaS2Adj + 12*lambdaS2*lambdaS2Adj + 20*lambdaS2Adj^2
```

## Interpretation and limitations

- UV RGEs include symmetry-allowed scalar quartics, bilinears, and the additional Yukawa where gauge invariance permits it.
- EFT1 RGEs describe the renormalisable SM+S1+S2 theory; the heavy fermion and its Yukawas have been removed.
- EFT1 higher-dimensional Wilson RGEs and final C5 have **not** been rederived with the full general interactions.
- Mixed scalar mass terms require mass-basis diagonalisation for B(-1) and C(-1) before general matching.
- The existing heavy-fermion flavour-mass parameterisation is not certified as the most general non-diagonal flavour mass basis.

# T3 E_-2: general-interaction symbolic results

Fields: S1=(3,-1), S2=(3,0), F=(2,-1/2) [SU(2) dimension, hypercharge].

Additional independent quartics: **6**; additional complex bilinears: **1**; additional Yukawa: **y2Conj**.

Heavy scalar mass basis: **CheckedNotMatched**. General sequential matching: **not performed**. Matchete/RGBeta physical coupling dictionary: **not certified**.

## Independent interaction inventory

- `Yukawa1_GenericCG`
- `Yukawa2_GenericCG`
- `Yukawa2ConjugateScalar_GenericCG`
- `Scalar1SelfInv1`
- `Scalar1SelfInv2`
- `Scalar2SelfInv1`
- `Scalar2SelfInv2`
- `HiggsPortal1Inv1`
- `HiggsPortal1Inv2`
- `HiggsPortal2Inv1`
- `HiggsPortal2Inv2`
- `ScalarCrossPortalInv1`
- `ScalarCrossPortalInv2`
- `ScalarCrossPortalInv3`
- `T3ScalarMix_GenericCG`
- `T3ExtraQ_0_0_0_0_0_4_1`
- `T3ExtraQ_0_0_0_0_1_3_1`
- `T3ExtraQ_0_0_1_1_0_2_1`
- `T3ExtraQ_0_0_1_1_0_2_2`
- `T3ExtraQ_0_2_0_1_0_1_1`
- `T3ExtraQ_1_1_0_0_0_2_1`
- `T3ExtraM_0_0_0_0_0_2_1`

## Complete BSM UV Lagrangian (native Matchete InputForm)

The expression below includes BSM free fields and interactions; it excludes the pure Standard Model Lagrangian.

```text
Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}]]*CG[T3ExtraCGM000002I1, {Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3]}]*Coupling[lambdaT3ExtraM000002I1, {}, 0] + Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$3, T3BSMd3]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd3]}, {}]]*CG[T3ExtraCGQ000004I1, {Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3], Index[d$$3, T3BSMd3], Index[d$$4, T3BSMd3]}]*Coupling[lambdaT3ExtraQ000004I1, {}, 0] + Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}]]*CG[T3ExtraCGQ020101I1, {Index[d$$1, SU2L[fund]], Index[d$$2, SU2L[fund]], Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3]}]*Coupling[lambdaT3ExtraQ020101I1, {}, 0] + Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}]]*CG[T3ExtraCGQ110002I1, {Bar[Index[d$$2, SU2L[fund]]], Index[d$$1, SU2L[fund]], Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3]}]*Coupling[lambdaT3ExtraQ110002I1, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}] - Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*Coupling[MS1, {}, 0]^2*Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}] + Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {Index[d$$1, Lorentz]}]]*Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {Index[d$$1, Lorentz]}] + (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3HiggsPortal1CGInv1, {Index[d$$1, SU2L[fund]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3]}]*Coupling[lambdaH1Inv1, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd3]}, {}])/2 + (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3HiggsPortal1CGInv1, {Index[d$$1, SU2L[fund]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$2, T3BSMd3], Index[d$$1, T3BSMd3]}]*Coupling[lambdaH1Inv1, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd3]}, {}])/2 + (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3HiggsPortal1CGInv2, {Index[d$$1, SU2L[fund]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3]}]*Coupling[lambdaH1Inv2, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd3]}, {}])/2 + (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3HiggsPortal1CGInv2, {Index[d$$1, SU2L[fund]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$2, T3BSMd3], Index[d$$1, T3BSMd3]}]*Coupling[lambdaH1Inv2, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd3]}, {}])/2 + Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3MixCG, {Bar[Index[d$$1, SU2L[fund]]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$2, T3BSMd3], Index[d$$1, T3BSMd3]}]*Coupling[lambdaT3, {}, 0]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd3]}, {}] + Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$3, T3BSMd3]}, {}]]*CG[T3ExtraCGQ001102I1, {Index[d$$4, T3BSMd3], Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3], Index[d$$3, T3BSMd3]}]*Coupling[lambdaT3ExtraQ001102I1, {}, 0]*Field[NewScalar1, Scalar, {Index[d$$4, T3BSMd3]}, {}] + Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$3, T3BSMd3]}, {}]]*CG[T3ExtraCGQ001102I2, {Index[d$$4, T3BSMd3], Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3], Index[d$$3, T3BSMd3]}]*Coupling[lambdaT3ExtraQ001102I2, {}, 0]*Field[NewScalar1, Scalar, {Index[d$$4, T3BSMd3]}, {}] + (Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*Bar[Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd3]}, {}]]*CG[T3Scalar1SelfCGInv1, {Index[d$$1, T3BSMd3], Index[d$$3, T3BSMd3], Index[d$$2, T3BSMd3], Index[d$$4, T3BSMd3]}]*Coupling[lambdaS1Inv1, {}, 0]*Field[NewScalar1, Scalar, {Index[d$$3, T3BSMd3]}, {}]*Field[NewScalar1, Scalar, {Index[d$$4, T3BSMd3]}, {}])/2 + (Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*Bar[Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd3]}, {}]]*CG[T3Scalar1SelfCGInv2, {Index[d$$1, T3BSMd3], Index[d$$3, T3BSMd3], Index[d$$2, T3BSMd3], Index[d$$4, T3BSMd3]}]*Coupling[lambdaS1Inv2, {}, 0]*Field[NewScalar1, Scalar, {Index[d$$3, T3BSMd3]}, {}]*Field[NewScalar1, Scalar, {Index[d$$4, T3BSMd3]}, {}])/2 - Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*Coupling[MS2, {}, 0]^2*Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}] + Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {Index[d$$1, Lorentz]}]]*Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {Index[d$$1, Lorentz]}] + Bar[Coupling[lambdaT3, {}, 0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]]*Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[Bar[T3MixCG], {Index[d$$1, SU2L[fund]], Index[d$$2, SU2L[fund]], Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3]}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}] + (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3HiggsPortal2CGInv1, {Index[d$$1, SU2L[fund]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3]}]*Coupling[lambdaH2Inv1, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}])/2 + (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3HiggsPortal2CGInv1, {Index[d$$1, SU2L[fund]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$2, T3BSMd3], Index[d$$1, T3BSMd3]}]*Coupling[lambdaH2Inv1, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}])/2 + (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3HiggsPortal2CGInv2, {Index[d$$1, SU2L[fund]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3]}]*Coupling[lambdaH2Inv2, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}])/2 + (Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3HiggsPortal2CGInv2, {Index[d$$1, SU2L[fund]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$2, T3BSMd3], Index[d$$1, T3BSMd3]}]*Coupling[lambdaH2Inv2, {}, 0]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}])/2 + Bar[Coupling[lambdaT3ExtraQ020101I1, {}, 0]]*CG[Bar[T3ExtraCGQ020101I1], {Bar[Index[d$$1, SU2L[fund]]], Bar[Index[d$$2, SU2L[fund]]], Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3]}]*Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}] + Bar[Coupling[lambdaT3ExtraM000002I1, {}, 0]]*CG[T3ExtraCGM000002I1, {Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3]}]*Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}] + Bar[Coupling[lambdaT3ExtraQ110002I1, {}, 0]]*Bar[Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]]*CG[T3ExtraCGQ110002I1, {Bar[Index[d$$2, SU2L[fund]]], Index[d$$1, SU2L[fund]], Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3]}]*Field[H, Scalar, {Index[d$$2, SU2L[fund]]}, {}]*Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}] + Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$3, T3BSMd3]}, {}]]*CG[T3ExtraCGQ000013I1, {Index[d$$4, T3BSMd3], Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3], Index[d$$3, T3BSMd3]}]*Coupling[lambdaT3ExtraQ000013I1, {}, 0]*Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd3]}, {}] + Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}]]*CG[T3CrossScalarCGInv1, {Index[d$$1, T3BSMd3], Index[d$$3, T3BSMd3], Index[d$$2, T3BSMd3], Index[d$$4, T3BSMd3]}]*Coupling[lambda12Inv1, {}, 0]*Field[NewScalar1, Scalar, {Index[d$$3, T3BSMd3]}, {}]*Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd3]}, {}] + Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}]]*CG[T3CrossScalarCGInv2, {Index[d$$1, T3BSMd3], Index[d$$3, T3BSMd3], Index[d$$2, T3BSMd3], Index[d$$4, T3BSMd3]}]*Coupling[lambda12Inv2, {}, 0]*Field[NewScalar1, Scalar, {Index[d$$3, T3BSMd3]}, {}]*Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd3]}, {}] + Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}]]*CG[T3CrossScalarCGInv3, {Index[d$$1, T3BSMd3], Index[d$$3, T3BSMd3], Index[d$$2, T3BSMd3], Index[d$$4, T3BSMd3]}]*Coupling[lambda12Inv3, {}, 0]*Field[NewScalar1, Scalar, {Index[d$$3, T3BSMd3]}, {}]*Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd3]}, {}] + (Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}]]*CG[T3Scalar2SelfCGInv1, {Index[d$$1, T3BSMd3], Index[d$$3, T3BSMd3], Index[d$$2, T3BSMd3], Index[d$$4, T3BSMd3]}]*Coupling[lambdaS2Inv1, {}, 0]*Field[NewScalar2, Scalar, {Index[d$$3, T3BSMd3]}, {}]*Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd3]}, {}])/2 + (Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*Bar[Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}]]*CG[T3Scalar2SelfCGInv2, {Index[d$$1, T3BSMd3], Index[d$$3, T3BSMd3], Index[d$$2, T3BSMd3], Index[d$$4, T3BSMd3]}]*Coupling[lambdaS2Inv2, {}, 0]*Field[NewScalar2, Scalar, {Index[d$$3, T3BSMd3]}, {}]*Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd3]}, {}])/2 + Bar[Coupling[lambdaT3ExtraQ001102I1, {}, 0]]*Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3ExtraCGQ001102I1, {Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3], Index[d$$3, T3BSMd3], Index[d$$4, T3BSMd3]}]*Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd3]}, {}]*Field[NewScalar2, Scalar, {Index[d$$3, T3BSMd3]}, {}]*Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd3]}, {}] + Bar[Coupling[lambdaT3ExtraQ001102I2, {}, 0]]*Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3ExtraCGQ001102I2, {Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3], Index[d$$3, T3BSMd3], Index[d$$4, T3BSMd3]}]*Field[NewScalar1, Scalar, {Index[d$$2, T3BSMd3]}, {}]*Field[NewScalar2, Scalar, {Index[d$$3, T3BSMd3]}, {}]*Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd3]}, {}] + Bar[Coupling[lambdaT3ExtraQ000013I1, {}, 0]]*Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3ExtraCGQ000013I1, {Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3], Index[d$$3, T3BSMd3], Index[d$$4, T3BSMd3]}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}]*Field[NewScalar2, Scalar, {Index[d$$3, T3BSMd3]}, {}]*Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd3]}, {}] + Bar[Coupling[lambdaT3ExtraQ000004I1, {}, 0]]*CG[T3ExtraCGQ000004I1, {Index[d$$1, T3BSMd3], Index[d$$2, T3BSMd3], Index[d$$3, T3BSMd3], Index[d$$4, T3BSMd3]}]*Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]*Field[NewScalar2, Scalar, {Index[d$$2, T3BSMd3]}, {}]*Field[NewScalar2, Scalar, {Index[d$$3, T3BSMd3]}, {}]*Field[NewScalar2, Scalar, {Index[d$$4, T3BSMd3]}, {}] - Coupling[MF, {}, 0]*NCM[Bar[Field[NewFermion, Fermion, {Index[d$$1, T3BSMd2], Index[d$$1, NFlavor]}, {}]], Field[NewFermion, Fermion, {Index[d$$1, T3BSMd2], Index[d$$1, NFlavor]}, {}]] + Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[T3Y2CG, {Index[d$$1, SU2L[fund]], Bar[Index[d$$1, T3BSMd2]], Index[d$$1, T3BSMd3]}]*Coupling[y2, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]*NCM[Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]], DiracProduct[Proj[1]], Field[NewFermion, Fermion, {Index[d$$1, T3BSMd2], Index[d$$1, NFlavor]}, {}]] + CG[T3Y2ConjCG, {Index[d$$1, SU2L[fund]], Bar[Index[d$$1, T3BSMd2]], Index[d$$1, T3BSMd3]}]*Coupling[y2Conj, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]*Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]*NCM[Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]], DiracProduct[Proj[1]], Field[NewFermion, Fermion, {Index[d$$1, T3BSMd2], Index[d$$1, NFlavor]}, {}]] + CG[T3Y1CG, {Index[d$$1, SU2L[fund]], Index[d$$1, T3BSMd2], Index[d$$1, T3BSMd3]}]*Coupling[y1, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]*Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]*NCM[Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]], DiracProduct[GammaCC, Proj[1]], Transp[Bar[Field[NewFermion, Fermion, {Index[d$$1, T3BSMd2], Index[d$$1, NFlavor]}, {}]]]] + I*NCM[Bar[Field[NewFermion, Fermion, {Index[d$$1, T3BSMd2], Index[d$$1, NFlavor]}, {}]], DiracProduct[GammaM[Index[d$$1, Lorentz]]], Field[NewFermion, Fermion, {Index[d$$1, T3BSMd2], Index[d$$1, NFlavor]}, {Index[d$$1, Lorentz]}]] + Bar[Coupling[y2Conj, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]]*Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[Bar[T3Y2ConjCG], {Bar[Index[d$$1, SU2L[fund]]], Index[d$$1, T3BSMd2], Index[d$$1, T3BSMd3]}]*NCM[Bar[Field[NewFermion, Fermion, {Index[d$$1, T3BSMd2], Index[d$$1, NFlavor]}, {}]], DiracProduct[Proj[-1]], Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]] + Bar[Coupling[y2, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]]*CG[Bar[T3Y2CG], {Bar[Index[d$$1, SU2L[fund]]], Index[d$$1, T3BSMd2], Index[d$$1, T3BSMd3]}]*Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]*NCM[Bar[Field[NewFermion, Fermion, {Index[d$$1, T3BSMd2], Index[d$$1, NFlavor]}, {}]], DiracProduct[Proj[-1]], Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]] + Bar[Coupling[y1, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]]*Bar[Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*CG[Bar[T3Y1CG], {Bar[Index[d$$1, SU2L[fund]]], Bar[Index[d$$1, T3BSMd2]], Index[d$$1, T3BSMd3]}]*NCM[Transp[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]], DiracProduct[GammaCC, Proj[-1]], Field[NewFermion, Fermion, {Index[d$$1, T3BSMd2], Index[d$$1, NFlavor]}, {}]]
```

## One-loop renormalisable beta functions

RGBeta's report convention is `16*pi^2*dX/dln(mu)`. These are candidate native RGBeta expressions, not independently mapped/certified Matchete coefficients.

### UV: 32 beta functions

**beta_gY**

```text
(59*gY^3)/6
```

**beta_g2**

```text
g2^3/6
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
(-9*g2^2*Matrix[ye][gen[$i], gen[$j]])/4 - (15*gY^2*Matrix[ye][gen[$i], gen[$j]])/4 + 3*Tr[yd . Trans[Bar[yd]]]*Matrix[ye][gen[$i], gen[$j]] + Tr[ye . Trans[Bar[ye]]]*Matrix[ye][gen[$i], gen[$j]] + 3*Tr[yu . Trans[Bar[yu]]]*Matrix[ye][gen[$i], gen[$j]] + (3*(Matrix[y1, Trans[Bar[y1]], ye][gen[$i], gen[$j]] + Matrix[y2, Trans[Bar[y2]], ye][gen[$i], gen[$j]] + Matrix[y2Conj, Trans[Bar[y2Conj]], ye][gen[$i], gen[$j]] + 2*Matrix[ye, Trans[Bar[ye]], ye][gen[$i], gen[$j]]))/4
```

**beta_y1**

```text
(-9*g2^2*Matrix[y1][gen[$i], heavy[$j]])/2 - (3*gY^2*Matrix[y1][gen[$i], heavy[$j]])/2 + Tr[y1 . Trans[Bar[y1]]]*Matrix[y1][gen[$i], heavy[$j]] + (6*Matrix[y1, Trans[Bar[y1]], y1][gen[$i], heavy[$j]] + 3*Matrix[y2, Trans[Bar[y2]], y1][gen[$i], heavy[$j]] + 3*Matrix[y2Conj, Trans[Bar[y2Conj]], y1][gen[$i], heavy[$j]] + 2*Matrix[ye, Trans[Bar[ye]], y1][gen[$i], heavy[$j]])/4
```

**beta_y2**

```text
(-9*g2^2*Matrix[y2][gen[$i], heavy[$j]])/2 - (3*gY^2*Matrix[y2][gen[$i], heavy[$j]])/2 + Tr[y2 . Trans[Bar[y2]]]*Matrix[y2][gen[$i], heavy[$j]] + Tr[y2Conj . Trans[Bar[y2Conj]]]*Matrix[y2][gen[$i], heavy[$j]] + 2*Tr[y2 . Trans[Bar[y2Conj]]]*Matrix[y2Conj][gen[$i], heavy[$j]] + (9*Matrix[y1, Trans[Bar[y1]], y2][gen[$i], heavy[$j]] + 18*Matrix[y2, Trans[Bar[y2]], y2][gen[$i], heavy[$j]] + 29*Matrix[y2, Trans[Bar[y2Conj]], y2Conj][gen[$i], heavy[$j]] + 29*Matrix[y2Conj, Trans[Bar[y2Conj]], y2][gen[$i], heavy[$j]] + 6*Matrix[ye, Trans[Bar[ye]], y2][gen[$i], heavy[$j]])/12
```

**beta_MF**

```text
(-9*g2^2*Matrix[MF][heavy[$i], heavy[$j]])/2 - (3*gY^2*Matrix[MF][heavy[$i], heavy[$j]])/2 + (3*(Matrix[MF, Trans[Bar[y2]], y2][heavy[$i], heavy[$j]] + Matrix[MF, Trans[Bar[y2Conj]], y2Conj][heavy[$i], heavy[$j]] + Matrix[Trans[y1], Bar[y1], MF][heavy[$i], heavy[$j]]))/4
```

**beta_mS1Sq**

```text
-12*g2^2*mS1Sq - 6*gY^2*mS1Sq + (6*lambda12 + 4*lambda12Cross)*mS2Sq + 4*muT3GeneralS2Sq*(3*Bar[lambdaT3GeneralQ001102I1] + Bar[lambdaT3GeneralQ001102I2]) + 12*lambdaT3GeneralQ001102I1*Bar[muT3GeneralS2Sq] + 4*lambdaT3GeneralQ001102I2*Bar[muT3GeneralS2Sq] + mS1Sq*(8*lambdaS1 + 4*lambdaS1Adj + 2*Tr[y1 . Trans[Bar[y1]]]) - 4*Tr[MF . Trans[Bar[MF]] . Trans[y1] . Bar[y1]]
```

**beta_mS2Sq**

```text
(6*lambda12 + 4*lambda12Cross)*mS1Sq - 12*g2^2*mS2Sq + 20*lambdaT3GeneralQ000013I1*Bar[muT3GeneralS2Sq] + 4*Bar[muT3GeneralS2Sq]*Tr[y2 . Trans[Bar[y2Conj]]] + muT3GeneralS2Sq*(20*Bar[lambdaT3GeneralQ000013I1] + 4*Tr[y2Conj . Trans[Bar[y2]]]) + mS2Sq*(8*lambdaS2 + 4*lambdaS2Adj + 2*Tr[y2 . Trans[Bar[y2]]] + 2*Tr[y2Conj . Trans[Bar[y2Conj]]]) - 4*Tr[MF . Trans[Bar[y2]] . y2 . Trans[Bar[MF]]] - 4*Tr[MF . Trans[Bar[y2Conj]] . y2Conj . Trans[Bar[MF]]]
```

**beta_lambdaH**

```text
(9*g2^4)/4 + (3*gY^4)/4 + gY^2*((3*g2^2)/2 - 3*lambdaH) - 9*g2^2*lambdaH + 12*lambdaH^2 + 6*lambdaH1^2 + 3*lambdaH1Adj^2 + 6*lambdaH2^2 + 3*lambdaH2Adj^2 + 4*lambdaT3*Bar[lambdaT3] + 4*lambdaT3GeneralQ020101I1*Bar[lambdaT3GeneralQ020101I1] + 24*lambdaT3GeneralQ110002I1*Bar[lambdaT3GeneralQ110002I1] + lambdaH*(12*Tr[yd . Trans[Bar[yd]]] + 4*Tr[ye . Trans[Bar[ye]]] + 12*Tr[yu . Trans[Bar[yu]]]) - 12*Tr[yd . Trans[Bar[yd]] . yd . Trans[Bar[yd]]] - 4*Tr[ye . Trans[Bar[ye]] . ye . Trans[Bar[ye]]] - 12*Tr[yu . Trans[Bar[yu]] . yu . Trans[Bar[yu]]]
```

**beta_lambdaS1**

```text
22*g2^4 + 12*gY^4 + 6*lambda12^2 + 8*lambda12*lambda12Cross + (20*lambda12Cross^2)/3 + 4*lambdaH1^2 + gY^2*(8*g2^2 - 12*lambdaS1) - 24*g2^2*lambdaS1 + 14*lambdaS1^2 + 16*lambdaS1Adj^2 + 8*lambdaT3GeneralQ001102I1*(3*Bar[lambdaT3GeneralQ001102I1] + Bar[lambdaT3GeneralQ001102I2]) + lambdaT3GeneralQ001102I2*(8*Bar[lambdaT3GeneralQ001102I1] + (20*Bar[lambdaT3GeneralQ001102I2])/3) + lambdaS1*(8*lambdaS1Adj + 4*Tr[y1 . Trans[Bar[y1]]]) - (8*Tr[y1 . Trans[Bar[y1]] . y1 . Trans[Bar[y1]]])/3
```

**beta_lambdaS2**

```text
22*g2^4 + 6*lambda12^2 + 8*lambda12*lambda12Cross + (20*lambda12Cross^2)/3 + 4*lambdaH2^2 - 24*g2^2*lambdaS2 + 14*lambdaS2^2 + 16*lambdaS2Adj^2 + (832*lambdaT3GeneralQ000004I1*Bar[lambdaT3GeneralQ000004I1])/3 + (8*lambdaT3GeneralQ001102I1*(3*Bar[lambdaT3GeneralQ001102I1] + Bar[lambdaT3GeneralQ001102I2]))/3 + lambdaT3GeneralQ001102I2*((8*Bar[lambdaT3GeneralQ001102I1])/3 + 4*Bar[lambdaT3GeneralQ001102I2]) + (16*lambdaT3GeneralQ110002I1*Bar[lambdaT3GeneralQ110002I1])/3 + (32*Bar[lambdaT3GeneralQ000013I1]*Tr[y2 . Trans[Bar[y2Conj]]])/3 + lambdaT3GeneralQ000013I1*((352*Bar[lambdaT3GeneralQ000013I1])/3 + (32*Tr[y2Conj . Trans[Bar[y2]]])/3) + lambdaS2*(8*lambdaS2Adj + 4*Tr[y2 . Trans[Bar[y2]]] + 4*Tr[y2Conj . Trans[Bar[y2Conj]]]) - (8*Tr[y2 . Trans[Bar[y2]] . y2 . Trans[Bar[y2]]])/3 - (248*Tr[y2 . Trans[Bar[y2]] . y2Conj . Trans[Bar[y2Conj]]])/45 - (248*Tr[y2 . Trans[Bar[y2Conj]] . y2Conj . Trans[Bar[y2]]])/45 - (8*Tr[y2Conj . Trans[Bar[y2Conj]] . y2Conj . Trans[Bar[y2Conj]]])/3
```

**beta_lambdaH1**

```text
6*g2^4 + 3*gY^4 - (33*g2^2*lambdaH1)/2 - (15*gY^2*lambdaH1)/2 + 6*lambdaH*lambdaH1 + 4*lambdaH1^2 + 6*lambdaH1Adj^2 + (6*lambda12 + 4*lambda12Cross)*lambdaH2 + 8*lambdaH1*lambdaS1 + 4*lambdaT3*Bar[lambdaT3] + 4*lambdaT3GeneralQ110002I1*(3*Bar[lambdaT3GeneralQ001102I1] + Bar[lambdaT3GeneralQ001102I2]) + 4*lambdaT3GeneralQ020101I1*Bar[lambdaT3GeneralQ020101I1] + 12*lambdaT3GeneralQ001102I1*Bar[lambdaT3GeneralQ110002I1] + 4*lambdaT3GeneralQ001102I2*Bar[lambdaT3GeneralQ110002I1] + lambdaH1*(4*lambdaS1Adj + 2*Tr[y1 . Trans[Bar[y1]]] + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]]) - 2*Tr[y1 . Trans[Bar[y1]] . ye . Trans[Bar[ye]]]
```

**beta_lambdaH2**

```text
6*g2^4 + (6*lambda12 + 4*lambda12Cross)*lambdaH1 - (33*g2^2*lambdaH2)/2 - (3*gY^2*lambdaH2)/2 + 6*lambdaH*lambdaH2 + 4*lambdaH2^2 + 6*lambdaH2Adj^2 + 8*lambdaH2*lambdaS2 + 4*lambdaT3*Bar[lambdaT3] + 4*lambdaT3GeneralQ020101I1*Bar[lambdaT3GeneralQ020101I1] + 20*lambdaT3GeneralQ000013I1*Bar[lambdaT3GeneralQ110002I1] + 4*Bar[lambdaT3GeneralQ110002I1]*Tr[y2 . Trans[Bar[y2Conj]]] + lambdaT3GeneralQ110002I1*(4*(5*Bar[lambdaT3GeneralQ000013I1] + 4*Bar[lambdaT3GeneralQ110002I1]) + 4*Tr[y2Conj . Trans[Bar[y2]]]) + lambdaH2*(4*lambdaS2Adj + 2*Tr[y2 . Trans[Bar[y2]]] + 2*Tr[y2Conj . Trans[Bar[y2Conj]]] + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]]) - 2*Tr[y2 . Trans[Bar[y2]] . ye . Trans[Bar[ye]]] - 2*Tr[y2Conj . Trans[Bar[y2Conj]] . ye . Trans[Bar[ye]]]
```

**beta_lambda12**

```text
(63*g2^4)/5 - 24*g2^2*lambda12 - 6*gY^2*lambda12 + 4*lambda12^2 + 16*lambda12Adj^2 + 4*lambda12Cross^2 + 4*lambdaH1*lambdaH2 + (8*lambda12 + 4*lambda12Cross)*lambdaS1 + 4*lambda12Cross*lambdaS1Adj + (8*lambda12 + 4*lambda12Cross)*lambdaS2 + 4*lambda12Cross*lambdaS2Adj + 2*lambdaT3*Bar[lambdaT3] + 4*lambdaT3GeneralQ001102I2*(Bar[lambdaT3GeneralQ000013I1] + Bar[lambdaT3GeneralQ001102I2]) + 4*lambdaT3GeneralQ000013I1*(5*Bar[lambdaT3GeneralQ001102I1] + Bar[lambdaT3GeneralQ001102I2]) + 2*lambdaT3GeneralQ020101I1*Bar[lambdaT3GeneralQ020101I1] + 4*Bar[lambdaT3GeneralQ001102I1]*Tr[y2 . Trans[Bar[y2Conj]]] + lambdaT3GeneralQ001102I1*(4*(5*Bar[lambdaT3GeneralQ000013I1] + 4*Bar[lambdaT3GeneralQ001102I1]) + 4*Tr[y2Conj . Trans[Bar[y2]]]) + lambda12*(4*lambdaS1Adj + 4*lambdaS2Adj + 2*Tr[y1 . Trans[Bar[y1]]] + 2*Tr[y2 . Trans[Bar[y2]]] + 2*Tr[y2Conj . Trans[Bar[y2Conj]]]) - (22*Tr[y1 . Trans[Bar[y1]] . y2 . Trans[Bar[y2]]])/15 - (22*Tr[y1 . Trans[Bar[y1]] . y2Conj . Trans[Bar[y2Conj]]])/15
```

**beta_lambdaT3**

```text
(-33*g2^2*lambdaT3)/2 - (9*gY^2*lambdaT3)/2 + 2*lambda12*lambdaT3 + 2*lambdaH*lambdaT3 + 4*lambdaH1*lambdaT3 + 4*lambdaH2*lambdaT3 + 4*lambdaT3GeneralQ001102I1*Bar[lambdaT3GeneralQ020101I1] - 2*lambdaT3GeneralQ001102I2*Bar[lambdaT3GeneralQ020101I1] + 8*lambdaT3GeneralQ110002I1*Bar[lambdaT3GeneralQ020101I1] + 2*Bar[lambdaT3GeneralQ020101I1]*Tr[y2 . Trans[Bar[y2Conj]]] + lambdaT3*((-4*lambda12Adj)/3 - 2*lambda12Cross - 2*lambdaH1Adj - (4*lambdaH2Adj)/3 + Tr[y1 . Trans[Bar[y1]]] + Tr[y2 . Trans[Bar[y2]]] + Tr[y2Conj . Trans[Bar[y2Conj]]] + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaH1Adj**

```text
gY^2*(-4*g2^2 - (15*lambdaH1Adj)/2) - (33*g2^2*lambdaH1Adj)/2 + 2*lambdaH*lambdaH1Adj + 8*lambdaH1*lambdaH1Adj + 12*lambda12Adj*lambdaH2Adj + 2*lambdaH1Adj*lambdaS1 - (4*lambdaT3*Bar[lambdaT3])/3 - (4*lambdaT3GeneralQ020101I1*Bar[lambdaT3GeneralQ020101I1])/3 + lambdaH1Adj*(16*lambdaS1Adj + 2*Tr[y1 . Trans[Bar[y1]]] + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]]) - (4*Tr[y1 . Trans[Bar[y1]] . ye . Trans[Bar[ye]]])/3
```

**beta_lambdaH2Adj**

```text
12*lambda12Adj*lambdaH1Adj - (33*g2^2*lambdaH2Adj)/2 - (3*gY^2*lambdaH2Adj)/2 + 2*lambdaH*lambdaH2Adj + 8*lambdaH2*lambdaH2Adj + 2*lambdaH2Adj*lambdaS2 - (8*lambdaT3*Bar[lambdaT3])/9 - (4*lambdaT3GeneralQ020101I1*Bar[lambdaT3GeneralQ020101I1])/3 + lambdaH2Adj*(16*lambdaS2Adj + 2*Tr[y2 . Trans[Bar[y2]]] + 2*Tr[y2Conj . Trans[Bar[y2Conj]]] + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]]) + (4*Tr[y2 . Trans[Bar[y2]] . ye . Trans[Bar[ye]]])/3 - (4*Tr[y2Conj . Trans[Bar[y2Conj]] . ye . Trans[Bar[ye]]])/3
```

**beta_lambdaS1Adj**

```text
-2*g2^4 + 12*lambda12Adj^2 - (4*lambda12Cross^2)/3 + 3*lambdaH1Adj^2 + gY^2*(8*g2^2 - 12*lambdaS1Adj) - 24*g2^2*lambdaS1Adj + 12*lambdaS1*lambdaS1Adj + 20*lambdaS1Adj^2 - (4*lambdaT3GeneralQ001102I2*Bar[lambdaT3GeneralQ001102I2])/3 + 4*lambdaS1Adj*Tr[y1 . Trans[Bar[y1]]] - (2*Tr[y1 . Trans[Bar[y1]] . y1 . Trans[Bar[y1]]])/3
```

**beta_lambdaS2Adj**

```text
-2*g2^4 + 12*lambda12Adj^2 - (4*lambda12Cross^2)/3 + 3*lambdaH2Adj^2 - 24*g2^2*lambdaS2Adj + 12*lambdaS2*lambdaS2Adj + 20*lambdaS2Adj^2 - (224*lambdaT3GeneralQ000004I1*Bar[lambdaT3GeneralQ000004I1])/3 - (4*lambdaT3GeneralQ001102I2*Bar[lambdaT3GeneralQ001102I1])/3 - (4*lambdaT3GeneralQ001102I1*(3*Bar[lambdaT3GeneralQ001102I1] + Bar[lambdaT3GeneralQ001102I2]))/3 - (8*lambdaT3GeneralQ110002I1*Bar[lambdaT3GeneralQ110002I1])/3 - (4*Bar[lambdaT3GeneralQ000013I1]*Tr[y2 . Trans[Bar[y2Conj]]])/3 + lambdaT3GeneralQ000013I1*((-44*Bar[lambdaT3GeneralQ000013I1])/3 - (4*Tr[y2Conj . Trans[Bar[y2]]])/3) + lambdaS2Adj*(4*Tr[y2 . Trans[Bar[y2]]] + 4*Tr[y2Conj . Trans[Bar[y2Conj]]]) - (2*Tr[y2 . Trans[Bar[y2]] . y2 . Trans[Bar[y2]]])/3 + (76*Tr[y2 . Trans[Bar[y2]] . y2Conj . Trans[Bar[y2Conj]]])/45 + (76*Tr[y2 . Trans[Bar[y2Conj]] . y2Conj . Trans[Bar[y2]]])/45 - (2*Tr[y2Conj . Trans[Bar[y2Conj]] . y2Conj . Trans[Bar[y2Conj]]])/3
```

**beta_lambda12Adj**

```text
-24*g2^2*lambda12Adj - 6*gY^2*lambda12Adj + 8*lambda12*lambda12Adj + 3*lambdaH1Adj*lambdaH2Adj + 2*lambda12Adj*lambdaS1 + 16*lambda12Adj*lambdaS1Adj + 2*lambda12Adj*lambdaS2 + 16*lambda12Adj*lambdaS2Adj - (2*lambdaT3*Bar[lambdaT3])/9 - (lambdaT3GeneralQ020101I1*Bar[lambdaT3GeneralQ020101I1])/3 + lambda12Adj*((16*lambda12Cross)/3 + 2*Tr[y1 . Trans[Bar[y1]]] + 2*Tr[y2 . Trans[Bar[y2]]] + 2*Tr[y2Conj . Trans[Bar[y2Conj]]]) + (2*Tr[y1 . Trans[Bar[y1]] . y2 . Trans[Bar[y2]]])/3 - (2*Tr[y1 . Trans[Bar[y1]] . y2Conj . Trans[Bar[y2Conj]]])/3
```

**beta_lambda12Cross**

```text
(51*g2^4)/10 - 24*g2^2*lambda12Cross - 6*gY^2*lambda12Cross + 8*lambda12*lambda12Cross + 10*lambda12Cross^2 + 2*lambda12Cross*lambdaS1 - 2*lambda12Cross*lambdaS1Adj + 2*lambda12Cross*lambdaS2 - 2*lambda12Cross*lambdaS2Adj - lambdaT3*Bar[lambdaT3] + 4*lambdaT3GeneralQ000013I1*Bar[lambdaT3GeneralQ001102I2] + 8*lambdaT3GeneralQ001102I1*Bar[lambdaT3GeneralQ001102I2] - lambdaT3GeneralQ020101I1*Bar[lambdaT3GeneralQ020101I1] + 2*Bar[lambdaT3GeneralQ001102I2]*Tr[y2 . Trans[Bar[y2Conj]]] + lambdaT3GeneralQ001102I2*(4*Bar[lambdaT3GeneralQ000013I1] + 8*Bar[lambdaT3GeneralQ001102I1] + 10*Bar[lambdaT3GeneralQ001102I2] + 2*Tr[y2Conj . Trans[Bar[y2]]]) + lambda12Cross*(2*Tr[y1 . Trans[Bar[y1]]] + 2*Tr[y2 . Trans[Bar[y2]]] + 2*Tr[y2Conj . Trans[Bar[y2Conj]]]) - (4*Tr[y1 . Trans[Bar[y1]] . y2 . Trans[Bar[y2]]])/5 - (4*Tr[y1 . Trans[Bar[y1]] . y2Conj . Trans[Bar[y2Conj]]])/5
```

**beta_y2Conj**

```text
2*Tr[y2Conj . Trans[Bar[y2]]]*Matrix[y2][gen[$i], heavy[$j]] - (9*g2^2*Matrix[y2Conj][gen[$i], heavy[$j]])/2 - (3*gY^2*Matrix[y2Conj][gen[$i], heavy[$j]])/2 + Tr[y2 . Trans[Bar[y2]]]*Matrix[y2Conj][gen[$i], heavy[$j]] + Tr[y2Conj . Trans[Bar[y2Conj]]]*Matrix[y2Conj][gen[$i], heavy[$j]] + (9*Matrix[y1, Trans[Bar[y1]], y2Conj][gen[$i], heavy[$j]] + 29*Matrix[y2, Trans[Bar[y2]], y2Conj][gen[$i], heavy[$j]] + 29*Matrix[y2Conj, Trans[Bar[y2]], y2][gen[$i], heavy[$j]] + 18*Matrix[y2Conj, Trans[Bar[y2Conj]], y2Conj][gen[$i], heavy[$j]] + 6*Matrix[ye, Trans[Bar[ye]], y2Conj][gen[$i], heavy[$j]])/12
```

**beta_lambdaT3GeneralQ000004I1**

```text
-24*g2^2*lambdaT3GeneralQ000004I1 + 12*lambdaS2*lambdaT3GeneralQ000004I1 - 8*lambdaS2Adj*lambdaT3GeneralQ000004I1 + 11*lambdaT3GeneralQ000013I1^2 + 3*lambdaT3GeneralQ001102I1^2 + 2*lambdaT3GeneralQ001102I1*lambdaT3GeneralQ001102I2 + lambdaT3GeneralQ001102I2^2 + 2*lambdaT3GeneralQ110002I1^2 + 2*lambdaT3GeneralQ000013I1*Tr[y2 . Trans[Bar[y2Conj]]] + lambdaT3GeneralQ000004I1*(4*Tr[y2 . Trans[Bar[y2]]] + 4*Tr[y2Conj . Trans[Bar[y2Conj]]]) - (23*Tr[y2 . Trans[Bar[y2Conj]] . y2 . Trans[Bar[y2Conj]]])/15
```

**beta_lambdaT3GeneralQ000013I1**

```text
-24*g2^2*lambdaT3GeneralQ000013I1 + lambda12*(6*lambdaT3GeneralQ001102I1 + 2*lambdaT3GeneralQ001102I2) + lambda12Cross*(4*lambdaT3GeneralQ001102I1 + 4*lambdaT3GeneralQ001102I2) + 4*lambdaH2*lambdaT3GeneralQ110002I1 + lambdaS2*(22*lambdaT3GeneralQ000013I1 + 2*Tr[y2 . Trans[Bar[y2Conj]]]) + lambdaT3GeneralQ000004I1*(88*Bar[lambdaT3GeneralQ000013I1] + 8*Tr[y2Conj . Trans[Bar[y2]]]) + lambdaT3GeneralQ000013I1*(4*Tr[y2 . Trans[Bar[y2]]] + 4*Tr[y2Conj . Trans[Bar[y2Conj]]]) - (46*Tr[y2 . Trans[Bar[y2]] . y2 . Trans[Bar[y2Conj]]])/15 - (46*Tr[y2 . Trans[Bar[y2Conj]] . y2Conj . Trans[Bar[y2Conj]]])/15
```

**beta_lambdaT3GeneralQ001102I1**

```text
-24*g2^2*lambdaT3GeneralQ001102I1 - 6*gY^2*lambdaT3GeneralQ001102I1 + 2*lambdaS2*lambdaT3GeneralQ001102I1 + lambdaS2Adj*(-4*lambdaT3GeneralQ001102I1 - 2*lambdaT3GeneralQ001102I2) + lambdaS1Adj*(4*lambdaT3GeneralQ001102I1 + 2*lambdaT3GeneralQ001102I2) + lambdaS1*(8*lambdaT3GeneralQ001102I1 + 2*lambdaT3GeneralQ001102I2) + lambda12Cross*(4*lambdaT3GeneralQ000013I1 + 4*lambdaT3GeneralQ001102I2) + 2*lambdaT3*lambdaT3GeneralQ020101I1 + 4*lambdaH1*lambdaT3GeneralQ110002I1 + 8*lambdaT3GeneralQ000004I1*(5*Bar[lambdaT3GeneralQ001102I1] + Bar[lambdaT3GeneralQ001102I2]) + lambda12*(10*lambdaT3GeneralQ000013I1 + 8*lambdaT3GeneralQ001102I1 + 2*Tr[y2 . Trans[Bar[y2Conj]]]) + lambdaT3GeneralQ001102I1*(2*Tr[y1 . Trans[Bar[y1]]] + 2*Tr[y2 . Trans[Bar[y2]]] + 2*Tr[y2Conj . Trans[Bar[y2Conj]]]) - (22*Tr[y1 . Trans[Bar[y1]] . y2 . Trans[Bar[y2Conj]]])/15
```

**beta_lambdaT3GeneralQ001102I2**

```text
-24*g2^2*lambdaT3GeneralQ001102I2 - 6*gY^2*lambdaT3GeneralQ001102I2 + 8*lambda12*lambdaT3GeneralQ001102I2 + 2*lambdaS1*lambdaT3GeneralQ001102I2 - 2*lambdaS1Adj*lambdaT3GeneralQ001102I2 + 2*lambdaS2*lambdaT3GeneralQ001102I2 + 2*lambdaS2Adj*lambdaT3GeneralQ001102I2 - 2*lambdaT3*lambdaT3GeneralQ020101I1 + 16*lambdaT3GeneralQ000004I1*Bar[lambdaT3GeneralQ001102I2] + lambda12Cross*(8*lambdaT3GeneralQ000013I1 + 16*lambdaT3GeneralQ001102I1 + 20*lambdaT3GeneralQ001102I2 + 4*Tr[y2 . Trans[Bar[y2Conj]]]) + lambdaT3GeneralQ001102I2*(2*Tr[y1 . Trans[Bar[y1]]] + 2*Tr[y2 . Trans[Bar[y2]]] + 2*Tr[y2Conj . Trans[Bar[y2Conj]]]) - (8*Tr[y1 . Trans[Bar[y1]] . y2 . Trans[Bar[y2Conj]]])/5
```

**beta_lambdaT3GeneralQ020101I1**

```text
(-33*g2^2*lambdaT3GeneralQ020101I1)/2 - (9*gY^2*lambdaT3GeneralQ020101I1)/2 + 2*lambda12*lambdaT3GeneralQ020101I1 - 2*lambda12Adj*lambdaT3GeneralQ020101I1 - 2*lambda12Cross*lambdaT3GeneralQ020101I1 + 2*lambdaH*lambdaT3GeneralQ020101I1 + 4*lambdaH1*lambdaT3GeneralQ020101I1 - 2*lambdaH1Adj*lambdaT3GeneralQ020101I1 + 4*lambdaH2*lambdaT3GeneralQ020101I1 - 2*lambdaH2Adj*lambdaT3GeneralQ020101I1 + 4*lambdaT3GeneralQ001102I1*Bar[lambdaT3] - 2*lambdaT3GeneralQ001102I2*Bar[lambdaT3] + 8*lambdaT3GeneralQ110002I1*Bar[lambdaT3] + 2*Bar[lambdaT3]*Tr[y2 . Trans[Bar[y2Conj]]] + lambdaT3GeneralQ020101I1*(Tr[y1 . Trans[Bar[y1]]] + Tr[y2 . Trans[Bar[y2]]] + Tr[y2Conj . Trans[Bar[y2Conj]]] + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaT3GeneralQ110002I1**

```text
lambdaH1*(6*lambdaT3GeneralQ001102I1 + 2*lambdaT3GeneralQ001102I2) + 4*lambdaT3*lambdaT3GeneralQ020101I1 - (33*g2^2*lambdaT3GeneralQ110002I1)/2 - (3*gY^2*lambdaT3GeneralQ110002I1)/2 + 6*lambdaH*lambdaT3GeneralQ110002I1 + 2*lambdaS2*lambdaT3GeneralQ110002I1 - 4*lambdaS2Adj*lambdaT3GeneralQ110002I1 + 40*lambdaT3GeneralQ000004I1*Bar[lambdaT3GeneralQ110002I1] + lambdaH2*(10*lambdaT3GeneralQ000013I1 + 8*lambdaT3GeneralQ110002I1 + 2*Tr[y2 . Trans[Bar[y2Conj]]]) + lambdaT3GeneralQ110002I1*(2*Tr[y2 . Trans[Bar[y2]]] + 2*Tr[y2Conj . Trans[Bar[y2Conj]]] + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]]) - 2*Tr[y2 . Trans[Bar[y2Conj]] . ye . Trans[Bar[ye]]]
```

**beta_muT3GeneralS2Sq**

```text
(6*lambdaT3GeneralQ001102I1 + 2*lambdaT3GeneralQ001102I2)*mS1Sq - 12*g2^2*muT3GeneralS2Sq + 2*lambdaS2*muT3GeneralS2Sq - 4*lambdaS2Adj*muT3GeneralS2Sq + 40*lambdaT3GeneralQ000004I1*Bar[muT3GeneralS2Sq] + mS2Sq*(10*lambdaT3GeneralQ000013I1 + 2*Tr[y2 . Trans[Bar[y2Conj]]]) + muT3GeneralS2Sq*(2*Tr[y2 . Trans[Bar[y2]]] + 2*Tr[y2Conj . Trans[Bar[y2Conj]]]) - 4*Tr[MF . Trans[Bar[y2Conj]] . y2 . Trans[Bar[MF]]]
```

### EFT1: 28 beta functions

**beta_gY**

```text
(47*gY^3)/6
```

**beta_g2**

```text
(-11*g2^3)/6
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
-12*g2^2*mS1Sq - 6*gY^2*mS1Sq + (8*lambdaS1 + 4*lambdaS1Adj)*mS1Sq + (6*lambda12 + 4*lambda12Cross)*mS2Sq + 4*muT3GeneralS2Sq*(3*Bar[lambdaT3GeneralQ001102I1] + Bar[lambdaT3GeneralQ001102I2]) + 12*lambdaT3GeneralQ001102I1*Bar[muT3GeneralS2Sq] + 4*lambdaT3GeneralQ001102I2*Bar[muT3GeneralS2Sq]
```

**beta_mS2Sq**

```text
(6*lambda12 + 4*lambda12Cross)*mS1Sq - 12*g2^2*mS2Sq + (8*lambdaS2 + 4*lambdaS2Adj)*mS2Sq + 20*muT3GeneralS2Sq*Bar[lambdaT3GeneralQ000013I1] + 20*lambdaT3GeneralQ000013I1*Bar[muT3GeneralS2Sq]
```

**beta_lambdaH**

```text
(9*g2^4)/4 + (3*gY^4)/4 + gY^2*((3*g2^2)/2 - 3*lambdaH) - 9*g2^2*lambdaH + 12*lambdaH^2 + 6*lambdaH1^2 + 3*lambdaH1Adj^2 + 6*lambdaH2^2 + 3*lambdaH2Adj^2 + 4*lambdaT3*Bar[lambdaT3] + 4*lambdaT3GeneralQ020101I1*Bar[lambdaT3GeneralQ020101I1] + 24*lambdaT3GeneralQ110002I1*Bar[lambdaT3GeneralQ110002I1] + lambdaH*(12*Tr[yd . Trans[Bar[yd]]] + 4*Tr[ye . Trans[Bar[ye]]] + 12*Tr[yu . Trans[Bar[yu]]]) - 12*Tr[yd . Trans[Bar[yd]] . yd . Trans[Bar[yd]]] - 4*Tr[ye . Trans[Bar[ye]] . ye . Trans[Bar[ye]]] - 12*Tr[yu . Trans[Bar[yu]] . yu . Trans[Bar[yu]]]
```

**beta_lambdaS1**

```text
22*g2^4 + 12*gY^4 + 6*lambda12^2 + 8*lambda12*lambda12Cross + (20*lambda12Cross^2)/3 + 4*lambdaH1^2 + gY^2*(8*g2^2 - 12*lambdaS1) - 24*g2^2*lambdaS1 + 14*lambdaS1^2 + 8*lambdaS1*lambdaS1Adj + 16*lambdaS1Adj^2 + 8*lambdaT3GeneralQ001102I1*(3*Bar[lambdaT3GeneralQ001102I1] + Bar[lambdaT3GeneralQ001102I2]) + lambdaT3GeneralQ001102I2*(8*Bar[lambdaT3GeneralQ001102I1] + (20*Bar[lambdaT3GeneralQ001102I2])/3)
```

**beta_lambdaS2**

```text
22*g2^4 + 6*lambda12^2 + 8*lambda12*lambda12Cross + (20*lambda12Cross^2)/3 + 4*lambdaH2^2 - 24*g2^2*lambdaS2 + 14*lambdaS2^2 + 8*lambdaS2*lambdaS2Adj + 16*lambdaS2Adj^2 + (832*lambdaT3GeneralQ000004I1*Bar[lambdaT3GeneralQ000004I1])/3 + (352*lambdaT3GeneralQ000013I1*Bar[lambdaT3GeneralQ000013I1])/3 + (8*lambdaT3GeneralQ001102I1*(3*Bar[lambdaT3GeneralQ001102I1] + Bar[lambdaT3GeneralQ001102I2]))/3 + lambdaT3GeneralQ001102I2*((8*Bar[lambdaT3GeneralQ001102I1])/3 + 4*Bar[lambdaT3GeneralQ001102I2]) + (16*lambdaT3GeneralQ110002I1*Bar[lambdaT3GeneralQ110002I1])/3
```

**beta_lambdaH1**

```text
6*g2^4 + 3*gY^4 - (33*g2^2*lambdaH1)/2 - (15*gY^2*lambdaH1)/2 + 6*lambdaH*lambdaH1 + 4*lambdaH1^2 + 6*lambdaH1Adj^2 + (6*lambda12 + 4*lambda12Cross)*lambdaH2 + 8*lambdaH1*lambdaS1 + 4*lambdaT3*Bar[lambdaT3] + 4*lambdaT3GeneralQ110002I1*(3*Bar[lambdaT3GeneralQ001102I1] + Bar[lambdaT3GeneralQ001102I2]) + 4*lambdaT3GeneralQ020101I1*Bar[lambdaT3GeneralQ020101I1] + 12*lambdaT3GeneralQ001102I1*Bar[lambdaT3GeneralQ110002I1] + 4*lambdaT3GeneralQ001102I2*Bar[lambdaT3GeneralQ110002I1] + lambdaH1*(4*lambdaS1Adj + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaH2**

```text
6*g2^4 + (6*lambda12 + 4*lambda12Cross)*lambdaH1 - (33*g2^2*lambdaH2)/2 - (3*gY^2*lambdaH2)/2 + 6*lambdaH*lambdaH2 + 4*lambdaH2^2 + 6*lambdaH2Adj^2 + 8*lambdaH2*lambdaS2 + 4*lambdaT3*Bar[lambdaT3] + 4*lambdaT3GeneralQ020101I1*Bar[lambdaT3GeneralQ020101I1] + 20*lambdaT3GeneralQ000013I1*Bar[lambdaT3GeneralQ110002I1] + 4*lambdaT3GeneralQ110002I1*(5*Bar[lambdaT3GeneralQ000013I1] + 4*Bar[lambdaT3GeneralQ110002I1]) + lambdaH2*(4*lambdaS2Adj + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambda12**

```text
(63*g2^4)/5 - 24*g2^2*lambda12 - 6*gY^2*lambda12 + 4*lambda12^2 + 16*lambda12Adj^2 + 4*lambda12Cross^2 + 4*lambdaH1*lambdaH2 + (8*lambda12 + 4*lambda12Cross)*lambdaS1 + 4*lambda12Cross*lambdaS1Adj + (8*lambda12 + 4*lambda12Cross)*lambdaS2 + 4*lambda12Cross*lambdaS2Adj + lambda12*(4*lambdaS1Adj + 4*lambdaS2Adj) + 2*lambdaT3*Bar[lambdaT3] + 4*lambdaT3GeneralQ001102I1*(5*Bar[lambdaT3GeneralQ000013I1] + 4*Bar[lambdaT3GeneralQ001102I1]) + 4*lambdaT3GeneralQ001102I2*(Bar[lambdaT3GeneralQ000013I1] + Bar[lambdaT3GeneralQ001102I2]) + 4*lambdaT3GeneralQ000013I1*(5*Bar[lambdaT3GeneralQ001102I1] + Bar[lambdaT3GeneralQ001102I2]) + 2*lambdaT3GeneralQ020101I1*Bar[lambdaT3GeneralQ020101I1]
```

**beta_lambdaT3**

```text
(-33*g2^2*lambdaT3)/2 - (9*gY^2*lambdaT3)/2 + 2*lambda12*lambdaT3 + 2*lambdaH*lambdaT3 + 4*lambdaH1*lambdaT3 + 4*lambdaH2*lambdaT3 + 4*lambdaT3GeneralQ001102I1*Bar[lambdaT3GeneralQ020101I1] - 2*lambdaT3GeneralQ001102I2*Bar[lambdaT3GeneralQ020101I1] + 8*lambdaT3GeneralQ110002I1*Bar[lambdaT3GeneralQ020101I1] + lambdaT3*((-4*lambda12Adj)/3 - 2*lambda12Cross - 2*lambdaH1Adj - (4*lambdaH2Adj)/3 + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaH1Adj**

```text
gY^2*(-4*g2^2 - (15*lambdaH1Adj)/2) - (33*g2^2*lambdaH1Adj)/2 + 2*lambdaH*lambdaH1Adj + 8*lambdaH1*lambdaH1Adj + 12*lambda12Adj*lambdaH2Adj + 2*lambdaH1Adj*lambdaS1 - (4*lambdaT3*Bar[lambdaT3])/3 - (4*lambdaT3GeneralQ020101I1*Bar[lambdaT3GeneralQ020101I1])/3 + lambdaH1Adj*(16*lambdaS1Adj + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaH2Adj**

```text
12*lambda12Adj*lambdaH1Adj - (33*g2^2*lambdaH2Adj)/2 - (3*gY^2*lambdaH2Adj)/2 + 2*lambdaH*lambdaH2Adj + 8*lambdaH2*lambdaH2Adj + 2*lambdaH2Adj*lambdaS2 - (8*lambdaT3*Bar[lambdaT3])/9 - (4*lambdaT3GeneralQ020101I1*Bar[lambdaT3GeneralQ020101I1])/3 + lambdaH2Adj*(16*lambdaS2Adj + 6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaS1Adj**

```text
-2*g2^4 + 12*lambda12Adj^2 - (4*lambda12Cross^2)/3 + 3*lambdaH1Adj^2 + gY^2*(8*g2^2 - 12*lambdaS1Adj) - 24*g2^2*lambdaS1Adj + 12*lambdaS1*lambdaS1Adj + 20*lambdaS1Adj^2 - (4*lambdaT3GeneralQ001102I2*Bar[lambdaT3GeneralQ001102I2])/3
```

**beta_lambdaS2Adj**

```text
-2*g2^4 + 12*lambda12Adj^2 - (4*lambda12Cross^2)/3 + 3*lambdaH2Adj^2 - 24*g2^2*lambdaS2Adj + 12*lambdaS2*lambdaS2Adj + 20*lambdaS2Adj^2 - (224*lambdaT3GeneralQ000004I1*Bar[lambdaT3GeneralQ000004I1])/3 - (44*lambdaT3GeneralQ000013I1*Bar[lambdaT3GeneralQ000013I1])/3 - (4*lambdaT3GeneralQ001102I2*Bar[lambdaT3GeneralQ001102I1])/3 - (4*lambdaT3GeneralQ001102I1*(3*Bar[lambdaT3GeneralQ001102I1] + Bar[lambdaT3GeneralQ001102I2]))/3 - (8*lambdaT3GeneralQ110002I1*Bar[lambdaT3GeneralQ110002I1])/3
```

**beta_lambda12Adj**

```text
-24*g2^2*lambda12Adj - 6*gY^2*lambda12Adj + 8*lambda12*lambda12Adj + (16*lambda12Adj*lambda12Cross)/3 + 3*lambdaH1Adj*lambdaH2Adj + 2*lambda12Adj*lambdaS1 + 16*lambda12Adj*lambdaS1Adj + 2*lambda12Adj*lambdaS2 + 16*lambda12Adj*lambdaS2Adj - (2*lambdaT3*Bar[lambdaT3])/9 - (lambdaT3GeneralQ020101I1*Bar[lambdaT3GeneralQ020101I1])/3
```

**beta_lambda12Cross**

```text
(51*g2^4)/10 - 24*g2^2*lambda12Cross - 6*gY^2*lambda12Cross + 8*lambda12*lambda12Cross + 10*lambda12Cross^2 + 2*lambda12Cross*lambdaS1 - 2*lambda12Cross*lambdaS1Adj + 2*lambda12Cross*lambdaS2 - 2*lambda12Cross*lambdaS2Adj - lambdaT3*Bar[lambdaT3] + 4*lambdaT3GeneralQ000013I1*Bar[lambdaT3GeneralQ001102I2] + 8*lambdaT3GeneralQ001102I1*Bar[lambdaT3GeneralQ001102I2] + lambdaT3GeneralQ001102I2*(4*Bar[lambdaT3GeneralQ000013I1] + 8*Bar[lambdaT3GeneralQ001102I1] + 10*Bar[lambdaT3GeneralQ001102I2]) - lambdaT3GeneralQ020101I1*Bar[lambdaT3GeneralQ020101I1]
```

**beta_lambdaT3GeneralQ000004I1**

```text
-24*g2^2*lambdaT3GeneralQ000004I1 + 12*lambdaS2*lambdaT3GeneralQ000004I1 - 8*lambdaS2Adj*lambdaT3GeneralQ000004I1 + 11*lambdaT3GeneralQ000013I1^2 + 3*lambdaT3GeneralQ001102I1^2 + 2*lambdaT3GeneralQ001102I1*lambdaT3GeneralQ001102I2 + lambdaT3GeneralQ001102I2^2 + 2*lambdaT3GeneralQ110002I1^2
```

**beta_lambdaT3GeneralQ000013I1**

```text
-24*g2^2*lambdaT3GeneralQ000013I1 + 22*lambdaS2*lambdaT3GeneralQ000013I1 + lambda12*(6*lambdaT3GeneralQ001102I1 + 2*lambdaT3GeneralQ001102I2) + lambda12Cross*(4*lambdaT3GeneralQ001102I1 + 4*lambdaT3GeneralQ001102I2) + 4*lambdaH2*lambdaT3GeneralQ110002I1 + 88*lambdaT3GeneralQ000004I1*Bar[lambdaT3GeneralQ000013I1]
```

**beta_lambdaT3GeneralQ001102I1**

```text
-24*g2^2*lambdaT3GeneralQ001102I1 - 6*gY^2*lambdaT3GeneralQ001102I1 + 2*lambdaS2*lambdaT3GeneralQ001102I1 + lambda12*(10*lambdaT3GeneralQ000013I1 + 8*lambdaT3GeneralQ001102I1) + lambdaS2Adj*(-4*lambdaT3GeneralQ001102I1 - 2*lambdaT3GeneralQ001102I2) + lambdaS1Adj*(4*lambdaT3GeneralQ001102I1 + 2*lambdaT3GeneralQ001102I2) + lambdaS1*(8*lambdaT3GeneralQ001102I1 + 2*lambdaT3GeneralQ001102I2) + lambda12Cross*(4*lambdaT3GeneralQ000013I1 + 4*lambdaT3GeneralQ001102I2) + 2*lambdaT3*lambdaT3GeneralQ020101I1 + 4*lambdaH1*lambdaT3GeneralQ110002I1 + 8*lambdaT3GeneralQ000004I1*(5*Bar[lambdaT3GeneralQ001102I1] + Bar[lambdaT3GeneralQ001102I2])
```

**beta_lambdaT3GeneralQ001102I2**

```text
-24*g2^2*lambdaT3GeneralQ001102I2 - 6*gY^2*lambdaT3GeneralQ001102I2 + 8*lambda12*lambdaT3GeneralQ001102I2 + 2*lambdaS1*lambdaT3GeneralQ001102I2 - 2*lambdaS1Adj*lambdaT3GeneralQ001102I2 + 2*lambdaS2*lambdaT3GeneralQ001102I2 + 2*lambdaS2Adj*lambdaT3GeneralQ001102I2 + lambda12Cross*(8*lambdaT3GeneralQ000013I1 + 16*lambdaT3GeneralQ001102I1 + 20*lambdaT3GeneralQ001102I2) - 2*lambdaT3*lambdaT3GeneralQ020101I1 + 16*lambdaT3GeneralQ000004I1*Bar[lambdaT3GeneralQ001102I2]
```

**beta_lambdaT3GeneralQ020101I1**

```text
(-33*g2^2*lambdaT3GeneralQ020101I1)/2 - (9*gY^2*lambdaT3GeneralQ020101I1)/2 + 2*lambda12*lambdaT3GeneralQ020101I1 - 2*lambda12Adj*lambdaT3GeneralQ020101I1 - 2*lambda12Cross*lambdaT3GeneralQ020101I1 + 2*lambdaH*lambdaT3GeneralQ020101I1 + 4*lambdaH1*lambdaT3GeneralQ020101I1 - 2*lambdaH1Adj*lambdaT3GeneralQ020101I1 + 4*lambdaH2*lambdaT3GeneralQ020101I1 - 2*lambdaH2Adj*lambdaT3GeneralQ020101I1 + 4*lambdaT3GeneralQ001102I1*Bar[lambdaT3] - 2*lambdaT3GeneralQ001102I2*Bar[lambdaT3] + 8*lambdaT3GeneralQ110002I1*Bar[lambdaT3] + lambdaT3GeneralQ020101I1*(6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_lambdaT3GeneralQ110002I1**

```text
lambdaH1*(6*lambdaT3GeneralQ001102I1 + 2*lambdaT3GeneralQ001102I2) + 4*lambdaT3*lambdaT3GeneralQ020101I1 - (33*g2^2*lambdaT3GeneralQ110002I1)/2 - (3*gY^2*lambdaT3GeneralQ110002I1)/2 + 6*lambdaH*lambdaT3GeneralQ110002I1 + 2*lambdaS2*lambdaT3GeneralQ110002I1 - 4*lambdaS2Adj*lambdaT3GeneralQ110002I1 + lambdaH2*(10*lambdaT3GeneralQ000013I1 + 8*lambdaT3GeneralQ110002I1) + 40*lambdaT3GeneralQ000004I1*Bar[lambdaT3GeneralQ110002I1] + lambdaT3GeneralQ110002I1*(6*Tr[yd . Trans[Bar[yd]]] + 2*Tr[ye . Trans[Bar[ye]]] + 6*Tr[yu . Trans[Bar[yu]]])
```

**beta_muT3GeneralS2Sq**

```text
(6*lambdaT3GeneralQ001102I1 + 2*lambdaT3GeneralQ001102I2)*mS1Sq + 10*lambdaT3GeneralQ000013I1*mS2Sq - 12*g2^2*muT3GeneralS2Sq + 2*lambdaS2*muT3GeneralS2Sq - 4*lambdaS2Adj*muT3GeneralS2Sq + 40*lambdaT3GeneralQ000004I1*Bar[muT3GeneralS2Sq]
```

## Interpretation and limitations

- UV RGEs include symmetry-allowed scalar quartics, bilinears, and the additional Yukawa where gauge invariance permits it.
- EFT1 RGEs describe the renormalisable SM+S1+S2 theory; the heavy fermion and its Yukawas have been removed.
- EFT1 higher-dimensional Wilson RGEs and final C5 have **not** been rederived with the full general interactions.
- Mixed scalar mass terms require mass-basis diagonalisation for B(-1) and C(-1) before general matching.
- The existing heavy-fermion flavour-mass parameterisation is not certified as the most general non-diagonal flavour mass basis.

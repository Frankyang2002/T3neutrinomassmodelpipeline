(* General D(-2) direct EFT1 LLSS -> LLHH contribution.
   Verified against all four Matchete PL seed terms and one-generation beta.
   Coefficient is beta_16pi2 * Log[muLow/muHigh] -- NO hbar factor here.
   The lower scalar matching must carry this pure-SM contribution SEPARATELY,
   not include it as a heavy operator in the new one-loop scalar Match.
   This file alone does NOT execute scalar matching.
*)
ClearAll[EFT1D2C12, EFT1D2C1Conj, EFT1D2WilsonRunningHeavyInsertion,
  EFT1D2WilsonRunningWeinbergCoefficient, EFT1D2WilsonRunningMetadata];
If[!KeyExistsQ[GetCouplings[], MF],
  DefineCoupling[MF, SelfConjugate -> True];
];
If[!KeyExistsQ[GetCouplings[], MuScalar],
  DefineCoupling[MuScalar, SelfConjugate -> True];
];
EFT1D2C12[p_, q_] := Module[{r = Unique["nfl$"]},
  Sum[(Bar[Coupling[y1, {p, Index[r, NFlavor]}, 0]]*Bar[Coupling[y2, {q, Index[r, NFlavor]}, 0]] + Bar[Coupling[y2, {p, Index[r, NFlavor]}, 0]]*Bar[Coupling[y1, {q, Index[r, NFlavor]}, 0]])/2, {r, 1, NF}]
];
EFT1D2C1Conj[p_, q_] := Module[{r = Unique["nfl$"]},
  Sum[(Bar[Coupling[y1, {p, Index[r, NFlavor]}, 0]]*Bar[Coupling[y2Conj, {q, Index[r, NFlavor]}, 0]] + Bar[Coupling[y2Conj, {p, Index[r, NFlavor]}, 0]]*Bar[Coupling[y1, {q, Index[r, NFlavor]}, 0]])/2, {r, 1, NF}]
];
EFT1D2WilsonRunningHeavyInsertion = 0;
EFT1D2WilsonRunningWeinbergCoefficient[p_, q_] :=
  Log[Coupling[MuScalar, {}, 0]/Coupling[MF, {}, 0]] *
  (((4/3)*3^(1/2)*Coupling[lambdaT3, {}, 0])*EFT1D2C12[p, q] + ((4/3)*3^(1/2)*Bar[Coupling[lambdaT3ExtraQ020101I1, {}, 0]])*EFT1D2C1Conj[p, q]) /
  Coupling[MF, {}, 0];
EFT1D2WilsonRunningMetadata = <|
  "Status" -> "PreparedNotYetMatched",
  "Model" -> "D_-2",
  "DirectWeinbergOnlyAtOneLoop" -> True,
  "FullLeptonFlavor" -> True,
  "HeavyOperatorInsertion" -> 0,
  "HbarFactorIncluded" -> False,
  "PhysicalMassEigenvalueRelationsImposed" -> False,
  "ScalarHardThresholdPerformed" -> False
|>;

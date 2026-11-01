(* Auto-generated direct one-loop Weinberg transport.

   The heavy LLSS self-running is O(hbar).  Feeding that correction through
   the scalar loop would be O(hbar^2), so the authoritative fixed-one-loop
   path sets the heavy running insertion to zero and carries only the direct
   LLSS -> Weinberg mixing.
*)

ClearAll[
  EFT1C12,
  EFT1DeltaWeinberg,
  EFT1WilsonRunningHeavyInsertion,
  EFT1WilsonRunningWeinbergCoefficient,
  EFT1WilsonRunningLogReplacement,
  EFT1WilsonFlavorRunningMetadata
];

If[!KeyExistsQ[GetCouplings[], EFT1RunLog],
  DefineCoupling[EFT1RunLog, SelfConjugate -> True];
];

EFT1WilsonRunningLogReplacement =
  Coupling[EFT1RunLog, {}, 0] ->
    Log[MS/Coupling[MF, {}, 0]];

EFT1C12[p_, q_] := Module[{r = Unique["nfl$"]},
  (
    Bar[Coupling[y1, {p, Index[r, NFlavor]}, 0]] *
    Bar[Coupling[y2, {q, Index[r, NFlavor]}, 0]] +
    Bar[Coupling[y2, {p, Index[r, NFlavor]}, 0]] *
    Bar[Coupling[y1, {q, Index[r, NFlavor]}, 0]]
  ) / (2*Coupling[MF, {}, 0])
];

EFT1DeltaWeinberg[p_, q_] :=
  ((4/3)*3^(1/2)*Coupling[lambdaT3, {}, 0]*Coupling[EFT1RunLog, {}, 0]) * EFT1C12[p, q];

EFT1WilsonRunningHeavyInsertion = 0;

EFT1WilsonRunningWeinbergCoefficient[p_, q_] :=
  EFT1DeltaWeinberg[p, q];

EFT1WilsonFlavorRunningMetadata = <|
  "MuHigh" -> "MF",
  "MuLow" -> "MS",
  "LogRatio" -> "log(MS/MF)",
  "HeavyOperatorCount" -> 0,
  "WeinbergOperatorCount" -> 1,
  "DirectWeinbergOnlyAtOneLoop" -> True,
  "MixedC12FlavorSymmetrized" -> True,
  "OneGenerationRegression" -> True,
  "EqualScaleRunningVanishes" -> True
|>;

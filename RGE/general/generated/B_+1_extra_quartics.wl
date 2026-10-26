(* AUTO-GENERATED RGBeta extension for B_+1.
   Invariant coverage checked against the independent Matchete count.
   Do not identify these coupling values with Matchete couplings:
   tensor basis orientation and normalisation remain uncertified.
   Do not activate native RGBeta gauge option: v1.2.0 on
   Windows rejects it at AddQuartic option validation.
   Gauge-structure checks are in the independent Python audit. *)
T3GeneralQuarticAddedSymbols = {};
T3GeneralQuarticExpectedCount = 1;
T3GeneralQuarticReusedCount = 0;

AddQuartic[lambdaT3GeneralQ022000I1,
  {Bar @ H, Bar @ H, S1, S1},
  GroupInvariant -> (((1/2) * (del[SU2L @ fund, #1, #3] * del[SU2L @ fund, #2, #4]) + (1/2) * (del[SU2L @ fund, #1, #4] * del[SU2L @ fund, #2, #3])) &),
  SelfConjugate -> False
];
AppendTo[T3GeneralQuarticAddedSymbols, lambdaT3GeneralQ022000I1];

If[Length[T3GeneralQuarticAddedSymbols] + T3GeneralQuarticReusedCount =!= T3GeneralQuarticExpectedCount,
  Print["ERROR: General RGBeta quartic coverage mismatch"]; Abort[]
];

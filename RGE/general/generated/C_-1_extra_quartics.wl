(* AUTO-GENERATED RGBeta extension for C_-1.
   Invariant coverage checked against the independent Matchete count.
   Do not identify these coupling values with Matchete couplings:
   tensor basis orientation and normalisation remain uncertified.
   Do not activate native RGBeta gauge option: v1.2.0 on
   Windows rejects it at AddQuartic option validation.
   Gauge-structure checks are in the independent Python audit. *)
T3GeneralQuarticAddedSymbols = {};
T3GeneralQuarticExpectedCount = 7;
T3GeneralQuarticReusedCount = 6;

AddQuartic[lambdaT3GeneralQ110101I1,
  {H, Bar @ H, Bar @ S1, Bar @ S2},
  GroupInvariant -> ((del[SU2L @ fund, #4, #1] * eps[SU2L @ fund, #2, #3]) &),
  SelfConjugate -> False
];
AppendTo[T3GeneralQuarticAddedSymbols, lambdaT3GeneralQ110101I1];

If[Length[T3GeneralQuarticAddedSymbols] + T3GeneralQuarticReusedCount =!= T3GeneralQuarticExpectedCount,
  Print["ERROR: General RGBeta quartic coverage mismatch"]; Abort[]
];

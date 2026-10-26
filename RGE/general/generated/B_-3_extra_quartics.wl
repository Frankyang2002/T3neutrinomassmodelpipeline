(* AUTO-GENERATED RGBeta extension for B_-3.
   Invariant coverage checked against the independent Matchete count.
   Do not identify these coupling values with Matchete couplings:
   tensor basis orientation and normalisation remain uncertified.
   Do not activate native RGBeta gauge option: v1.2.0 on
   Windows rejects it at AddQuartic option validation.
   Gauge-structure checks are in the independent Python audit. *)
T3GeneralQuarticAddedSymbols = {};
T3GeneralQuarticExpectedCount = 1;
T3GeneralQuarticReusedCount = 0;

AddQuartic[lambdaT3GeneralQ020002I1,
  {Bar @ H, Bar @ H, Bar @ S2, Bar @ S2},
  GroupInvariant -> (((1/2) * (eps[SU2L @ fund, #1, #3] * eps[SU2L @ fund, #2, #4]) + (1/2) * (eps[SU2L @ fund, #1, #4] * eps[SU2L @ fund, #2, #3])) &),
  SelfConjugate -> False
];
AppendTo[T3GeneralQuarticAddedSymbols, lambdaT3GeneralQ020002I1];

If[Length[T3GeneralQuarticAddedSymbols] + T3GeneralQuarticReusedCount =!= T3GeneralQuarticExpectedCount,
  Print["ERROR: General RGBeta quartic coverage mismatch"]; Abort[]
];

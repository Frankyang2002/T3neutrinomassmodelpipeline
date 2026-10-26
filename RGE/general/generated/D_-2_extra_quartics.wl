(* AUTO-GENERATED RGBeta extension for D_-2.
   Invariant coverage checked against the independent Matchete count.
   Do not identify these coupling values with Matchete couplings:
   tensor basis orientation and normalisation remain uncertified.
   Do not activate native RGBeta gauge option: v1.2.0 on
   Windows rejects it at AddQuartic option validation.
   Gauge-structure checks are in the independent Python audit. *)
T3GeneralQuarticAddedSymbols = {};
T3GeneralQuarticExpectedCount = 5;
T3GeneralQuarticReusedCount = 0;

AddQuartic[lambdaT3GeneralQ000004I1,
  {Bar @ S2, Bar @ S2, Bar @ S2, Bar @ S2},
  GroupInvariant -> (1 &),
  SelfConjugate -> False
];
AppendTo[T3GeneralQuarticAddedSymbols, lambdaT3GeneralQ000004I1];

AddQuartic[lambdaT3GeneralQ000013I1,
  {S2, Bar @ S2, Bar @ S2, Bar @ S2},
  GroupInvariant -> (1 &),
  SelfConjugate -> False
];
AppendTo[T3GeneralQuarticAddedSymbols, lambdaT3GeneralQ000013I1];

AddQuartic[lambdaT3GeneralQ001102I1,
  {S1, Bar @ S1, Bar @ S2, Bar @ S2},
  GroupInvariant -> (del[SU2L @ S2, #1, #2] &),
  SelfConjugate -> False
];
AppendTo[T3GeneralQuarticAddedSymbols, lambdaT3GeneralQ001102I1];

AddQuartic[lambdaT3GeneralQ020101I1,
  {Bar @ H, Bar @ H, Bar @ S1, Bar @ S2},
  GroupInvariant -> (delS2[SU2L, #3, #1, #2] &),
  SelfConjugate -> False
];
AppendTo[T3GeneralQuarticAddedSymbols, lambdaT3GeneralQ020101I1];

AddQuartic[lambdaT3GeneralQ110002I1,
  {H, Bar @ H, Bar @ S2, Bar @ S2},
  GroupInvariant -> (del[SU2L @ fund, #2, #1] &),
  SelfConjugate -> False
];
AppendTo[T3GeneralQuarticAddedSymbols, lambdaT3GeneralQ110002I1];

If[Length[T3GeneralQuarticAddedSymbols] + T3GeneralQuarticReusedCount =!= T3GeneralQuarticExpectedCount,
  Print["ERROR: General RGBeta quartic coverage mismatch"]; Abort[]
];

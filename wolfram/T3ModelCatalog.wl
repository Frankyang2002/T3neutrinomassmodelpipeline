(* T3ModelCatalog.wl
   Representation-level catalogue for the T3 one-loop topology.

   Conventions:
     Q = T3 + Y
     Y(S1) = alpha/2
     Y(S2) = (alpha+2)/2
     Y(F)  = (alpha+1)/2

   The historical A-E labels are retained as aliases, but the real model
   constructor is T3ModelFromDimensions[d1,d2,dF,alpha].
*)
ClearAll[
  T3YukawaAllowedQ,
  T3ScalarMixingAllowedQ,
  T3DimensionsAllowedQ,
  T3ModelFromDimensions,
  T3ModelFromClass
];

(* Because L is an SU(2) doublet,
       2 x dF = (dF-1) + (dF+1),
   with the dF-1 term absent for dF=1. *)
T3YukawaAllowedQ[dS_Integer?Positive, dF_Integer?Positive] :=
  Abs[dS - dF] === 1;

(* The two identical Higgs fields are bosons, so their SU(2) pair is in the
   symmetric triplet channel.  Therefore S1 x S2^dagger must contain a
   triplet.  SU(2) tensor products are multiplicity-free. *)
T3ScalarMixingAllowedQ[d1_Integer?Positive, d2_Integer?Positive] := Module[
  {j1, j2},
  j1 = (d1 - 1)/2;
  j2 = (d2 - 1)/2;
  TrueQ[
    Abs[j1 - j2] <= 1 <= j1 + j2 &&
    IntegerQ[j1 + j2]
  ]
];

T3DimensionsAllowedQ[d1_Integer?Positive, d2_Integer?Positive, dF_Integer?Positive] :=
  TrueQ[
    T3YukawaAllowedQ[d1, dF] &&
    T3YukawaAllowedQ[d2, dF] &&
    T3ScalarMixingAllowedQ[d1, d2]
  ];

(* Matchete flavor indices require dimension > 1.  Three heavy-fermion
   generations are used by default, matching the previous implementation. *)
T3ModelFromDimensions[
  d1_Integer?Positive,
  d2_Integer?Positive,
  dF_Integer?Positive,
  alpha_Integer,
  multiplicity_: 3
] := Module[{y1, y2, yF},
  If[!T3DimensionsAllowedQ[d1, d2, dF], Return[$Failed]];
  If[!IntegerQ[multiplicity] || multiplicity < 2, Return[$Failed]];

  y1 = alpha/2;
  y2 = (alpha + 2)/2;
  yF = (alpha + 1)/2;

  <|
    "Class" -> "T3-d" <> ToString[d1] <> "-d" <> ToString[d2] <> "-F" <> ToString[dF],
    "Alpha" -> alpha,
    "Lepton" -> <|"Type" -> "Fermion", "SU2" -> 2, "Y" -> -1/2|>,
    "Higgs" -> <|"Type" -> "ComplexScalar", "SU2" -> 2, "Y" -> 1/2|>,
    "Scalar1" -> <|"Type" -> "ComplexScalar", "SU2" -> d1, "Y" -> y1, "MassSymbol" -> MS1|>,
    "Scalar2" -> <|"Type" -> "ComplexScalar", "SU2" -> d2, "Y" -> y2, "MassSymbol" -> MS2|>,
    "Fermion" -> <|"Type" -> "Fermion", "SU2" -> dF, "Y" -> yF, "Multiplicity" -> multiplicity, "MassSymbol" -> MF|>
  |>
];

(* Backwards-compatible aliases for the five singlet/doublet/triplet classes. *)
T3ModelFromClass[class_String, alpha_Integer, multiplicity_: 3] := Module[{dims},
  dims = Switch[
    ToUpperCase[class],
    "A", {1, 3, 2},
    "B", {2, 2, 1},
    "C", {2, 2, 3},
    "D", {3, 1, 2},
    "E", {3, 3, 2},
    _, Return[$Failed]
  ];
  Association[
    T3ModelFromDimensions[dims[[1]], dims[[2]], dims[[3]], alpha, multiplicity],
    "Class" -> "T3-" <> ToUpperCase[class]
  ]
];

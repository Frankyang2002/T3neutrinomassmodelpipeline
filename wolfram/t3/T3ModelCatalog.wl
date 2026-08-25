(* T3ModelCatalog.wl
   Representation-level catalogue for the T3 one-loop topology.

   Hypercharge convention:
     Q = T3 + Y,
     Y(S1) = alpha/2,
     Y(S2) = (alpha + 2)/2,
     Y(F)  = (alpha + 1)/2.

   The historical A-E labels are retained as aliases.  Generic models are
   constructed from the SU(2) dimensions (dS1, dS2, dF) and alpha.
*)
ClearAll[
  T3YukawaAllowedQ,
  T3ScalarMixingAllowedQ,
  T3DimensionsAllowedQ,
  T3ModelFromDimensions,
  T3ModelFromClass,
  T3Hypercharges,
  T3ClassDimensions
];

(* Historical T3 classes from the singlet/doublet/triplet classification. *)
T3ClassDimensions = <|
  "A" -> {1, 3, 2},
  "B" -> {2, 2, 1},
  "C" -> {2, 2, 3},
  "D" -> {3, 1, 2},
  "E" -> {3, 3, 2}
|>;

(* Hypercharges fixed by the T3 interaction structure. *)
T3Hypercharges[alpha_Integer] := <|
  "Scalar1" -> alpha/2,
  "Scalar2" -> (alpha + 2)/2,
  "Fermion" -> (alpha + 1)/2
|>;

(* Since L is a doublet, a Yukawa vertex L S F can only couple scalar
   representations with dS = dF +/- 1. *)
T3YukawaAllowedQ[dS_Integer?Positive, dF_Integer?Positive] :=
  Abs[dS - dF] === 1;

(* The identical Higgs fields occupy the symmetric SU(2) triplet channel.
   Hence S1 x S2^dagger must contain j=1. *)
T3ScalarMixingAllowedQ[d1_Integer?Positive, d2_Integer?Positive] := Module[
  {j1 = (d1 - 1)/2, j2 = (d2 - 1)/2},
  TrueQ[
    Abs[j1 - j2] <= 1 <= j1 + j2 &&
    IntegerQ[j1 + j2]
  ]
];

T3DimensionsAllowedQ[
  d1_Integer?Positive,
  d2_Integer?Positive,
  dF_Integer?Positive
] := TrueQ[
  T3YukawaAllowedQ[d1, dF] &&
  T3YukawaAllowedQ[d2, dF] &&
  T3ScalarMixingAllowedQ[d1, d2]
];

(* Matchete flavour indices require multiplicity > 1.  The default of three
   heavy-fermion generations preserves the existing project convention. *)
T3ModelFromDimensions[
  d1_Integer?Positive,
  d2_Integer?Positive,
  dF_Integer?Positive,
  alpha_Integer,
  multiplicity_: 3
] := Module[{hypercharges},
  If[!T3DimensionsAllowedQ[d1, d2, dF], Return[$Failed]];
  If[!IntegerQ[multiplicity] || multiplicity < 2, Return[$Failed]];

  hypercharges = T3Hypercharges[alpha];

  <|
    "Class" -> "T3-d" <> ToString[d1] <> "-d" <> ToString[d2] <> "-F" <> ToString[dF],
    "Alpha" -> alpha,
    "Lepton" -> <|
      "Type" -> "Fermion",
      "SU2" -> 2,
      "Y" -> -1/2
    |>,
    "Higgs" -> <|
      "Type" -> "ComplexScalar",
      "SU2" -> 2,
      "Y" -> 1/2
    |>,
    "Scalar1" -> <|
      "Type" -> "ComplexScalar",
      "SU2" -> d1,
      "Y" -> hypercharges["Scalar1"],
      "MassSymbol" -> MS1
    |>,
    "Scalar2" -> <|
      "Type" -> "ComplexScalar",
      "SU2" -> d2,
      "Y" -> hypercharges["Scalar2"],
      "MassSymbol" -> MS2
    |>,
    "Fermion" -> <|
      "Type" -> "Fermion",
      "SU2" -> dF,
      "Y" -> hypercharges["Fermion"],
      "Multiplicity" -> multiplicity,
      "MassSymbol" -> MF
    |>
  |>
];

(* Backwards-compatible aliases for the five historical classes. *)
T3ModelFromClass[class_String, alpha_Integer, multiplicity_: 3] := Module[
  {label = ToUpperCase[class], dims, model},

  dims = Lookup[T3ClassDimensions, label, Missing["UnknownClass"]];
  If[MissingQ[dims], Return[$Failed]];

  model = T3ModelFromDimensions[
    dims[[1]],
    dims[[2]],
    dims[[3]],
    alpha,
    multiplicity
  ];
  If[model === $Failed, Return[$Failed]];

  Association[model, "Class" -> "T3-" <> label]
];

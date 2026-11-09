(* Scalar-potential invariant basis and candidate construction for T3 models.
   Extracted from LagrangianBuilder.wl; physics behaviour is unchanged. *)

(* All non-singlet irreps have one Matchete representation index regardless of
   their dimension. *)
ScalarNorm[which_Integer, d_Integer?Positive, i_] := If[
  d === 1,
  ScalarFieldValue[which, True] ScalarFieldValue[which, False],
  ScalarFieldValue[which, True, i] ScalarFieldValue[which, False, i]
];

ScalarSelf[which_, d_, coupling_] := Module[{i, j},
  coupling[]/2 ScalarNorm[which, d, i] ScalarNorm[which, d, j]
];

HiggsPortal[which_, d_, coupling_] := Module[{i, j},
  coupling[] Bar[H[i]] H[i] ScalarNorm[which, d, j]
];

CrossScalarPortal[d1_, d2_] := Module[{i, j},
  lambda12[] ScalarNorm[1, d1, i] ScalarNorm[2, d2, j]
];

(* Convert an invariant tensor to a flat vector so linear independence can be
   tested even when the SU(2) invariant is just the scalar 1. *)
TensorVector[tensor_] := If[
  ArrayDepth[tensor] === 0,
  {tensor},
  Flatten[Normal[tensor]]
];

(* Exchange two representation-index slots of an invariant tensor. *)
SwapTensorSlots[tensor_, first_Integer, second_Integer] := Module[{permutation},
  permutation = Range[ArrayDepth[tensor]];
  permutation[[{first, second}]] = permutation[[{second, first}]];
  Transpose[tensor, permutation]
];

(* Project an invariant tensor onto all requested identical-boson symmetric
   subspaces.  Each pair refers to positions after singlet slots are removed. *)
SymmetrizeTensorPairs[tensor_, pairs_List] := Fold[
  (#1 + SwapTensorSlots[#1, #2[[1]], #2[[2]]])/2 &,
  tensor,
  pairs
];

(* Keep a linearly independent basis from a list of invariant tensors. *)
IndependentTensorBasis[tensors_List] := Module[
  {selected = {}, vectors = {}, vector, newRank},

  Scan[
    Function[tensor,
      vector = TensorVector[tensor];
      newRank = MatrixRank[Append[vectors, vector]];
      If[newRank > Length[vectors],
        AppendTo[selected, tensor];
        AppendTo[vectors, vector];
      ];
    ],
    tensors
  ];

  selected
];

(* Return every physical SU(2) invariant tensor for a field structure.
   symmetricPairs is a list such as {{1,3},{2,4}} for
   S^\dagger S S^\dagger S. *)
PhysicalSU2InvariantBasis[
  dims_List,
  objectConjugated_List,
  symmetricPairs_: {},
  smPositions_: {}
] := Module[
  {keep, keptDims, keptConj, algebraReps, cgReps, tensors,
   positionMap, pairsAfterDrop, projected, basis},

  If[Length[dims] =!= Length[objectConjugated], Return[$Failed]];
  If[!AllTrue[dims, IntegerQ[#] && Positive[#] &], Return[$Failed]];

  Scan[EnsureBSMRepresentation, DeleteDuplicates[dims]];

  keep = Flatten@Position[dims, _?(# > 1 &)];
  keptDims = dims[[keep]];
  keptConj = objectConjugated[[keep]];

  If[keptDims === {},
    Return[
      <|
        "Basis" -> {1},
        "CGReps" -> {},
        "Keep" -> {},
        "RawCount" -> 1,
        "PhysicalCount" -> 1
      |>
    ]
  ];

  algebraReps = MapThread[
    Function[{d, conjugated},
      Module[{label = SU2DynkinLabel[d]},
        If[EvenQ[d] && TrueQ[conjugated], CRep[label], label]
      ]
    ],
    {keptDims, keptConj}
  ];

  cgReps = MapThread[
    Function[{originalPosition, d, conjugated},
      Module[{baseRep},
        baseRep = If[
          MemberQ[smPositions, originalPosition],
          SU2L[fund],
          BSMRepresentationName[d]
        ];

        If[EvenQ[d] && !TrueQ[conjugated], Bar[baseRep], baseRep]
      ]
    ],
    {keep, keptDims, keptConj}
  ];

  tensors = Check[InvariantTensors[SU[2], algebraReps], $Failed];

  If[tensors === $Failed || !ListQ[tensors] || tensors === {},
    Return[$Failed]
  ];

  positionMap = AssociationThread[keep -> Range[Length[keep]]];

  pairsAfterDrop = Cases[
    symmetricPairs,
    {first_Integer, second_Integer} /;
      KeyExistsQ[positionMap, first] &&
      KeyExistsQ[positionMap, second] :>
        {positionMap[first], positionMap[second]}
  ];

  projected = SymmetrizeTensorPairs[#, pairsAfterDrop] & /@ tensors;
  basis = IndependentTensorBasis[projected];

  If[basis === {}, Return[$Failed]];

  <|
    "Basis" -> basis,
    "CGReps" -> cgReps,
    "Keep" -> keep,
    "RawCount" -> Length[tensors],
    "PhysicalCount" -> Length[basis]
  |>
];

(* Define a dynamically named real Matchete coupling. *)
DefineDynamicRealCoupling[coupling_Symbol] := Module[{defined},
  defined = Check[
    Apply[DefineCoupling, {coupling, SelfConjugate -> True}];
    coupling,
    $Failed
  ];
  defined
];

(* Define one CG tensor and one independent coupling for every physical
   invariant in a scalar-potential field structure. *)
DefineSU2InvariantFamily[
  cgBase_String,
  couplingBase_String,
  dims_List,
  objectConjugated_List,
  symmetricPairs_: {},
  smPositions_: {}
] := Module[
  {data, basis, cgReps, keep, cgNames, couplings, cgName, couplingName,
   cgDefined, couplingDefined},

  data = PhysicalSU2InvariantBasis[
    dims,
    objectConjugated,
    symmetricPairs,
    smPositions
  ];

  If[data === $Failed, Return[$Failed]];

  basis = data["Basis"];
  cgReps = data["CGReps"];
  keep = data["Keep"];

  cgNames = Table[
    If[cgReps === {},
      None,
      cgName = Symbol[cgBase <> "Inv" <> ToString[index]];
      cgDefined = Check[
        Apply[DefineCG, {cgName, cgReps, basis[[index]]}];
        cgName,
        $Failed
      ];
      If[cgDefined === $Failed, Return[$Failed]];
      cgName
    ],
    {index, Length[basis]}
  ];

  couplings = Table[
    couplingName = Symbol[couplingBase <> "Inv" <> ToString[index]];
    couplingDefined = DefineDynamicRealCoupling[couplingName];
    If[couplingDefined === $Failed, Return[$Failed]];
    couplingName,
    {index, Length[basis]}
  ];

  Print[
    "  ", couplingBase,
    ": raw SU(2) invariants=", data["RawCount"],
    ", physical invariants=", data["PhysicalCount"]
  ];

  <|
    "CGs" -> cgNames,
    "Couplings" -> couplings,
    "Keep" -> keep,
    "Count" -> Length[basis]
  |>
];

(* Return a scalar field with or without its SU(2) index. *)
ScalarPotentialFieldValue[
  which_Integer,
  conjugated_,
  d_Integer?Positive,
  index_
] := If[
  d === 1,
  ScalarFieldValue[which, conjugated],
  ScalarFieldValue[which, conjugated, index]
];

(* Build all independent (S^\dagger S)(S^\dagger S) contractions. *)
BuildScalarSelfInvariantCandidates[
  which_Integer,
  d_Integer?Positive] := Module[
  {family, indices, labels, cg, coupling, fieldProduct, baseName},

  baseName = "lambdaS" <> ToString[which];

  family = DefineSU2InvariantFamily[
    "T3Scalar" <> ToString[which] <> "SelfCG",
    baseName,
    {d, d, d, d},
    {True, False, True, False},
    {{1, 3}, {2, 4}},
    {}
  ];

  If[family === $Failed, Return[{}]];

  Table[
    indices = Array[Unique["s"] &, 4];
    labels = indices[[family["Keep"]]];
    cg = family["CGs"][[invariant]];
    coupling = family["Couplings"][[invariant]];

    fieldProduct =
      ScalarPotentialFieldValue[which, True, d, indices[[1]]] *
      ScalarPotentialFieldValue[which, False, d, indices[[2]]] *
      ScalarPotentialFieldValue[which, True, d, indices[[3]]] *
      ScalarPotentialFieldValue[which, False, d, indices[[4]]];

    <|
      "Name" -> "Scalar" <> ToString[which] <> "SelfInv" <> ToString[invariant],
      "Class" -> "ScalarSelf",
      "Expression" ->
        coupling[]/2 fieldProduct If[cg === None, 1, CG[cg, labels]]
    |>,
    {invariant, family["Count"]}
  ]
];

(* Build all independent (H^\dagger H)(S^\dagger S) contractions. *)
BuildHiggsPortalInvariantCandidates[
  which_Integer,
  d_Integer?Positive] := Module[
  {family, indices, labels, cg, coupling, fieldProduct, baseName},

  baseName = "lambdaH" <> ToString[which];

  family = DefineSU2InvariantFamily[
    "T3HiggsPortal" <> ToString[which] <> "CG",
    baseName,
    {2, 2, d, d},
    {True, False, True, False},
    {},
    {1, 2}
  ];

  If[family === $Failed, Return[{}]];

  Table[
    indices = Array[Unique["p"] &, 4];
    labels = indices[[family["Keep"]]];
    cg = family["CGs"][[invariant]];
    coupling = family["Couplings"][[invariant]];

    fieldProduct =
      Bar[H[indices[[1]]]] H[indices[[2]]] *
      ScalarPotentialFieldValue[which, True, d, indices[[3]]] *
      ScalarPotentialFieldValue[which, False, d, indices[[4]]];

    <|
      "Name" -> "HiggsPortal" <> ToString[which] <> "Inv" <> ToString[invariant],
      "Class" -> "Portal",
      "Expression" ->
        1/2 PlusHc[
          coupling[] fieldProduct If[cg === None, 1, CG[cg, labels]]
        ]
    |>,
    {invariant, family["Count"]}
  ]
];

(* Build all independent (S1^\dagger S1)(S2^\dagger S2) contractions. *)
BuildCrossScalarInvariantCandidates[
  d1_Integer?Positive,
  d2_Integer?Positive] := Module[
  {family, indices, labels, cg, coupling, fieldProduct},

  family = DefineSU2InvariantFamily[
    "T3CrossScalarCG",
    "lambda12",
    {d1, d1, d2, d2},
    {True, False, True, False},
    {},
    {}
  ];

  If[family === $Failed, Return[{}]];

  Table[
    indices = Array[Unique["c"] &, 4];
    labels = indices[[family["Keep"]]];
    cg = family["CGs"][[invariant]];
    coupling = family["Couplings"][[invariant]];

    fieldProduct =
      ScalarPotentialFieldValue[1, True, d1, indices[[1]]] *
      ScalarPotentialFieldValue[1, False, d1, indices[[2]]] *
      ScalarPotentialFieldValue[2, True, d2, indices[[3]]] *
      ScalarPotentialFieldValue[2, False, d2, indices[[4]]];

    <|
      "Name" -> "ScalarCrossPortalInv" <> ToString[invariant],
      "Class" -> "Portal",
      "Expression" ->
        coupling[] fieldProduct If[cg === None, 1, CG[cg, labels]]
    |>,
    {invariant, family["Count"]}
  ]
];

(* Build the complete scalar-potential invariant basis for the two T3 scalars. *)
BuildGeneralScalarPotentialCandidates[model_Association] := Module[{d1, d2},
  d1 = model["Scalar1", "SU2"];
  d2 = model["Scalar2", "SU2"];
Join[
    BuildScalarSelfInvariantCandidates[1, d1],
    BuildScalarSelfInvariantCandidates[2, d2],
    BuildHiggsPortalInvariantCandidates[1, d1],
    BuildHiggsPortalInvariantCandidates[2, d2],
    BuildCrossScalarInvariantCandidates[d1, d2]
  ]
];

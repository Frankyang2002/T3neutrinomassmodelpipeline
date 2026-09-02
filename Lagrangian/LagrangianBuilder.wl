(* TLDR: 
- Take in SU2 and U1 representations and change into Matchete fields. 
- Create all gauge invariant interactions
- Return UV Lagrangian

Note that we use invariant tensor function from group magic in matchete. 
We can explicitly tell it the symmetric fields so that we get the tensor that we can SU2 contracts fields correctly
It returns explicit CG tensors in a list so like HH with 2x2=3+1 would be only the CG for the symmetric part {H1,H2,H3} 
*)
ClearAll[
  SU2DynkinLabel, SU2RepresentationName, EnsureSU2Representation,
  BSMRepresentationName, EnsureBSMRepresentation, BSMIndexType,
  SU2IndexType, SU2CGIndexRepresentation, DefineSU2InvariantCG,
  DefineT3InvariantCGs, T3LegacyModelQ,
  DefineT3TwoScalarFields, DefineT3TwoScalarCouplings,
  ScalarFieldValue, ScalarNorm, ScalarSelf, HiggsPortal, CrossScalarPortal,
  TensorVector, SwapTensorSlots, SymmetrizeTensorPairs, IndependentTensorBasis,
  PhysicalSU2InvariantBasis, DefineDynamicRealCoupling, DefineSU2InvariantFamily,
  ScalarPotentialFieldValue, BuildScalarSelfInvariantCandidates,
  BuildHiggsPortalInvariantCandidates, BuildCrossScalarInvariantCandidates,
  BuildGeneralScalarPotentialCandidates,
  BuildT3YukawaCandidates, BuildT3MixingCandidates,
  BuildT3YukawaCandidatesLegacy, BuildT3MixingCandidatesLegacy,
  BuildT3InteractionCandidates, ValidateT3Candidate,
  SelectFirstValidByGroup, T3IngredientsPresentQ, BuildT3Lagrangian
];

(* ------------------------------------------------- *)
(* Getting our SU2 and BSM representations           *)
(* ------------------------------------------------- *)

(* For SU(2), the irrep of dimension d has highest-weight Dynkin label {d-1}*)
SU2DynkinLabel[d_Integer?Positive] := {d - 1};

(* We just get the representation name, generalised for higher dim *)
SU2RepresentationName[1] := None;
SU2RepresentationName[2] := fund;
SU2RepresentationName[3] := adj;
SU2RepresentationName[d_Integer?Positive] /; d >= 4 :=
  Symbol["T3SU2d" <> ToString[d]];

(* We already know our first 3 have SU2 Reps, but we need to check it for higher reps *)
EnsureSU2Representation[1] := True;
EnsureSU2Representation[2] := True;
EnsureSU2Representation[3] := True;
EnsureSU2Representation[d_Integer?Positive] /; d >= 4 := Module[
  {rep, existing, after, status},
  rep = SU2RepresentationName[d]; (* Create Representation Name *)
  existing = Quiet@Check[GetRepresentations[], <||>]; (* Check if representation is already registered *)
  If[AssociationQ[existing] && KeyExistsQ[existing, rep], Return[True]]; (* If registered do nothing *)

  (* If not registered, we pair this representation with the dynkin label d and register *)
  status = Check[
    DefineRepresentation[
      rep,
      SU2L,
      SU2DynkinLabel[d],
      IndexAlphabet -> {"u", "v", "w", "x", "y", "z"}
    ];
    True,
    $Failed
  ];
  If[status === $Failed, Return[$Failed]];

  (* We check if the representation is registered *)
  after = Quiet@Check[GetRepresentations[], <||>];
  If[AssociationQ[after] && KeyExistsQ[after, rep], True, $Failed]
];

(* For matchete our SM fields remain normal but BSM fields obtain its own representation to make CG work, 
we have 2 representations, one for SM doublet basis and another for BSM representation basis
This is also true for d=2 and d=3*)


(* Get BSMrep name *)
BSMRepresentationName[1] := None;
BSMRepresentationName[d_Integer?Positive] /; d >= 2 :=
  Symbol["T3BSMd" <> ToString[d]];

(* We log our representation with its respective Dynkin label *)
EnsureBSMRepresentation[1] := True;
EnsureBSMRepresentation[d_Integer?Positive] /; d >= 2 := Module[
  {rep, before, after, defineStatus},

  rep = BSMRepresentationName[d];

  before = Quiet@Check[GetRepresentations[], <||>];
  If[AssociationQ[before] && KeyExistsQ[before, rep], Return[True]];

  defineStatus = Check[
    DefineRepresentation[
      rep,
      SU2L,
      SU2DynkinLabel[d],
      IndexAlphabet -> {"u", "v", "w", "x", "y", "z"}
    ];
    True,
    $Failed
  ];

  If[defineStatus === $Failed,
    Print["  FAILED to define SU(2) representation d=", d,
      " (Dynkin ", SU2DynkinLabel[d], ")"];
    Return[$Failed]
  ];

  (* Never claim success unless Matchete actually registered the irrep. *)
  after = Quiet@Check[GetRepresentations[], <||>];
  If[AssociationQ[after] && KeyExistsQ[after, rep],
    Print["  Defined BSM SU(2) representation ", rep,
      " with d=", d, ", Dynkin=", SU2DynkinLabel[d]];
    True,
    Print["  ERROR: DefineRepresentation returned without registering ", rep];
    $Failed
  ]
];

(* Dimension is now an index that relates our dimension to our logged dynkin labels *)
BSMIndexType[1] := None;
BSMIndexType[d_Integer?Positive] := Module[{ok = EnsureBSMRepresentation[d]},
  If[ok === $Failed, Return[Missing["UnsupportedSU2", d]]];
  BSMRepresentationName[d]
];

(* Same thing, but for generic SU2, like for d=2 we have fund, for d=3 we have adj *)
SU2IndexType[1] := None;
SU2IndexType[d_Integer?Positive] := Module[{ok = EnsureSU2Representation[d]},
  If[ok === $Failed, Return[Missing["UnsupportedSU2", d]]];
  SU2L[SU2RepresentationName[d]]
];

(* All our classified ABCDE graphs have 3 or less dimensions, we check if any of them are these classified ones *)
T3LegacyModelQ[model_Association] :=
  Max[model["Scalar1", "SU2"], model["Scalar2", "SU2"], model["Fermion", "SU2"]] <= 3;

(* Return BSM representation from dynkin index *)
SU2CGIndexRepresentation[d_Integer?Positive, objectConjugated_] := Module[{rep},
  If[d === 1, Return[None]];
  rep = BSMRepresentationName[d];
  (* For custom representations, DefineRepresentation registers the symbol
     rep itself (e.g. T3BSMd4) as the representation/index type. *)
  rep
];


(* Given field and SU(2) rep, get invariant Clebsch Gordon tensor to contract indices into SU(2) singlet from Matchete 
Our input is the Clebsch Gordon name, like T3Y1CG, whether fields in terms such as LFS is conjugated, so True,True,False would mean L is conj, F is conj, and S isnt
symmetricPositions tells us if we have symmetry in term, like for HHSS^\dagger, the 2 S terms are symmetric, thus removing antisymmetric spaces
smPositions tells us what position a field is a standard model, where we can use in built Matchete representation for them instead of BSM
*)

(*Essentially: Dimensions -> Dynkin Reps -> Invariant Tensors -> CG that can be used in Matchete*)
DefineSU2InvariantCG[
  cgName_Symbol,
  dims_List,
  objectConjugated_List,
  symmetricPositions_: {},
  smPositions_: {}
] := Module[
  {keep, keptDims, keptConj, algebraReps, cgReps, tensors, symmetryAfterDrop, positionMap},


  If[Length[dims] =!= Length[objectConjugated], Return[$Failed]]; (* For each field we need to say if it is conjugated or not *)
  If[!AllTrue[dims, IntegerQ[#] && Positive[#] &], Return[$Failed]]; (* Only integer positive dimensions *)

  (* For each field we have, we need to make sure we have a representation for its dimension, also deleting repeated *)
  Scan[EnsureBSMRepresentation, DeleteDuplicates[dims]]; 

  (* Remove singlets as they dont have indices, create new arrays with what is left over *)
  keep = Flatten@Position[dims, _?(# > 1 &)]; 
  keptDims = dims[[keep]];
  keptConj = objectConjugated[[keep]];

  (* No CG if we have a singlet *)
  If[keptDims === {}, Return[None]];

  (* Matchete distinguishes the orientation of pseudoreal (half integer isospin and even dimensions) SU(2) indices.
       - a conjugated pseudoreal FIELD is represented by CRep[label] in
         InvariantTensors and by an UNBARRED representation in DefineCG;
       - an unconjugated pseudoreal FIELD uses the plain Dynkin label in
         InvariantTensors and the BARRED representation in DefineCG.
     For odd-dimensional full integer isospin SU(2) irreps the representation is real, so no CRep/Bar distinction is required.
    Objectconjugated/keptConj tells us if the field object is gauge-conjugated *)

  (* The following uses this, where unconjugated becomes {3} and conjugated will be CRep[{3}], this is what we need to do for Machete *)
  algebraReps = MapThread[
    Function[{d, conjugated},
      Module[{label = SU2DynkinLabel[d]},
        If[EvenQ[d] && TrueQ[conjugated],
          CRep[label],
          label
        ]
      ]
    ],
    {keptDims, keptConj}
  ];

  (* Matchete's Invariant Tensors returns tensors with indices that transform in conjugate representation 
    So tensor index transforms opposite to the algebra field
  *)
  cgReps = MapThread[
    Function[{originalPosition, d, conjugated},
      Module[{baseRep},
        baseRep = If[MemberQ[smPositions, originalPosition],
          (* SM doublets use Matchete's built-in fundamental representation. *)
          SU2L[fund],
          SU2CGIndexRepresentation[d, conjugated]
        ];

        If[EvenQ[d] && !TrueQ[conjugated],
          Bar[baseRep],
          baseRep
        ]
      ]
    ],
    {keep, keptDims, keptConj}
  ];

  (* If we removes a singlet, we need to shift all our positions to fit the new positions of our arrays, like symmetry positions *)
  positionMap = AssociationThread[keep -> Range[Length[keep]]];


  (* We want to prevent antisymmetric tensor products when we have symmetric fields, this is used later *)
  symmetryAfterDrop = Select[
    Lookup[positionMap, #, Missing["Dropped"]] & /@ symmetricPositions,
    IntegerQ
  ];

  Print["  CG ", SymbolName[cgName], ": dims=", dims,
    ", oriented algebra reps=", InputForm[algebraReps],
    ", CG reps=", InputForm[cgReps],
    ", symmetric positions=", symmetryAfterDrop];

  (* We obtain Invariant tensors, noting that symmetry is tracked *)
  tensors = Check[
    If[symmetryAfterDrop === {},
      InvariantTensors[SU[2], algebraReps],
      InvariantTensors[SU[2], algebraReps, SymmetricIndices -> symmetryAfterDrop]
    ],
    $Failed
  ];

  Print["    InvariantTensors result head: ", Head[tensors],
    If[ListQ[tensors], ", count=" <> ToString[Length[tensors]], ""]];

  (* If no tensors we failed *)
  If[tensors === $Failed || !ListQ[tensors] || Length[tensors] === 0,
    Print["    FAILED while generating invariant tensor for ", SymbolName[cgName]];
    Return[$Failed]
  ];

  (* If we have more than one invariant tensor in large representation, we just use the first basis *)
  If[Length[tensors] > 1,
    Print["WARNING: ", SymbolName[cgName], " has ", Length[tensors],
      " invariant tensors; using the first basis tensor."]
  ];

  Module[{defined},
    defined = Check[
      DefineCG[cgName, cgReps, First[tensors]]; (* < We get the CG coefficients from invariant tensors *)
      cgName,
      $Failed
    ];
    Print["    DefineCG result: ", InputForm[defined]];
    defined
  ]
];

(* Register the three topology-defining invariant tensors for our 3 fields in the model.  Their field order is the same order used in the interaction. *)
(* We get invariant CG from each our interactions vertices, and define it for the model *)
DefineT3InvariantCGs[model_Association] := Module[
  {d1, d2, dF, y1cg, y2cg, mixcg},
  
  d1 = model["Scalar1", "SU2"];
  d2 = model["Scalar2", "SU2"];
  dF = model["Fermion", "SU2"];

  (* Bar[L] F^c S1 *)
  y1cg = DefineSU2InvariantCG[
    T3Y1CG, {2, dF, d1}, {True, True, False}, {}, {1}
  ];

  (* Bar[L] F S2^dagger *)
  y2cg = DefineSU2InvariantCG[
    T3Y2CG, {2, dF, d2}, {True, False, True}, {}, {1}
  ];

  (* H H S1 S2^dagger.  The two H fields are identical bosons, so force the
     symmetric H-H channel. *)
  mixcg = DefineSU2InvariantCG[
    T3MixCG,
    {2, 2, d1, d2},
    {False, False, False, True},
    {1, 2},
    {1, 2}
  ];

  If[MemberQ[{y1cg, y2cg, mixcg}, $Failed], Return[$Failed]];
  <|"Yukawa1" -> y1cg, "Yukawa2" -> y2cg, "Mixing" -> mixcg|>
];

(* We define couplings for matchete *)
DefineT3TwoScalarCouplings[] := Module[{},
  DefineCoupling[y1, Indices -> {Flavor, NFlavor}, SelfConjugate -> False];
  DefineCoupling[y2, Indices -> {Flavor, NFlavor}, SelfConjugate -> False];
  DefineCoupling[lambdaS1, SelfConjugate -> True]; (* (S_1^\daggerS_1)^2 *)
  DefineCoupling[lambdaS2, SelfConjugate -> True]; (* (S_2^\daggerS_2)^2 *)
  DefineCoupling[lambdaH1, SelfConjugate -> True]; (* (H^\daggerH)(S_1^\daggerS_1) *)
  DefineCoupling[lambdaH2, SelfConjugate -> True]; (* (H^\daggerH)(S_2^\daggerS_2) *)
  DefineCoupling[lambda12, SelfConjugate -> True]; (* (S_1^\daggerS_1)(S_2^\daggerS_2) *)
  DefineCoupling[lambdaT3, SelfConjugate -> False]; (* (HH)(S_1S_2^\dagger) +h.c. *)
  True
];

(* We convert into matchete fields *)
DefineT3TwoScalarFields[model_Association] := Module[
  {f, s1, s2, fidx, s1idx, s2idx, findices, s1indices, s2indices, selfConj},

  f = model["Fermion"];
  s1 = model["Scalar1"];
  s2 = model["Scalar2"];

  (* Define every higher SU(2) representation before it is used as an index. *)
  If[AnyTrue[
      EnsureBSMRepresentation /@ DeleteDuplicates[{f["SU2"], s1["SU2"], s2["SU2"]}],
      # === $Failed &
    ],
    Return[$Failed]
  ];

  fidx = If[T3LegacyModelQ[model], SU2IndexType[f["SU2"]], BSMIndexType[f["SU2"]]];
  s1idx = If[T3LegacyModelQ[model], SU2IndexType[s1["SU2"]], BSMIndexType[s1["SU2"]]];
  s2idx = If[T3LegacyModelQ[model], SU2IndexType[s2["SU2"]], BSMIndexType[s2["SU2"]]];
  If[AnyTrue[{fidx, s1idx, s2idx}, MissingQ], Return[$Failed]];

  DefineFlavorIndex[
    NFlavor,
    Lookup[f, "Multiplicity", 3],
    IndexAlphabet -> {"r", "s", "t"}
  ];

  findices = If[fidx === None, {NFlavor}, {fidx, NFlavor}];
  s1indices = If[s1idx === None, {}, {s1idx}];
  s2indices = If[s2idx === None, {}, {s2idx}];

  (* Integer-isospin SU(2) irreps (odd dimension) are real.  A neutral
     fermion in such a representation may therefore be self-conjugate. *)
  selfConj = TrueQ[PossibleZeroQ[f["Y"]]] && OddQ[f["SU2"]];

  If[selfConj,
    DefineField[
      NewFermion, Fermion,
      SelfConjugate -> True,
      Indices -> findices,
      Mass -> {Heavy, Lookup[f, "MassSymbol", MF]}
    ],
    DefineField[
      NewFermion, Fermion,
      SelfConjugate -> False,
      Indices -> findices,
      Charges -> {U1Y[f["Y"]]},
      Mass -> {Heavy, Lookup[f, "MassSymbol", MF]}
    ]
  ];

  DefineField[
    NewScalar1, Scalar,
    SelfConjugate -> False,
    Indices -> s1indices,
    Charges -> {U1Y[s1["Y"]]},
    Mass -> {Heavy, Lookup[s1, "MassSymbol", MS1]}
  ];

  DefineField[
    NewScalar2, Scalar,
    SelfConjugate -> False,
    Indices -> s2indices,
    Charges -> {U1Y[s2["Y"]]},
    Mass -> {Heavy, Lookup[s2, "MassSymbol", MS2]}
  ];

  True
];

ScalarFieldValue[1, False, ___] := NewScalar1[];
ScalarFieldValue[1, True,  ___] := Bar[NewScalar1[]];
ScalarFieldValue[2, False, ___] := NewScalar2[];
ScalarFieldValue[2, True,  ___] := Bar[NewScalar2[]];
ScalarFieldValue[1, False, i_] := NewScalar1[i];
ScalarFieldValue[1, True,  i_] := Bar[NewScalar1[i]];
ScalarFieldValue[2, False, i_] := NewScalar2[i];
ScalarFieldValue[2, True,  i_] := Bar[NewScalar2[i]];

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
  smPositions_: {},
  legacyQ_: False
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
          If[TrueQ[legacyQ],
            SU2IndexType[d],
            SU2CGIndexRepresentation[d, conjugated]
          ]
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
  smPositions_: {},
  legacyQ_: False
] := Module[
  {data, basis, cgReps, keep, cgNames, couplings, cgName, couplingName,
   cgDefined, couplingDefined},

  data = PhysicalSU2InvariantBasis[
    dims,
    objectConjugated,
    symmetricPairs,
    smPositions,
    legacyQ
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
  d_Integer?Positive,
  legacyQ_: False
] := Module[
  {family, indices, labels, cg, coupling, fieldProduct, baseName},

  baseName = "lambdaS" <> ToString[which];

  family = DefineSU2InvariantFamily[
    "T3Scalar" <> ToString[which] <> "SelfCG",
    baseName,
    {d, d, d, d},
    {True, False, True, False},
    {{1, 3}, {2, 4}},
    {},
    legacyQ
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
  d_Integer?Positive,
  legacyQ_: False
] := Module[
  {family, indices, labels, cg, coupling, fieldProduct, baseName},

  baseName = "lambdaH" <> ToString[which];

  family = DefineSU2InvariantFamily[
    "T3HiggsPortal" <> ToString[which] <> "CG",
    baseName,
    {2, 2, d, d},
    {True, False, True, False},
    {},
    {1, 2},
    legacyQ
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
  d2_Integer?Positive,
  legacyQ_: False
] := Module[
  {family, indices, labels, cg, coupling, fieldProduct},

  family = DefineSU2InvariantFamily[
    "T3CrossScalarCG",
    "lambda12",
    {d1, d1, d2, d2},
    {True, False, True, False},
    {},
    {},
    legacyQ
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
BuildGeneralScalarPotentialCandidates[model_Association] := Module[{d1, d2, legacyQ},
  d1 = model["Scalar1", "SU2"];
  d2 = model["Scalar2", "SU2"];
  legacyQ = T3LegacyModelQ[model];

  Join[
    BuildScalarSelfInvariantCandidates[1, d1, legacyQ],
    BuildScalarSelfInvariantCandidates[2, d2, legacyQ],
    BuildHiggsPortalInvariantCandidates[1, d1, legacyQ],
    BuildHiggsPortalInvariantCandidates[2, d2, legacyQ],
    BuildCrossScalarInvariantCandidates[d1, d2, legacyQ]
  ]
];


(* ------------------------------------------------------------------------- *)
(* Legacy d<=3 contractions
   -------------------------------------------------------------------------
   These are the exact contractions used by the previously validated A-E
   pipeline.  Keeping them as the low-dimensional regression oracle prevents
   the arbitrary-irrep implementation from silently changing established
   physics conventions. *)
BuildT3YukawaCandidatesLegacy[model_Association, which_Integer] := Module[
  {s, ds, df, y, scalarBar, useCConj, i,j,k,I,p,r, f1, baseName},
  s = model[If[which===1,"Scalar1","Scalar2"]];
  ds=s["SU2"]; df=model["Fermion","SU2"];
  useCConj = which === 1; scalarBar = which === 2;
  y = If[which===1,y1,y2]; baseName = "Yukawa"<>ToString[which];
  f1[args___] := If[useCConj, CConj[NewFermion[args]], NewFermion[args]];

  Switch[{ds,df},
    {1,2},
      {
        <|"Name"->baseName<>"_Direct", "Class"->baseName, "AlternativeGroup"->baseName,
          "Expression"->PlusHc[y[p,r] NCM[Bar[l[i,p]], f1[i,r]] ScalarFieldValue[which,scalarBar]]|>,
        <|"Name"->baseName<>"_Eps", "Class"->baseName, "AlternativeGroup"->baseName,
          "Expression"->PlusHc[y[p,r] NCM[Bar[l[i,p]], f1[j,r]] CG[eps[SU2L],{i,j}] ScalarFieldValue[which,scalarBar]]|>
      },
    {2,1},
      {
        <|"Name"->baseName<>"_Eps", "Class"->baseName, "AlternativeGroup"->baseName,
          "Expression"->PlusHc[y[p,r] NCM[Bar[l[i,p]],f1[r]] CG[eps[SU2L],{i,j}] ScalarFieldValue[which,scalarBar,j]]|>,
        <|"Name"->baseName<>"_Delta", "Class"->baseName, "AlternativeGroup"->baseName,
          "Expression"->PlusHc[y[p,r] NCM[Bar[l[i,p]],f1[r]] ScalarFieldValue[which,scalarBar,i]]|>
      },
    {2,3},
      {
        <|"Name"->baseName<>"_GenEps", "Class"->baseName, "AlternativeGroup"->baseName,
          "Expression"->PlusHc[y[p,r] NCM[Bar[l[i,p]],f1[I,r]]
            CG[gen[SU2L[fund]],{I,i,k}] CG[eps[SU2L],{k,j}]
            ScalarFieldValue[which,scalarBar,j]]|>,
        <|"Name"->baseName<>"_Gen", "Class"->baseName, "AlternativeGroup"->baseName,
          "Expression"->PlusHc[y[p,r] NCM[Bar[l[i,p]],f1[I,r]]
            CG[gen[SU2L[fund]],{I,i,j}] ScalarFieldValue[which,scalarBar,j]]|>
      },
    {3,2},
      {
        <|"Name"->baseName<>"_Gen", "Class"->baseName, "AlternativeGroup"->baseName,
          "Expression"->PlusHc[y[p,r] NCM[Bar[l[i,p]],f1[j,r]]
            CG[gen[SU2L[fund]],{I,i,j}] ScalarFieldValue[which,scalarBar,I]]|>,
        <|"Name"->baseName<>"_EpsGen", "Class"->baseName, "AlternativeGroup"->baseName,
          "Expression"->PlusHc[y[p,r] NCM[Bar[l[i,p]],f1[j,r]]
            CG[eps[SU2L],{j,k}] CG[gen[SU2L[fund]],{I,i,k}]
            ScalarFieldValue[which,scalarBar,I]]|>
      },
    _, $Failed
  ]
];

BuildT3MixingCandidatesLegacy[model_Association] := Module[
  {d1=model["Scalar1","SU2"], d2=model["Scalar2","SU2"], i,j,k,A,B,C},
  Switch[{d1,d2},
    {2,2}, {
      <|"Name"->"T3ScalarMix_Doublets", "Class"->"T3ScalarMix",
        "AlternativeGroup"->"T3ScalarMix",
        "Expression"->PlusHc[lambdaT3[] H[i] Bar[NewScalar2[i]]
          CG[Bar[eps[SU2L]],{Bar[j],Bar[k]}] H[j] NewScalar1[k]]|>
    },
    {1,3}, {
      <|"Name"->"T3ScalarMix_1x3", "Class"->"T3ScalarMix",
        "AlternativeGroup"->"T3ScalarMix",
        "Expression"->PlusHc[lambdaT3[] NewScalar1[] Bar[NewScalar2[A]] H[i] H[j]
          CG[gen[SU2L[fund]],{A,k,Bar[i]}]
          CG[Bar[eps[SU2L]],{Bar[k],Bar[j]}]]|>
    },
    {3,1}, {
      <|"Name"->"T3ScalarMix_3x1", "Class"->"T3ScalarMix",
        "AlternativeGroup"->"T3ScalarMix",
        "Expression"->PlusHc[lambdaT3[] NewScalar1[A] Bar[NewScalar2[]] H[i] H[j]
          CG[gen[SU2L[fund]],{A,k,Bar[i]}]
          CG[Bar[eps[SU2L]],{Bar[k],Bar[j]}]]|>
    },
    {3,3}, {
      <|"Name"->"T3ScalarMix_3x3", "Class"->"T3ScalarMix",
        "AlternativeGroup"->"T3ScalarMix",
        "Expression"->PlusHc[lambdaT3[] NewScalar1[A] Bar[NewScalar2[B]] H[i] H[j]
          CG[gen[SU2L[fund]],{C,k,Bar[i]}]
          CG[Bar[eps[SU2L]],{Bar[k],Bar[j]}]
          CG[fStruct[SU2L],{A,B,C}]]|>
    },
    _, {}
  ]
];

(* Build the topology Yukawa with the generic Matchete-generated CG. Remember its yukawa *)
BuildT3YukawaCandidates[model_Association, which_Integer] := Module[
  {s, dS, dF, y, scalarBar, useCConj, baseName, p, r, i, j, k, f, scalar, labels, cg},
  (* Eg: 2x1x2 -> {i,k}, and 2x4x3 -> {i,j,k}  due to singlets getting out*)
  s = model[If[which === 1, "Scalar1", "Scalar2"]];
  dS = s["SU2"];
  dF = model["Fermion", "SU2"];
  useCConj = which === 1;
  scalarBar = which === 2;
  y = If[which === 1, y1, y2];
  baseName = "Yukawa" <> ToString[which];
  cg = If[which === 1, T3Y1CG, T3Y2CG];

  (* Checks if our fermion and scalar have an SU(2) index, so not singlet *)
  f = If[dF === 1,
    If[useCConj, CConj[NewFermion[r]], NewFermion[r]],
    If[useCConj, CConj[NewFermion[j, r]], NewFermion[j, r]]
  ];

  scalar = If[dS === 1,
    ScalarFieldValue[which, scalarBar],
    ScalarFieldValue[which, scalarBar, k]
  ];

  labels = Join[
    {i},
    If[dF === 1, {}, {j}],
    If[dS === 1, {}, {k}]
  ];

  {
    <|
      "Name" -> baseName <> "_GenericCG",
      "Class" -> baseName,
      "AlternativeGroup" -> baseName,
      "Expression" -> PlusHc[
        y[p, r] NCM[Bar[l[i, p]], f] scalar CG[cg, labels]
      ]
    |>
  }
];

(* Now H H S1 S2^dagger + h.c. vertex  The invariant tensor is
   generated for the exact requested representations and is symmetric in the
   two Higgs indices. *)
BuildT3MixingCandidates[model_Association] := Module[
  {d1, d2, i, j, a, b, s1, s2, labels},

  d1 = model["Scalar1", "SU2"];
  d2 = model["Scalar2", "SU2"];

  s1 = If[d1 === 1, NewScalar1[], NewScalar1[a]];
  s2 = If[d2 === 1, Bar[NewScalar2[]], Bar[NewScalar2[b]]];
  labels = Join[{i, j}, If[d1 === 1, {}, {a}], If[d2 === 1, {}, {b}]];

  {
    <|
      "Name" -> "T3ScalarMix_GenericCG",
      "Class" -> "T3ScalarMix",
      "AlternativeGroup" -> "T3ScalarMix",
      "Expression" -> PlusHc[
        lambdaT3[] H[i] H[j] s1 s2 CG[T3MixCG, labels]
      ]
    |>
  }
];

(* We check with Matchete if the lagrangian with the interaction is valid *)
ValidateT3Candidate[c_, LSM_, LFree_] := Module[{res},
  Print["Checking candidate: ", c["Name"]];
  (* Relabel indices is to prevent dummy index clashes, 
  like for A_i B_i + C_i D_i  are independent despite both using i 
  We prevent reusing index labels and throughout the file*)
  res = CheckAbort[
    Check[CheckLagrangian[(LSM + LFree + c["Expression"]) // RelabelIndices], $Failed],
    $Aborted
  ];
  Print["  result: ", InputForm[res]];
  Association[c, "Valid" -> TrueQ[res], "ValidationResult" -> res]
];

(* Keep the first valid contraction in each alternative group.  Candidates
   without a group are independent interactions and are all retained. *)
SelectFirstValidByGroup[list_List] := Module[{seen = <||>},
  Select[
    list,
    Function[candidate,
      Module[{group = Lookup[candidate, "AlternativeGroup", None]},
        If[group === None,
          True,
          If[KeyExistsQ[seen, group], False, seen[group] = True; True]
        ]
      ]
    ]
  ]
];

(* Assemble all interaction candidates without validating them.*)
BuildT3InteractionCandidates[model_Association, legacyQ_] := Module[{d1, d2},
  d1 = model["Scalar1", "SU2"];
  d2 = model["Scalar2", "SU2"];

  Join[
    If[legacyQ,
      BuildT3YukawaCandidatesLegacy[model, 1],
      BuildT3YukawaCandidates[model, 1]
    ],
    If[legacyQ,
      BuildT3YukawaCandidatesLegacy[model, 2],
      BuildT3YukawaCandidates[model, 2]
    ],
    BuildGeneralScalarPotentialCandidates[model],
    If[legacyQ,
      BuildT3MixingCandidatesLegacy[model],
      BuildT3MixingCandidates[model]
    ]
  ]
];

(* A genuine T3 contribution requires both lepton-fermion-scalar Yukawas and
   the quartic that connects S1 and S2 to the two Higgs legs. 
   No weinberg without these 3 elements *)
T3IngredientsPresentQ[valid_List] := And @@ (
  Function[class,
    AnyTrue[valid, Lookup[#, "Class", ""] === class &]
  ] /@ {"Yukawa1", "Yukawa2", "T3ScalarMix"}
);

(* We do everything her  *)
BuildT3Lagrangian[model_Association] := Module[
  {LSM, LFree, cgStatus, legacyQ, candidates, validated, valid, rejected,
   LInt, LBSM, LUV, fullValidation, ingredientsPresent},

  ResetAll[];
  
  (* We check if we use the classes or different models *)
  legacyQ = T3LegacyModelQ[model];

  Print["Build stage 1/4: load SM"];
  LSM = LoadModel["SM"];
  Print["  LoadModel result: ", If[LSM === $Failed, "$Failed", "OK"]];
  If[LSM === $Failed, Return[$Failed]];

  Print["Build stage 2/4: define BSM fields"];
  If[DefineT3TwoScalarFields[model] === $Failed,
    Print["  field definition result: $Failed"];
    Return[$Failed],
    Print["  field definition result: True"]
  ];

  Print["Build stage 3/4: define couplings"];
  If[Check[DefineT3TwoScalarCouplings[], $Failed] === $Failed,
    Print["  coupling definition result: $Failed"];
    Return[$Failed],
    Print["  coupling definition result: True"]
  ];

  If[legacyQ,
    Print["Build stage 4/4: legacy d<=3 contractions (no custom CG setup)"];
    cgStatus = <|"Mode" -> "Legacy"|>,
    Print["Build stage 4/4: define generic SU(2) CG tensors"];
    cgStatus = DefineT3InvariantCGs[model];
    Print["  CG definition result: ", InputForm[cgStatus]];
    If[cgStatus === $Failed, Return[$Failed]]
  ];

  (* Get free Lagrangian *)
  LFree = FreeLag[NewFermion, NewScalar1, NewScalar2] // RelabelIndices;

  candidates = BuildT3InteractionCandidates[model, legacyQ];
  candidates = Association[#, "Expression" -> RelabelIndices[#["Expression"]]] & /@ candidates;

  (* We validate each interaction individually, seeing if they are compatible *)
  validated = ValidateT3Candidate[#, LSM, LFree] & /@ candidates;
  valid = SelectFirstValidByGroup @ Select[validated, TrueQ[#["Valid"]] &];
  rejected = Select[validated, !TrueQ[#["Valid"]] &];

  (* Build interaction lagrangian and then the UV *)
  LInt = Total[Lookup[valid, "Expression", {}]] // Expand // RelabelIndices;
  LBSM = (LFree + LInt) // Expand // RelabelIndices;
  LUV = (LSM + LBSM) // Expand // RelabelIndices;
  fullValidation = CheckAbort[Check[CheckLagrangian[LUV], $Failed], $Aborted];

  (* Check if it has all the vertices for T3 diagram *)
  ingredientsPresent = T3IngredientsPresentQ[valid];

  <|
    "Status" -> If[TrueQ[fullValidation], "Success", "FullValidationFailed"],
    "Model" -> model,
    "LSM" -> LSM,
    "LFree" -> LFree,
    "LInt" -> LInt,
    "LBSM" -> LBSM,
    "LUV" -> LUV,
    "AllowedInteractions" -> Lookup[valid, "Name", {}],
    "RejectedInteractions" -> Lookup[rejected, "Name", {}],
    "T3IngredientsPresent" -> ingredientsPresent,
    "WeinbergIngredientsPresent" -> ingredientsPresent,
    "FullValidation" -> fullValidation
  |>
];

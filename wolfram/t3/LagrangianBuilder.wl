(* LagrangianBuilder.wl
   Generic T3 UV-Lagrangian builder for arbitrary finite-dimensional SU(2)
   irreducible representations supported by Matchete.

   The old implementation explicitly wrote contractions for dimensions 1,2,3.
   This version instead:
     1. defines any required SU(2) irrep using its Dynkin label {d-1};
     2. asks Matchete for the exact invariant tensor with InvariantTensors;
     3. registers that tensor as a CG coefficient with DefineCG;
     4. uses the resulting CG in the Yukawa and scalar-mixing vertices.

   This keeps Matchete as the final judge of both covariance and matching.
*)
ClearAll[
  SU2DynkinLabel, SU2RepresentationName, EnsureSU2Representation,
  BSMRepresentationName, EnsureBSMRepresentation, BSMIndexType,
  SU2IndexType, SU2CGIndexRepresentation, DefineSU2InvariantCG,
  DefineT3InvariantCGs, T3LegacyModelQ,
  DefineT3TwoScalarFields, DefineT3TwoScalarCouplings,
  ScalarFieldValue, ScalarNorm, ScalarSelf, HiggsPortal, CrossScalarPortal,
  BuildT3YukawaCandidates, BuildT3MixingCandidates,
  BuildT3YukawaCandidatesLegacy, BuildT3MixingCandidatesLegacy,
  BuildT3InteractionCandidates, ValidateT3Candidate,
  SelectFirstValidByGroup, T3IngredientsPresentQ, BuildT3Lagrangian
];

(* For SU(2), the irrep of dimension d has highest-weight Dynkin label {d-1}. *)
SU2DynkinLabel[d_Integer?Positive] := {d - 1};

(* Reuse the SM's built-in fundamental and adjoint names.  Higher irreps get
   deterministic symbols so the same dimension is defined only once/session. *)
SU2RepresentationName[1] := None;
SU2RepresentationName[2] := fund;
SU2RepresentationName[3] := adj;
SU2RepresentationName[d_Integer?Positive] /; d >= 4 :=
  Symbol["T3SU2d" <> ToString[d]];

EnsureSU2Representation[1] := True;
EnsureSU2Representation[2] := True;
EnsureSU2Representation[3] := True;
EnsureSU2Representation[d_Integer?Positive] /; d >= 4 := Module[
  {rep, existing, after, status},
  rep = SU2RepresentationName[d];
  existing = Quiet@Check[GetRepresentations[], <||>];
  If[AssociationQ[existing] && KeyExistsQ[existing, rep], Return[True]];

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

  after = Quiet@Check[GetRepresentations[], <||>];
  If[AssociationQ[after] && KeyExistsQ[after, rep], True, $Failed]
];

(* Generic-path BSM irreps use their own representation symbols, including
   d=2 and d=3.  This is deliberate: InvariantTensors and DefineRepresentation
   then use the same canonical GroupMagic basis.  The SM L and H remain in the
   built-in SU2L[fund] representation. *)
BSMRepresentationName[1] := None;
BSMRepresentationName[d_Integer?Positive] /; d >= 2 :=
  Symbol["T3BSMd" <> ToString[d]];

EnsureBSMRepresentation[1] := True;
EnsureBSMRepresentation[d_Integer?Positive] /; d >= 2 := Module[
  {rep, before, after, defineStatus},

  rep = BSMRepresentationName[d];

  before = Quiet@Check[GetRepresentations[], <||>];
  If[AssociationQ[before] && KeyExistsQ[before, rep], Return[True]];

  (* Matchete's DefineRepresentation takes the already-defined gauge-group
     name.  For the SM SU(2)_L group this is SU2L. *)
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

BSMIndexType[1] := None;
BSMIndexType[d_Integer?Positive] := Module[{ok = EnsureBSMRepresentation[d]},
  If[ok === $Failed, Return[Missing["UnsupportedSU2", d]]];
  BSMRepresentationName[d]
];

SU2IndexType[1] := None;
SU2IndexType[d_Integer?Positive] := Module[{ok = EnsureSU2Representation[d]},
  If[ok === $Failed, Return[Missing["UnsupportedSU2", d]]];
  SU2L[SU2RepresentationName[d]]
];

(* The five historical A-E assignments use only d<=3 and retain the exact
   already-regression-tested Matchete contractions.  Any model containing a
   higher irrep uses the generic canonical-BSM path. *)
T3LegacyModelQ[model_Association] :=
  Max[model["Scalar1", "SU2"], model["Scalar2", "SU2"], model["Fermion", "SU2"]] <= 3;

(* Return the base Matchete representation/index type for a custom BSM irrep.
   Pseudoreal Bar[...] orientation is applied centrally in
   DefineSU2InvariantCG, where it can be kept consistent with CRep[...]. *)
SU2CGIndexRepresentation[d_Integer?Positive, objectConjugated_] := Module[{rep},
  If[d === 1, Return[None]];
  rep = BSMRepresentationName[d];

  (* For custom representations, DefineRepresentation registers the symbol
     rep itself (e.g. T3BSMd4) as the representation/index type. *)
  rep
];

(* Define one exact SU(2) invariant tensor.

   dims gives the representation dimension of each FIELD object appearing in
   the interaction. objectConjugated marks whether that field object itself is
   conjugated (Bar[...] or charge conjugated for gauge purposes).

   InvariantTensors[alg,reps] returns tensors whose indices transform in the
   conjugate representations, exactly what is needed to contract the fields.
   Singlets carry no Matchete index and are removed before defining the CG.
*)
DefineSU2InvariantCG[
  cgName_Symbol,
  dims_List,
  objectConjugated_List,
  symmetricPositions_: {},
  smPositions_: {}
] := Module[
  {keep, keptDims, keptConj, algebraReps, cgReps, tensors, symmetryAfterDrop, positionMap},

  If[Length[dims] =!= Length[objectConjugated], Return[$Failed]];
  If[!AllTrue[dims, IntegerQ[#] && Positive[#] &], Return[$Failed]];

  Scan[EnsureBSMRepresentation, DeleteDuplicates[dims]];

  keep = Flatten@Position[dims, _?(# > 1 &)];
  keptDims = dims[[keep]];
  keptConj = objectConjugated[[keep]];

  (* No non-trivial indices: the interaction is already an SU(2) singlet. *)
  If[keptDims === {}, Return[None]];

  (* Matchete distinguishes the orientation of pseudoreal SU(2) indices.

     The orientation rule was verified explicitly for 2 x 4 x 3:
       - a conjugated pseudoreal FIELD is represented by CRep[label] in
         InvariantTensors and by an UNBARRED representation in DefineCG;
       - an unconjugated pseudoreal FIELD uses the plain Dynkin label in
         InvariantTensors and the BARRED representation in DefineCG.

     These choices are complementary.  For odd-dimensional SU(2) irreps the
     representation is real, so no CRep/Bar distinction is required.

     objectConjugated therefore has a direct physics meaning here: it records
     whether the field object appearing in the interaction is gauge-conjugated. *)
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

  (* Translate symmetry positions from the original field list to the reduced
     list after SU(2) singlets have been removed. *)
  positionMap = AssociationThread[keep -> Range[Length[keep]]];
  symmetryAfterDrop = Select[
    Lookup[positionMap, #, Missing["Dropped"]] & /@ symmetricPositions,
    IntegerQ
  ];

  Print["  CG ", SymbolName[cgName], ": dims=", dims,
    ", oriented algebra reps=", InputForm[algebraReps],
    ", CG reps=", InputForm[cgReps],
    ", symmetric positions=", symmetryAfterDrop];

  (* Do not Quiet this during generic-representation bring-up: Matchete's
     diagnostic message is essential if a representation/contraction fails. *)
  tensors = Check[
    If[symmetryAfterDrop === {},
      InvariantTensors[SU[2], algebraReps],
      InvariantTensors[SU[2], algebraReps, SymmetricIndices -> symmetryAfterDrop]
    ],
    $Failed
  ];

  Print["    InvariantTensors result head: ", Head[tensors],
    If[ListQ[tensors], ", count=" <> ToString[Length[tensors]], ""]];

  If[tensors === $Failed || !ListQ[tensors] || Length[tensors] === 0,
    Print["    FAILED while generating invariant tensor for ", SymbolName[cgName]];
    Return[$Failed]
  ];

  (* Every T3 invariant used here is one-dimensional once the identical-Higgs
     symmetry is imposed.  If Matchete ever finds more than one tensor, keep
     the first as a definite basis choice and expose the multiplicity. *)
  If[Length[tensors] > 1,
    Print["WARNING: ", SymbolName[cgName], " has ", Length[tensors],
      " invariant tensors; using the first basis tensor."]
  ];

  Module[{defined},
    defined = Check[
      DefineCG[cgName, cgReps, First[tensors]];
      cgName,
      $Failed
    ];
    Print["    DefineCG result: ", InputForm[defined]];
    defined
  ]
];

(* Register the three topology-defining invariant tensors for the current
   model.  Their field order is the same order used in the interaction. *)
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

DefineT3TwoScalarCouplings[] := Module[{},
  DefineCoupling[y1, Indices -> {Flavor, NFlavor}, SelfConjugate -> False];
  DefineCoupling[y2, Indices -> {Flavor, NFlavor}, SelfConjugate -> False];
  DefineCoupling[lambdaS1, SelfConjugate -> True];
  DefineCoupling[lambdaS2, SelfConjugate -> True];
  DefineCoupling[lambdaH1, SelfConjugate -> True];
  DefineCoupling[lambdaH2, SelfConjugate -> True];
  DefineCoupling[lambda12, SelfConjugate -> True];
  DefineCoupling[lambdaT3, SelfConjugate -> False];
  True
];

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

(* Build the topology Yukawa with the generic Matchete-generated CG. *)
BuildT3YukawaCandidates[model_Association, which_Integer] := Module[
  {s, dS, dF, y, scalarBar, useCConj, baseName, p, r, i, j, k, f, scalar, labels, cg},

  s = model[If[which === 1, "Scalar1", "Scalar2"]];
  dS = s["SU2"];
  dF = model["Fermion", "SU2"];
  useCConj = which === 1;
  scalarBar = which === 2;
  y = If[which === 1, y1, y2];
  baseName = "Yukawa" <> ToString[which];
  cg = If[which === 1, T3Y1CG, T3Y2CG];

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

(* Essential T3 scalar vertex H H S1 S2^dagger + h.c.  The invariant tensor is
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

ValidateT3Candidate[c_, LSM_, LFree_] := Module[{res},
  Print["Checking candidate: ", c["Name"]];
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

(* Assemble all interaction candidates without validating them.  Keeping this
   separate from BuildT3Lagrangian makes the physics content easy to inspect:
   two topology Yukawas, the scalar potential terms, and the T3-closing quartic. *)
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
    {
      <|"Name" -> "Scalar1Self", "Class" -> "ScalarSelf",
        "Expression" -> ScalarSelf[1, d1, lambdaS1]|>,
      <|"Name" -> "Scalar2Self", "Class" -> "ScalarSelf",
        "Expression" -> ScalarSelf[2, d2, lambdaS2]|>,
      <|"Name" -> "HiggsPortal1", "Class" -> "Portal",
        "Expression" -> HiggsPortal[1, d1, lambdaH1]|>,
      <|"Name" -> "HiggsPortal2", "Class" -> "Portal",
        "Expression" -> HiggsPortal[2, d2, lambdaH2]|>,
      <|"Name" -> "ScalarCrossPortal", "Class" -> "Portal",
        "Expression" -> CrossScalarPortal[d1, d2]|>
    },
    If[legacyQ,
      BuildT3MixingCandidatesLegacy[model],
      BuildT3MixingCandidates[model]
    ]
  ]
];

(* A genuine T3 contribution requires both lepton-fermion-scalar Yukawas and
   the quartic that connects S1 and S2 to the two Higgs legs. *)
T3IngredientsPresentQ[valid_List] := And @@ (
  Function[class,
    AnyTrue[valid, Lookup[#, "Class", ""] === class &]
  ] /@ {"Yukawa1", "Yukawa2", "T3ScalarMix"}
);

BuildT3Lagrangian[model_Association] := Module[
  {LSM, LFree, cgStatus, legacyQ, candidates, validated, valid, rejected,
   LInt, LBSM, LUV, fullValidation, ingredientsPresent},

  ResetAll[];
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

  (* Custom CGs can only be registered after LoadModel has created SU2L and
     after the BSM representations have been defined.  The historical d<=3
     models deliberately keep their already-regression-tested contractions. *)
  If[legacyQ,
    Print["Build stage 4/4: legacy d<=3 contractions (no custom CG setup)"];
    cgStatus = <|"Mode" -> "Legacy"|>,
    Print["Build stage 4/4: define generic SU(2) CG tensors"];
    cgStatus = DefineT3InvariantCGs[model];
    Print["  CG definition result: ", InputForm[cgStatus]];
    If[cgStatus === $Failed, Return[$Failed]]
  ];

  LFree = FreeLag[NewFermion, NewScalar1, NewScalar2] // RelabelIndices;

  candidates = BuildT3InteractionCandidates[model, legacyQ];
  candidates = Association[#, "Expression" -> RelabelIndices[#["Expression"]]] & /@ candidates;

  validated = ValidateT3Candidate[#, LSM, LFree] & /@ candidates;
  valid = SelectFirstValidByGroup @ Select[validated, TrueQ[#["Valid"]] &];
  rejected = Select[validated, !TrueQ[#["Valid"]] &];

  LInt = Total[Lookup[valid, "Expression", {}]] // Expand // RelabelIndices;
  LBSM = (LFree + LInt) // Expand // RelabelIndices;
  LUV = (LSM + LBSM) // Expand // RelabelIndices;
  fullValidation = CheckAbort[Check[CheckLagrangian[LUV], $Failed], $Aborted];

  ingredientsPresent = T3IngredientsPresentQ[valid];

  <|
    "Status" -> If[TrueQ[fullValidation], "Success", "FullValidationFailed"],
    "Model" -> model,
    "LSM" -> LSM,
    "LFree" -> LFree,
    "LBSM" -> LBSM,
    "LUV" -> LUV,
    "AllowedInteractions" -> Lookup[valid, "Name", {}],
    "RejectedInteractions" -> Lookup[rejected, "Name", {}],
    "T3IngredientsPresent" -> ingredientsPresent,
    "WeinbergIngredientsPresent" -> ingredientsPresent,
    "FullValidation" -> fullValidation
  |>
];

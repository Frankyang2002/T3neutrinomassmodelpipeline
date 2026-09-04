(* T3 BSM field definitions, coupling definitions, and scalar accessors.
   Extracted from LagrangianBuilder.wl; physics behaviour is unchanged. *)

DefineT3Couplings[] := Module[{},
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

DefineT3Fields[model_Association] := Module[
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

  fidx = BSMIndexType[f["SU2"]];
  s1idx = BSMIndexType[s1["SU2"]];
  s2idx = BSMIndexType[s2["SU2"]];
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

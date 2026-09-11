(* T3 BSM field definitions, coupling definitions, and scalar accessors.
   The heavy/light assignment is now configurable so sequential threshold
   matching can integrate only the requested T3 fields at each stage. *)

If[!ValueQ[$T3HeavyFields],
  $T3HeavyFields = {"F", "S1", "S2"}
];

SetT3HeavyFields[fields_List] := Module[{allowed},
  allowed = {"F", "S1", "S2"};

  If[!SubsetQ[allowed, fields],
    Print["ERROR: unknown T3 heavy-field label(s): ", Complement[fields, allowed]];
    Return[$Failed]
  ];

  $T3HeavyFields = DeleteDuplicates[fields];
  True
];

T3MassSpec[label_String, mass_] := If[
  MemberQ[$T3HeavyFields, label],
  {Heavy, mass},
  {Light, mass}
];

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

DefineT3FieldObjects[
  model_Association,
  defineFlavor_: True,
  activeFields_List : {"F", "S1", "S2"}
] := Module[
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

  If[TrueQ[defineFlavor],
    DefineFlavorIndex[
      NFlavor,
      Lookup[f, "Multiplicity", 3],
      IndexAlphabet -> {"r", "s", "t"}
    ]
  ];

  findices = If[fidx === None, {NFlavor}, {fidx, NFlavor}];
  s1indices = If[s1idx === None, {}, {s1idx}];
  s2indices = If[s2idx === None, {}, {s2idx}];

  (* Integer-isospin SU(2) irreps (odd dimension) are real. A neutral
     fermion in such a representation may therefore be self-conjugate. *)
  selfConj = TrueQ[PossibleZeroQ[f["Y"]]] && OddQ[f["SU2"]];

  If[MemberQ[activeFields, "F"],
    If[selfConj,
      DefineField[
        NewFermion, Fermion,
        SelfConjugate -> True,
        Indices -> findices,
        Mass -> T3MassSpec["F", Lookup[f, "MassSymbol", MF]]
      ],
      DefineField[
        NewFermion, Fermion,
        SelfConjugate -> False,
        Indices -> findices,
        Charges -> {U1Y[f["Y"]]},
        Mass -> T3MassSpec["F", Lookup[f, "MassSymbol", MF]]
      ]
    ]
  ];

  If[MemberQ[activeFields, "S1"],
    DefineField[
      NewScalar1, Scalar,
      SelfConjugate -> False,
      Indices -> s1indices,
      Charges -> {U1Y[s1["Y"]]},
      Mass -> T3MassSpec["S1", Lookup[s1, "MassSymbol", MS1]]
    ]
  ];

  If[MemberQ[activeFields, "S2"],
    DefineField[
      NewScalar2, Scalar,
      SelfConjugate -> False,
      Indices -> s2indices,
      Charges -> {U1Y[s2["Y"]]},
      Mass -> T3MassSpec["S2", Lookup[s2, "MassSymbol", MS2]]
    ]
  ];

  True
];

DefineT3Fields[model_Association] := DefineT3FieldObjects[model, True];

(* Re-register only the three T3 fields with a new heavy/light hierarchy.
   The symbolic Field[...] objects already present in an intermediate EFT keep
   the same names, while Match consults the updated Matchete field registry. *)

(* Once a heavy field has been integrated out, its mass parameter remains
   inside Wilson coefficients such as 1/MF.  Removing the field can also
   remove Matchete's knowledge that its mass is a real/self-conjugate
   parameter.  Re-register those masses explicitly so Hermitian conjugation
   continues to identify \bar M = M in later EFT stages. *)
PreserveIntegratedT3MassCouplings[activeFields_List] := Module[{},
  If[!MemberQ[activeFields, "F"],
    Quiet @ Check[
      DefineCoupling[MF, SelfConjugate -> True],
      Null
    ]
  ];

  If[!MemberQ[activeFields, "S1"],
    Quiet @ Check[
      DefineCoupling[MS1, SelfConjugate -> True],
      Null
    ]
  ];

  If[!MemberQ[activeFields, "S2"],
    Quiet @ Check[
      DefineCoupling[MS2, SelfConjugate -> True],
      Null
    ]
  ];

  True
];

ReclassifyT3Fields[
  model_Association,
  heavyFields_List,
  activeFields_List : {"F", "S1", "S2"}
] := Module[{status},

  status = SetT3HeavyFields[heavyFields];
  If[status === $Failed, Return[$Failed]];

  (* Remove the old registry entries first.  Crucially, fields already
     integrated out at an earlier threshold must NOT be registered again.
     CheckLagrangian tests canonical normalization against the currently
     registered field content.  Re-registering F after F has been removed
     makes Matchete expect an F kinetic term that is no longer present. *)
  Quiet @ Check[RemoveField[NewFermion], Null];
  Quiet @ Check[RemoveField[NewScalar1], Null];
  Quiet @ Check[RemoveField[NewScalar2], Null];

  PreserveIntegratedT3MassCouplings[activeFields];

  DefineT3FieldObjects[model, False, activeFields]
];

ScalarFieldValue[1, False, ___] := NewScalar1[];
ScalarFieldValue[1, True,  ___] := Bar[NewScalar1[]];
ScalarFieldValue[2, False, ___] := NewScalar2[];
ScalarFieldValue[2, True,  ___] := Bar[NewScalar2[]];

ScalarFieldValue[1, False, i_] := NewScalar1[i];
ScalarFieldValue[1, True,  i_] := Bar[NewScalar1[i]];
ScalarFieldValue[2, False, i_] := NewScalar2[i];
ScalarFieldValue[2, True,  i_] := Bar[NewScalar2[i]];

(*
  T3RGETensorExport.wl

  Machine-readable bridge from the generalized T3 model construction to the
  Python RGE tensor adapter.

  This file is deliberately separate from LagrangianBuilder.wl.  It relies on
  the generalized invariant-basis helpers already defined by the builder and
  does not modify the production Lagrangian construction.

  Exported in this stage:
    - scalar multiplet metadata;
    - complete BSM scalar-quartic component terms in the complex basis;
    - the T3 H H S1 S2^\dagger mixing quartic and its hermitian conjugate;
    - raw Yukawa invariant components, with conjugation/orientation metadata.

  The raw Yukawa entries are intentionally not assigned global Weyl indices
  here.  The translation of Matchete fermion conventions into the y_ija
  convention of Eq. (4.85) must be fixed once and regression-tested against
  the legacy scotogenic calculation.
*)

ClearAll[
  T3RGEExactString,
  T3RGEFactor,
  T3RGEFullComponentTuple,
  T3RGEComponentEntries,
  T3RGEIndexedCouplingNames,
  T3RGEInvariantFamilyMetadata,
  T3RGEQuarticEntriesFromBasis,
  T3RGEHermitianQuarticEntriesFromBasis,
  T3RGEHiggsSelfEntries,
  T3RGEScalarSelfEntries,
  T3RGEHiggsPortalEntries,
  T3RGECrossScalarEntries,
  T3RGEMixingEntries,
  T3RGERawYukawaEntries,
  T3RGEExportAssociation,
  ExportT3RGETensors
];


(* ------------------------------------------------------------------------- *)
(* Exact JSON-safe representations                                           *)
(* ------------------------------------------------------------------------- *)

T3RGEExactString[value_] := ToString[
  InputForm[FullSimplify[value]],
  CharacterEncoding -> "ASCII"
];

T3RGEFactor[name_String, component_Integer, conjugated_] := <|
  "scalar_name" -> name,
  "component" -> component,
  "conjugated" -> TrueQ[conjugated]
|>;


(* Expected physical SU(2) invariant multiplicities for the T3 scalar
   potential.  SU(2) tensor products are multiplicity-free, while identical
   bosons restrict the self-quartic to the symmetric pair channels. *)
T3RGEIndexedCouplingNames[base_String, count_Integer?NonNegative] :=
  Table[base <> "Inv" <> ToString[index], {index, count}];

T3RGEInvariantFamilyMetadata[model_Association] := Module[
  {d1, d2, scalarSelf1, scalarSelf2, portal1, portal2, cross},

  d1 = model["Scalar1", "SU2"];
  d2 = model["Scalar2", "SU2"];
  scalarSelf1 = Ceiling[d1/2];
  scalarSelf2 = Ceiling[d2/2];
  portal1 = Min[2, d1];
  portal2 = Min[2, d2];
  cross = Min[d1, d2];

  {
    <|"family" -> "HiggsSelf", "kind" -> "quartic",
      "expected_physical_count" -> 1, "couplings" -> {"lambdaH"}|>,
    <|"family" -> "Scalar1Self", "kind" -> "quartic",
      "expected_physical_count" -> scalarSelf1,
      "couplings" -> T3RGEIndexedCouplingNames["lambdaS1", scalarSelf1]|>,
    <|"family" -> "Scalar2Self", "kind" -> "quartic",
      "expected_physical_count" -> scalarSelf2,
      "couplings" -> T3RGEIndexedCouplingNames["lambdaS2", scalarSelf2]|>,
    <|"family" -> "HiggsPortal1", "kind" -> "quartic",
      "expected_physical_count" -> portal1,
      "couplings" -> T3RGEIndexedCouplingNames["lambdaH1", portal1]|>,
    <|"family" -> "HiggsPortal2", "kind" -> "quartic",
      "expected_physical_count" -> portal2,
      "couplings" -> T3RGEIndexedCouplingNames["lambdaH2", portal2]|>,
    <|"family" -> "ScalarCross", "kind" -> "quartic",
      "expected_physical_count" -> cross,
      "couplings" -> T3RGEIndexedCouplingNames["lambda12", cross]|>,
    <|"family" -> "T3Mixing", "kind" -> "quartic",
      "expected_physical_count" -> 1, "couplings" -> {"lambdaT3"}|>,
    <|"family" -> "Yukawa1", "kind" -> "yukawa",
      "expected_physical_count" -> 1, "couplings" -> {"y1"}|>,
    <|"family" -> "Yukawa2", "kind" -> "yukawa",
      "expected_physical_count" -> 1, "couplings" -> {"y2"}|>
  }
];


(* ------------------------------------------------------------------------- *)
(* Tensor component enumeration                                              *)
(* ------------------------------------------------------------------------- *)

(* ArrayRules appends a default rule such as {_,_,_}->0.  Keep only rules
   whose left-hand side is an explicit integer tuple. *)
T3RGEComponentEntries[tensor_] := Module[{rules},
  If[!ArrayQ[tensor],
    Return[{{} -> tensor}]
  ];

  rules = ArrayRules[Normal[tensor]];

  Cases[
    rules,
    (indices_List -> value_) /;
      VectorQ[indices, IntegerQ] && value =!= 0 :>
        (indices -> value)
  ]
];


(* Reinsert the component value 1 for representation positions that were
   dropped because the field is an SU(2) singlet. *)
T3RGEFullComponentTuple[
  keptComponents_List,
  keep_List,
  numberOfFields_Integer
] := Module[{full, position},
  full = ConstantArray[1, numberOfFields];

  Do[
    position = keep[[slot]];
    full[[position]] = keptComponents[[slot]],
    {slot, Length[keep]}
  ];

  full
];


(* ------------------------------------------------------------------------- *)
(* Generic quartic component export                                          *)
(* ------------------------------------------------------------------------- *)

T3RGEQuarticEntriesFromBasis[
  basisTensor_,
  keep_List,
  factorSpecifications_List,
  coefficient_
] := Module[
  {rules, numberOfFields, fullComponents},

  numberOfFields = Length[factorSpecifications];
  rules = T3RGEComponentEntries[basisTensor];

  Map[
    Function[rule,
      fullComponents = T3RGEFullComponentTuple[
        First[rule],
        keep,
        numberOfFields
      ];

      <|
        "coefficient" -> T3RGEExactString[coefficient Last[rule]],
        "factors" -> MapThread[
          T3RGEFactor[
            #1[[1]],
            #2,
            #1[[2]]
          ] &,
          {factorSpecifications, fullComponents}
        ]
      |>
    ],
    rules
  ]
];


T3RGEHermitianQuarticEntriesFromBasis[
  basisTensor_,
  keep_List,
  factorSpecifications_List,
  coefficient_
] := Module[{forward, backwardSpecifications, backward},
  forward = T3RGEQuarticEntriesFromBasis[
    basisTensor,
    keep,
    factorSpecifications,
    coefficient
  ];

  backwardSpecifications = {
    #[[1]],
    !TrueQ[#[[2]]]
  } & /@ factorSpecifications;

  backward = T3RGEQuarticEntriesFromBasis[
    Conjugate[basisTensor],
    keep,
    backwardSpecifications,
    Conjugate[coefficient]
  ];

  Join[forward, backward]
];


(* ------------------------------------------------------------------------- *)
(* Scalar-potential sectors                                                  *)
(* ------------------------------------------------------------------------- *)

T3RGEHiggsSelfEntries[] := Module[
  {data, basis, keep, factors},

  data = PhysicalSU2InvariantBasis[
    {2, 2, 2, 2},
    {True, False, True, False},
    {{1, 3}, {2, 4}},
    {1, 2, 3, 4},
    False
  ];

  If[data === $Failed, Return[{}]];

  basis = data["Basis"];
  keep = data["Keep"];
  factors = {
    {"H", True},
    {"H", False},
    {"H", True},
    {"H", False}
  };

  If[Length[basis] =!= 1,
    Print[
      "ERROR: expected one physical Higgs self-quartic invariant, found ",
      Length[basis], "."
    ];
    Return[{}]
  ];

  (* Project convention: V contains lambdaH/2 (H^dagger H)^2. *)
  T3RGEQuarticEntriesFromBasis[
    First[basis],
    keep,
    factors,
    Symbol["lambdaH"]/2
  ]
];

T3RGEScalarSelfEntries[
  which_Integer,
  d_Integer?Positive,
  legacyQ_: False
] := Module[
  {data, basis, keep, couplingBase, coupling, factors},

  data = PhysicalSU2InvariantBasis[
    {d, d, d, d},
    {True, False, True, False},
    {{1, 3}, {2, 4}},
    {},
    legacyQ
  ];

  If[data === $Failed, Return[{}]];

  basis = data["Basis"];
  keep = data["Keep"];
  couplingBase = "lambdaS" <> ToString[which];

  factors = {
    {"S" <> ToString[which], True},
    {"S" <> ToString[which], False},
    {"S" <> ToString[which], True},
    {"S" <> ToString[which], False}
  };

  Flatten[
    Table[
      coupling = Symbol[couplingBase <> "Inv" <> ToString[invariant]];

      T3RGEQuarticEntriesFromBasis[
        basis[[invariant]],
        keep,
        factors,
        coupling/2
      ],
      {invariant, Length[basis]}
    ],
    1
  ]
];


T3RGEHiggsPortalEntries[
  which_Integer,
  d_Integer?Positive,
  legacyQ_: False
] := Module[
  {data, basis, keep, couplingBase, coupling, factors},

  data = PhysicalSU2InvariantBasis[
    {2, 2, d, d},
    {True, False, True, False},
    {},
    {1, 2},
    legacyQ
  ];

  If[data === $Failed, Return[{}]];

  basis = data["Basis"];
  keep = data["Keep"];
  couplingBase = "lambdaH" <> ToString[which];

  factors = {
    {"H", True},
    {"H", False},
    {"S" <> ToString[which], True},
    {"S" <> ToString[which], False}
  };

  Flatten[
    Table[
      coupling = Symbol[couplingBase <> "Inv" <> ToString[invariant]];

      (* Production Lagrangian uses 1/2 PlusHc[...] for this sector. *)
      T3RGEHermitianQuarticEntriesFromBasis[
        basis[[invariant]],
        keep,
        factors,
        coupling/2
      ],
      {invariant, Length[basis]}
    ],
    1
  ]
];


T3RGECrossScalarEntries[
  d1_Integer?Positive,
  d2_Integer?Positive,
  legacyQ_: False
] := Module[
  {data, basis, keep, coupling, factors},

  data = PhysicalSU2InvariantBasis[
    {d1, d1, d2, d2},
    {True, False, True, False},
    {},
    {},
    legacyQ
  ];

  If[data === $Failed, Return[{}]];

  basis = data["Basis"];
  keep = data["Keep"];

  factors = {
    {"S1", True},
    {"S1", False},
    {"S2", True},
    {"S2", False}
  };

  Flatten[
    Table[
      coupling = Symbol["lambda12Inv" <> ToString[invariant]];

      T3RGEQuarticEntriesFromBasis[
        basis[[invariant]],
        keep,
        factors,
        coupling
      ],
      {invariant, Length[basis]}
    ],
    1
  ]
];


(* The topology-defining H H S1 S2^\dagger invariant has multiplicity one
   for valid T3 assignments, but use the same general invariant-basis logic
   rather than assuming an explicit epsilon/generator form. *)
T3RGEMixingEntries[
  d1_Integer?Positive,
  d2_Integer?Positive,
  legacyQ_: False
] := Module[
  {data, basis, keep, factors},

  data = PhysicalSU2InvariantBasis[
    {2, 2, d1, d2},
    {False, False, False, True},
    {{1, 2}},
    {1, 2},
    legacyQ
  ];

  If[data === $Failed, Return[{}]];

  basis = data["Basis"];
  keep = data["Keep"];

  If[Length[basis] =!= 1,
    Print[
      "ERROR: T3 H H S1 S2^dagger sector has ", Length[basis],
      " physical invariants; exactly one is required by the current topology."
    ];
    Return[{}]
  ];

  factors = {
    {"H", False},
    {"H", False},
    {"S1", False},
    {"S2", True}
  };

  Flatten[
    T3RGEHermitianQuarticEntriesFromBasis[
      #,
      keep,
      factors,
      lambdaT3
    ] & /@ basis,
    1
  ]
];


(* ------------------------------------------------------------------------- *)
(* Raw Yukawa invariant export                                               *)
(* ------------------------------------------------------------------------- *)

(* Export the invariant tensor and field orientation without choosing the
   global left-handed Weyl basis.  This keeps the group-theory result exact
   while making the remaining convention-sensitive step explicit. *)
T3RGERawYukawaEntries[
  model_Association,
  which_Integer,
  legacyQ_: False
] := Module[
  {dF, dS, data, basis, keep, tensor, rules, fullComponents,
   scalarConjugated, fermionConjugated, couplingName},

  dF = model["Fermion", "SU2"];
  dS = model[
    If[which === 1, "Scalar1", "Scalar2"],
    "SU2"
  ];

  If[which === 1,
    data = PhysicalSU2InvariantBasis[
      {2, dF, dS},
      {True, True, False},
      {},
      {1},
      legacyQ
    ];
    scalarConjugated = False;
    fermionConjugated = True;
    couplingName = "y1",
    data = PhysicalSU2InvariantBasis[
      {2, dF, dS},
      {True, False, True},
      {},
      {1},
      legacyQ
    ];
    scalarConjugated = True;
    fermionConjugated = False;
    couplingName = "y2"
  ];

  If[data === $Failed, Return[{}]];

  basis = data["Basis"];
  keep = data["Keep"];

  If[Length[basis] =!= 1,
    Print[
      "ERROR: Yukawa", which,
      " invariant multiplicity is ", Length[basis],
      "; exactly one is required by the current T3 topology."
    ];
    Return[{}]
  ];

  Flatten[
    Table[
      tensor = basis[[invariant]];
      rules = T3RGEComponentEntries[tensor];

      Map[
        Function[rule,
          fullComponents = T3RGEFullComponentTuple[
            First[rule],
            keep,
            3
          ];

          <|
            "interaction" -> "Yukawa" <> ToString[which],
            "invariant" -> invariant,
            "coupling" -> couplingName,
            "cg_coefficient" -> T3RGEExactString[Last[rule]],
            "lepton" -> <|
              "name" -> "L",
              "component" -> fullComponents[[1]],
              "conjugated_representation" -> True
            |>,
            "fermion" -> <|
              "name" -> "F",
              "component" -> fullComponents[[2]],
              "conjugated_representation" -> fermionConjugated
            |>,
            "scalar" -> <|
              "name" -> "S" <> ToString[which],
              "component" -> fullComponents[[3]],
              "conjugated" -> scalarConjugated
            |>
          |>
        ],
        rules
      ],
      {invariant, Length[basis]}
    ],
    1
  ]
];


(* ------------------------------------------------------------------------- *)
(* Complete exchange association                                             *)
(* ------------------------------------------------------------------------- *)

T3RGEExportAssociation[model_Association] := Module[
  {d1, d2, dF, y1, y2, yF, legacyQ, quartics, mixingQuartics,
   rawYukawa1, rawYukawa2, rawYukawas},

  d1 = model["Scalar1", "SU2"];
  d2 = model["Scalar2", "SU2"];
  dF = model["Fermion", "SU2"];

  y1 = model["Scalar1", "Y"];
  y2 = model["Scalar2", "Y"];
  yF = model["Fermion", "Y"];

  legacyQ = T3LegacyModelQ[model];

  mixingQuartics = T3RGEMixingEntries[d1, d2, legacyQ];
  rawYukawa1 = T3RGERawYukawaEntries[model, 1, legacyQ];
  rawYukawa2 = T3RGERawYukawaEntries[model, 2, legacyQ];

  If[mixingQuartics === {} || rawYukawa1 === {} || rawYukawa2 === {},
    Print[
      "ERROR: topology invariant export failed; refusing to write a partial ",
      "T3 tensor exchange."
    ];
    Return[$Failed]
  ];

  quartics = Join[
    T3RGEHiggsSelfEntries[],
    T3RGEScalarSelfEntries[1, d1, legacyQ],
    T3RGEScalarSelfEntries[2, d2, legacyQ],
    T3RGEHiggsPortalEntries[1, d1, legacyQ],
    T3RGEHiggsPortalEntries[2, d2, legacyQ],
    T3RGECrossScalarEntries[d1, d2, legacyQ],
    mixingQuartics
  ];

  rawYukawas = Join[rawYukawa1, rawYukawa2];

  <|
    "schema_version" -> 1,
    "model_name" -> Lookup[model, "Class", "T3-general"],
    "scalar_basis_convention" -> "(R1,I1,R2,I2,...)",
    "complex_scalar_convention" -> "Phi=(R+i I)/Sqrt[2]",
    "quartic_tensor_convention" ->
      "lambda_abcd = d^4 V4/(dphi_a dphi_b dphi_c dphi_d)",
    "invariant_metadata_version" -> 1,
    "invariant_families" -> T3RGEInvariantFamilyMetadata[model],
    "scalars" -> {
      <|
        "name" -> "H",
        "su2_dimension" -> 2,
        "hypercharge" -> "1/2"
      |>,
      <|
        "name" -> "S1",
        "su2_dimension" -> d1,
        "hypercharge" -> T3RGEExactString[y1]
      |>,
      <|
        "name" -> "S2",
        "su2_dimension" -> d2,
        "hypercharge" -> T3RGEExactString[y2]
      |>
    },
    "raw_fermions" -> {
      <|
        "name" -> "L",
        "su2_dimension" -> 2,
        "hypercharge" -> "-1/2",
        "note" -> "Map to global LH Weyl basis in Python adapter."
      |>,
      <|
        "name" -> "F",
        "su2_dimension" -> dF,
        "hypercharge" -> T3RGEExactString[yF],
        "multiplicity" -> Lookup[model["Fermion"], "Multiplicity", 3],
        "self_conjugate" -> (
          TrueQ[PossibleZeroQ[yF]] && OddQ[dF]
        ),
        "note" -> "Conjugation depends on Yukawa orientation; see raw_yukawa_components."
      |>
    },
    "quartic_components" -> quartics,
    "raw_yukawa_components" -> rawYukawas,
    "wilson_components" -> {},
    "status" -> <|
      "quartic_component_export" -> "Ready",
      "raw_yukawa_export" -> "Ready",
      "weyl_yukawa_mapping" -> "ConstructInPython",
      "matched_wilson_export" -> "ConstructInPython"
    |>
  |>
];


ExportT3RGETensors[
  model_Association,
  outputPath_String
] := Module[{data, result},
  data = T3RGEExportAssociation[model];

  If[data === $Failed, Return[$Failed]];

  result = Quiet@Check[
    Export[outputPath, data, "RawJSON"],
    $Failed
  ];

  If[result === $Failed,
    Print["ERROR: could not export T3 RGE tensor exchange to ", outputPath];
    Return[$Failed]
  ];

  Print[
    "Exported T3 RGE tensor exchange: ",
    Length[data["quartic_components"]], " quartic component terms, ",
    Length[data["raw_yukawa_components"]], " raw Yukawa component terms."
  ];

  outputPath
];

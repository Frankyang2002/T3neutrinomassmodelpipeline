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
Get[FileNameJoin[{DirectoryName[$InputFileName], "T3RGEComponentExport.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "T3RGEPotentialExport.wl"}]];


T3RGERawYukawaEntries[
  model_Association,
  which_Integer
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
      {1}
    ];
    scalarConjugated = False;
    fermionConjugated = True;
    couplingName = "y1",
    data = PhysicalSU2InvariantBasis[
      {2, dF, dS},
      {True, False, True},
      {},
      {1}
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
  {d1, d2, dF, y1, y2, yF, quartics, mixingQuartics,
   rawYukawa1, rawYukawa2, rawYukawas},

  d1 = model["Scalar1", "SU2"];
  d2 = model["Scalar2", "SU2"];
  dF = model["Fermion", "SU2"];

  y1 = model["Scalar1", "Y"];
  y2 = model["Scalar2", "Y"];
  yF = model["Fermion", "Y"];
  mixingQuartics = T3RGEMixingEntries[d1, d2];
  rawYukawa1 = T3RGERawYukawaEntries[model, 1];
  rawYukawa2 = T3RGERawYukawaEntries[model, 2];

  If[mixingQuartics === {} || rawYukawa1 === {} || rawYukawa2 === {},
    Print[
      "ERROR: topology invariant export failed; refusing to write a partial ",
      "T3 tensor exchange."
    ];
    Return[$Failed]
  ];

  quartics = Join[
    T3RGEHiggsSelfEntries[],
    T3RGEScalarSelfEntries[1, d1],
    T3RGEScalarSelfEntries[2, d2],
    T3RGEHiggsPortalEntries[1, d1],
    T3RGEHiggsPortalEntries[2, d2],
    T3RGECrossScalarEntries[d1, d2],
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
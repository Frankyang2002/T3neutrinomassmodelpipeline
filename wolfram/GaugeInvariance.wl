(*
  GaugeInvariance.wl
  Resolves symbolic interaction templates against model field data and checks
  U(1) hypercharge conservation and whether the SU(2) tensor product contains
  a singlet.
*)
ClearAll[
  FieldHypercharge,
  InferFieldHypercharge,
  InferFieldSU2Representations,
  U1InvariantQ,
  SU2Isospin,
  SU2Combine,
  SU2ProductRepresentations,
  SU2SingletQ,
  ResolveTemplateFields,
  TemplateInvarianceReport,
  TemplateInvariantQ,
  ModelFieldData,
  T3InteractionTemplates,
  AllowedInteractionReports
];
(* ============================================================ *)
(* U(1) hypercharge                                              *)
(* ============================================================ *)
FieldHypercharge[field_Association] := Module[
  {hypercharge, conjugated},
  hypercharge = Lookup[field, "Y", Missing["Hypercharge"]];
  conjugated = TrueQ[Lookup[field, "Conjugated", False]];
  If[MissingQ[hypercharge],
    Return[$Failed]
  ];
  If[conjugated, -hypercharge, hypercharge]
];
U1InvariantQ[fields_List] := Module[
  {charges},
  charges = FieldHypercharge /@ fields;
  If[MemberQ[charges, $Failed],
    Return[False]
  ];
  TrueQ[Simplify[Total[charges] == 0]]
];

(*
  Determine the hypercharge of one field directly from an interaction
  template.  Each ordinary field contributes +Y and each conjugated field
  contributes -Y.  The target field is chosen so that the complete sum is 0.
*)
InferFieldHypercharge[
  template_Association,
  fieldData_Association,
  targetFieldName_String
] := Module[
  {
    templateName,
    references,
    targetPositions,
    targetPosition,
    targetReference,
    targetSign,
    otherReferences,
    otherCharges
  },

  templateName = Lookup[template, "Name", "<unnamed>"];
  references = Lookup[template, "Fields", Missing["Fields"]];

  If[MissingQ[references] || !ListQ[references],
    Print[
      "ERROR: Template ", templateName,
      " has no valid Fields list."
    ];
    Return[$Failed]
  ];

  targetPositions = Flatten@Position[
    Lookup[references, "Name", Missing["FieldName"]],
    targetFieldName
  ];

  If[Length[targetPositions] =!= 1,
    Print[
      "ERROR: Template ", templateName,
      " must contain exactly one occurrence of target field ",
      targetFieldName, ". Found ", Length[targetPositions], "."
    ];
    Return[$Failed]
  ];

  targetPosition = First[targetPositions];
  targetReference = references[[targetPosition]];
  targetSign = If[
    TrueQ[Lookup[targetReference, "Conjugated", False]],
    -1,
    1
  ];

  otherReferences = Delete[references, targetPosition];
  otherCharges = Map[
    Function[fieldReference,
      Module[{fieldName, baseField, conjugated},
        fieldName = Lookup[
          fieldReference,
          "Name",
          Missing["FieldName"]
        ];

        If[MissingQ[fieldName] || !KeyExistsQ[fieldData, fieldName],
          Print[
            "ERROR: Unknown or missing field in template ",
            templateName, ": ", fieldName
          ];
          Return[$Failed]
        ];

        baseField = fieldData[fieldName];
        conjugated = TrueQ[
          Lookup[fieldReference, "Conjugated", False]
        ];

        FieldHypercharge[
          Join[
            baseField,
            <|"Conjugated" -> conjugated|>
          ]
        ]
      ]
    ],
    otherReferences
  ];

  If[MemberQ[otherCharges, $Failed],
    Print[
      "ERROR: Could not read all non-target hypercharges in template ",
      templateName, "."
    ];
    Return[$Failed]
  ];

  Simplify[-Total[otherCharges]/targetSign]
];

(* ============================================================ *)
(* SU(2) representation products                                 *)
(* ============================================================ *)
SU2Isospin[dimension_Integer?Positive] :=
  (dimension - 1)/2;
SU2Combine[
  rep1_Integer?Positive,
  rep2_Integer?Positive
] := Module[
  {j1, j2, resultingIsospins},
  j1 = SU2Isospin[rep1];
  j2 = SU2Isospin[rep2];
  resultingIsospins =
    Range[
      Abs[j1 - j2],
      j1 + j2,
      1
    ];
  2 resultingIsospins + 1
];
SU2ProductRepresentations[{}] := {};
SU2ProductRepresentations[
  {representation_Integer?Positive}
] := {representation};
SU2ProductRepresentations[
  representations_List
] /; Length[representations] >= 2 := Module[
  {previousProducts, finalRepresentation},
  previousProducts =
    SU2ProductRepresentations[Most[representations]];
  finalRepresentation = Last[representations];
  DeleteDuplicates[
    Flatten[
      SU2Combine[#, finalRepresentation] & /@ previousProducts
    ]
  ]
];
SU2SingletQ[representations_List] := Module[
  {products},
  If[representations === {},
    Return[False]
  ];
  If[
    !AllTrue[
      representations,
      IntegerQ[#] && Positive[#] &
    ],
    Return[False]
  ];
  products = SU2ProductRepresentations[representations];
  MemberQ[products, 1]
];
(*
  Infer all SU(2) representation dimensions for one target field that can
  make an interaction template contain a singlet.

  For SU(2), conjugation does not change the irrep dimension.  If the tensor
  product of all non-target fields contains irreps {r1,r2,...}, then the
  target may transform as any of those same irreps, because j \otimes j
  always contains J=0.

  Example: for L x Scalar x Fermion with L a doublet and Scalar dimension d,
  the allowed Fermion dimensions are precisely those in 2 x d.  Thus
      d=1 -> {2}
      d=2 -> {1,3}
      d=3 -> {2,4}.
*)
InferFieldSU2Representations[
  template_Association,
  fieldData_Association,
  targetFieldName_String
] := Module[
  {
    templateName,
    references,
    targetPositions,
    targetPosition,
    otherReferences,
    otherRepresentations
  },

  templateName = Lookup[template, "Name", "<unnamed>"];
  references = Lookup[template, "Fields", Missing["Fields"]];

  If[MissingQ[references] || !ListQ[references],
    Print[
      "ERROR: Template ", templateName,
      " has no valid Fields list."
    ];
    Return[$Failed]
  ];

  targetPositions = Flatten@Position[
    Lookup[references, "Name", Missing["FieldName"]],
    targetFieldName
  ];

  If[Length[targetPositions] =!= 1,
    Print[
      "ERROR: Template ", templateName,
      " must contain exactly one occurrence of target field ",
      targetFieldName, ". Found ", Length[targetPositions], "."
    ];
    Return[$Failed]
  ];

  targetPosition = First[targetPositions];
  otherReferences = Delete[references, targetPosition];

  otherRepresentations = Map[
    Function[fieldReference,
      Module[{fieldName, baseField, representation},
        fieldName = Lookup[
          fieldReference,
          "Name",
          Missing["FieldName"]
        ];

        If[MissingQ[fieldName] || !KeyExistsQ[fieldData, fieldName],
          Print[
            "ERROR: Unknown or missing field in template ",
            templateName, ": ", fieldName
          ];
          Return[$Failed]
        ];

        baseField = fieldData[fieldName];
        representation = Lookup[
          baseField,
          "SU2",
          Missing["SU2Representation"]
        ];

        If[
          MissingQ[representation] ||
          !IntegerQ[representation] ||
          representation < 1,
          Print[
            "ERROR: Invalid SU(2) representation for field ",
            fieldName, " in template ", templateName, "."
          ];
          Return[$Failed]
        ];

        representation
      ]
    ],
    otherReferences
  ];

  If[MemberQ[otherRepresentations, $Failed],
    Return[$Failed]
  ];

  Sort@DeleteDuplicates@SU2ProductRepresentations[otherRepresentations]
];

(* ============================================================ *)
(* Template resolution                                          *)
(* ============================================================ *)
ResolveTemplateFields[
  template_Association,
  fieldData_Association
] := Module[
  {templateName, references},
  templateName = Lookup[template, "Name", "<unnamed>"];
  references = Lookup[template, "Fields", Missing["Fields"]];
  If[MissingQ[references] || !ListQ[references],
    Print[
      "Template ",
      templateName,
      " has no valid Fields list."
    ];
    Return[$Failed]
  ];
  Map[
    Function[fieldReference,
      Module[
        {fieldName, baseField},
        fieldName =
          Lookup[
            fieldReference,
            "Name",
            Missing["FieldName"]
          ];
        If[MissingQ[fieldName],
          Print[
            "Template ",
            templateName,
            " contains a field without a Name."
          ];
          Return[$Failed]
        ];
        If[!KeyExistsQ[fieldData, fieldName],
          Print[
            "Unknown field ",
            fieldName,
            " in template ",
            templateName,
            "."
          ];
          Return[$Failed]
        ];
        baseField = fieldData[fieldName];
        Join[
          baseField,
          <|
            "Name" -> fieldName,
            "Conjugated" ->
              TrueQ[
                Lookup[
                  fieldReference,
                  "Conjugated",
                  False
                ]
              ]
          |>
        ]
      ]
    ],
    references
  ]
];
(* ============================================================ *)
(* Diagnostic report                                            *)
(* ============================================================ *)
TemplateInvarianceReport[
  template_Association,
  fieldData_Association
] := Module[
  {
    templateName,
    resolvedFields,
    su2Representations,
    hypercharges,
    hyperchargeSum,
    su2Products,
    u1Invariant,
    su2Invariant
  },
  templateName = Lookup[template, "Name", "<unnamed>"];
  resolvedFields =
    ResolveTemplateFields[
      template,
      fieldData
    ];
  If[
    resolvedFields === $Failed ||
    MemberQ[resolvedFields, $Failed],
    Return[
      <|
        "Name" -> templateName,
        "ResolvedFields" -> {},
        "HyperchargeSum" -> Missing["NotAvailable"],
        "SU2ProductRepresentations" -> {},
        "U1Invariant" -> False,
        "SU2Invariant" -> False,
        "Invariant" -> False
      |>
    ]
  ];
  hypercharges = FieldHypercharge /@ resolvedFields;
  hyperchargeSum =
    If[
      MemberQ[hypercharges, $Failed],
      Missing["NotAvailable"],
      Simplify[Total[hypercharges]]
    ];
  su2Representations =
    Lookup[
      resolvedFields,
      "SU2",
      Missing["SU2Representation"]
    ];
  su2Products =
    If[
      AnyTrue[su2Representations, MissingQ],
      {},
      SU2ProductRepresentations[su2Representations]
    ];
  u1Invariant =
    !MissingQ[hyperchargeSum] &&
    TrueQ[Simplify[hyperchargeSum == 0]];
  su2Invariant = MemberQ[su2Products, 1];
  <|
    "Name" -> templateName,
    "ResolvedFields" -> resolvedFields,
    "HyperchargeSum" -> hyperchargeSum,
    "SU2ProductRepresentations" -> su2Products,
    "U1Invariant" -> u1Invariant,
    "SU2Invariant" -> su2Invariant,
    "Invariant" -> (u1Invariant && su2Invariant)
  |>
];
TemplateInvariantQ[
  template_Association,
  fieldData_Association
] :=
  TrueQ[
    TemplateInvarianceReport[
      template,
      fieldData
    ]["Invariant"]
  ];
(* ============================================================ *)
(* Model-specific field data                                    *)
(* ============================================================ *)
ModelFieldData[model_Association] := <|
  "L" ->
    Lookup[
      model,
      "Lepton",
      Lookup[
        model,
        "L",
        <|
          "Type" -> "Fermion",
          "SU2" -> 2,
          "Y" -> -1/2
        |>
      ]
    ],
  "H" ->
    Lookup[
      model,
      "Higgs",
      Lookup[
        model,
        "H",
        <|
          "Type" -> "ComplexScalar",
          "SU2" -> 2,
          "Y" -> 1/2
        |>
      ]
    ],
  "Fermion" ->
    Lookup[
      model,
      "Fermion",
      Missing["Fermion"]
    ],
  "Scalar" ->
    Lookup[
      model,
      "Scalar",
      Missing["Scalar"]
    ]
|>;
(* ============================================================ *)
(* Candidate interactions                                       *)
(* ============================================================ *)
T3InteractionTemplates = {
  <|
    "Name" -> "Yukawa",
    "Fields" -> {
      <|"Name" -> "L", "Conjugated" -> True|>,
      <|"Name" -> "Scalar", "Conjugated" -> True|>,
      <|"Name" -> "Fermion", "Conjugated" -> False|>
    }
  |>,
  <|
    "Name" -> "ScalarSelf",
    "Fields" -> {
      <|"Name" -> "Scalar", "Conjugated" -> True|>,
      <|"Name" -> "Scalar", "Conjugated" -> False|>,
      <|"Name" -> "Scalar", "Conjugated" -> True|>,
      <|"Name" -> "Scalar", "Conjugated" -> False|>
    }
  |>,
  <|
    "Name" -> "Portal3",
    "Fields" -> {
      <|"Name" -> "H", "Conjugated" -> True|>,
      <|"Name" -> "H", "Conjugated" -> False|>,
      <|"Name" -> "Scalar", "Conjugated" -> True|>,
      <|"Name" -> "Scalar", "Conjugated" -> False|>
    }
  |>,
  <|
    "Name" -> "Portal4",
    "Fields" -> {
      <|"Name" -> "H", "Conjugated" -> True|>,
      <|"Name" -> "Scalar", "Conjugated" -> False|>,
      <|"Name" -> "Scalar", "Conjugated" -> True|>,
      <|"Name" -> "H", "Conjugated" -> False|>
    }
  |>,
  <|
    "Name" -> "Portal5",
    "Fields" -> {
      <|"Name" -> "H", "Conjugated" -> True|>,
      <|"Name" -> "Scalar", "Conjugated" -> False|>,
      <|"Name" -> "H", "Conjugated" -> True|>,
      <|"Name" -> "Scalar", "Conjugated" -> False|>
    }
  |>
};
AllowedInteractionReports[model_Association] := Module[
  {fieldData, reports},
  fieldData = ModelFieldData[model];
  If[
    MissingQ[fieldData["Fermion"]] ||
    MissingQ[fieldData["Scalar"]],
    Print["ERROR: Model must define Fermion and Scalar."];
    Return[$Failed]
  ];
  reports =
    TemplateInvarianceReport[
      #,
      fieldData
    ] & /@ T3InteractionTemplates;
  Select[
    reports,
    TrueQ[Lookup[#, "Invariant", False]] &
  ]
];

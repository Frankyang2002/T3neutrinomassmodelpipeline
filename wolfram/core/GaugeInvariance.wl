(*
  GaugeInvariance.wl

  Generic electroweak gauge-invariance utilities for symbolic interaction
  templates.  This file owns the shared U(1)_Y and SU(2) machinery used by
  the higher-level T3 invariance checks.

  Conventions
  -----------
  - "Y" stores hypercharge in Q = T3 + Y.
  - "Conjugated" reverses U(1) hypercharge but does not change an SU(2)
    irrep dimension.
  - "SU2" stores the representation dimension d = 2 j + 1.
  - Public predicates return strict True/False; malformed input fails safely.
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
  GaugeInvarianceReport,
  TemplateInvarianceReport,
  TemplateInvariantQ,
  ModelFieldData,
  T3InteractionTemplates,
  AllowedInteractionReports
];

(* ============================================================ *)
(* U(1)_Y                                                       *)
(* ============================================================ *)

FieldHypercharge[field_Association] := Module[
  {hypercharge},
  hypercharge = Lookup[field, "Y", Missing["Hypercharge"]];
  If[MissingQ[hypercharge], Return[$Failed]];

  If[
    TrueQ[Lookup[field, "Conjugated", False]],
    -hypercharge,
    hypercharge
  ]
];

U1InvariantQ[fields_List] := Module[
  {charges = FieldHypercharge /@ fields},
  If[MemberQ[charges, $Failed], Return[False]];
  TrueQ[Simplify[Total[charges] == 0]]
];

(*
  Infer the hypercharge of one target field from charge conservation.
  A conjugated occurrence contributes -Y, so the sign of the target in the
  template must be retained when solving Sum[Y_i] = 0.
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
    targetReference,
    targetSign,
    otherCharges
  },

  templateName = Lookup[template, "Name", "<unnamed>"];
  references = Lookup[template, "Fields", Missing["Fields"]];

  If[MissingQ[references] || !ListQ[references],
    Print["ERROR: Template ", templateName, " has no valid Fields list."];
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

  targetReference = references[[First[targetPositions]]];
  targetSign = If[TrueQ[Lookup[targetReference, "Conjugated", False]], -1, 1];

  otherCharges = Map[
    Function[fieldReference,
      Module[{fieldName, baseField},
        fieldName = Lookup[fieldReference, "Name", Missing["FieldName"]];

        If[MissingQ[fieldName] || !KeyExistsQ[fieldData, fieldName],
          Print[
            "ERROR: Unknown or missing field in template ",
            templateName, ": ", fieldName
          ];
          Return[$Failed]
        ];

        baseField = fieldData[fieldName];
        FieldHypercharge@Join[
          baseField,
          <|"Conjugated" -> TrueQ[Lookup[fieldReference, "Conjugated", False]]|>
        ]
      ]
    ],
    Delete[references, First[targetPositions]]
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
(* SU(2) representation products                                *)
(* ============================================================ *)

SU2Isospin[dimension_Integer?Positive] := (dimension - 1)/2;

SU2Combine[rep1_Integer?Positive, rep2_Integer?Positive] := Module[
  {j1 = SU2Isospin[rep1], j2 = SU2Isospin[rep2]},
  2 Range[Abs[j1 - j2], j1 + j2, 1] + 1
];

SU2ProductRepresentations[{}] := {};
SU2ProductRepresentations[{representation_Integer?Positive}] := {representation};
SU2ProductRepresentations[representations_List] /; Length[representations] >= 2 := Module[
  {previousProducts, finalRepresentation},
  previousProducts = SU2ProductRepresentations[Most[representations]];
  finalRepresentation = Last[representations];

  DeleteDuplicates@Flatten[
    SU2Combine[#, finalRepresentation] & /@ previousProducts
  ]
];

SU2SingletQ[representations_List] := Module[
  {},
  If[representations === {}, Return[False]];
  If[
    !AllTrue[representations, IntegerQ[#] && Positive[#] &],
    Return[False]
  ];
  MemberQ[SU2ProductRepresentations[representations], 1]
];

(*
  Infer all target irrep dimensions that allow the complete product to contain
  an SU(2) singlet.  SU(2) irreps are self-conjugate, so conjugation does not
  alter the dimension.  If the non-target fields combine to irrep r, choosing
  the target in the same irrep r permits r x r to contain the singlet.
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
    otherRepresentations
  },

  templateName = Lookup[template, "Name", "<unnamed>"];
  references = Lookup[template, "Fields", Missing["Fields"]];

  If[MissingQ[references] || !ListQ[references],
    Print["ERROR: Template ", templateName, " has no valid Fields list."];
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

  otherRepresentations = Map[
    Function[fieldReference,
      Module[{fieldName, representation},
        fieldName = Lookup[fieldReference, "Name", Missing["FieldName"]];

        If[MissingQ[fieldName] || !KeyExistsQ[fieldData, fieldName],
          Print[
            "ERROR: Unknown or missing field in template ",
            templateName, ": ", fieldName
          ];
          Return[$Failed]
        ];

        representation = Lookup[
          fieldData[fieldName],
          "SU2",
          Missing["SU2Representation"]
        ];

        If[
          MissingQ[representation] || !IntegerQ[representation] || representation < 1,
          Print[
            "ERROR: Invalid SU(2) representation for field ",
            fieldName, " in template ", templateName, "."
          ];
          Return[$Failed]
        ];

        representation
      ]
    ],
    Delete[references, First[targetPositions]]
  ];

  If[MemberQ[otherRepresentations, $Failed], Return[$Failed]];
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
    Print["Template ", templateName, " has no valid Fields list."];
    Return[$Failed]
  ];

  Map[
    Function[fieldReference,
      Module[{fieldName, baseField},
        fieldName = Lookup[fieldReference, "Name", Missing["FieldName"]];

        If[MissingQ[fieldName],
          Print["Template ", templateName, " contains a field without a Name."];
          Return[$Failed]
        ];

        If[!KeyExistsQ[fieldData, fieldName],
          Print["Unknown field ", fieldName, " in template ", templateName, "."];
          Return[$Failed]
        ];

        baseField = fieldData[fieldName];
        Join[
          baseField,
          <|
            "Name" -> fieldName,
            "Conjugated" -> TrueQ[Lookup[fieldReference, "Conjugated", False]]
          |>
        ]
      ]
    ],
    references
  ]
];

(* ============================================================ *)
(* Gauge diagnostic report                                      *)
(* ============================================================ *)

(*
  GaugeInvarianceReport is the canonical electroweak report.  Keeping it
  separate from TemplateInvarianceReport lets InvarianceChecker.wl extend the
  report with Z2 without duplicating any U(1) or SU(2) calculation.
*)
GaugeInvarianceReport[
  template_Association,
  fieldData_Association
] := Module[
  {
    templateName,
    resolvedFields,
    hypercharges,
    hyperchargeSum,
    su2Representations,
    su2Products,
    u1Invariant,
    su2Invariant
  },

  templateName = Lookup[template, "Name", "<unnamed>"];
  resolvedFields = ResolveTemplateFields[template, fieldData];

  If[resolvedFields === $Failed || MemberQ[resolvedFields, $Failed],
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
  hyperchargeSum = If[
    MemberQ[hypercharges, $Failed],
    Missing["NotAvailable"],
    Simplify[Total[hypercharges]]
  ];

  su2Representations = Lookup[
    resolvedFields,
    "SU2",
    Missing["SU2Representation"]
  ];
  su2Products = If[
    AnyTrue[su2Representations, MissingQ],
    {},
    SU2ProductRepresentations[su2Representations]
  ];

  u1Invariant = !MissingQ[hyperchargeSum] && TrueQ[Simplify[hyperchargeSum == 0]];
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

(* Backward-compatible gauge-only public API. *)
TemplateInvarianceReport[template_Association, fieldData_Association] :=
  GaugeInvarianceReport[template, fieldData];

TemplateInvariantQ[template_Association, fieldData_Association] :=
  TrueQ[TemplateInvarianceReport[template, fieldData]["Invariant"]];

(* ============================================================ *)
(* Existing T3 model helpers                                    *)
(* ============================================================ *)

ModelFieldData[model_Association] := <|
  "L" -> Lookup[
    model,
    "Lepton",
    Lookup[
      model,
      "L",
      <|"Type" -> "Fermion", "SU2" -> 2, "Y" -> -1/2|>
    ]
  ],
  "H" -> Lookup[
    model,
    "Higgs",
    Lookup[
      model,
      "H",
      <|"Type" -> "ComplexScalar", "SU2" -> 2, "Y" -> 1/2|>
    ]
  ],
  "Fermion" -> Lookup[model, "Fermion", Missing["Fermion"]],
  "Scalar" -> Lookup[model, "Scalar", Missing["Scalar"]]
|>;

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

  If[MissingQ[fieldData["Fermion"]] || MissingQ[fieldData["Scalar"]],
    Print["ERROR: Model must define Fermion and Scalar."];
    Return[$Failed]
  ];

  reports = TemplateInvarianceReport[#, fieldData] & /@ T3InteractionTemplates;
  Select[reports, TrueQ[Lookup[#, "Invariant", False]] &]
];

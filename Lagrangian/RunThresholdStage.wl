(* RunThresholdStage.wl

   Execute ONE sequential T3 threshold in a completely fresh Matchete kernel.

   CLI:
     result.wxf EFTOrder LoopOrder dS1 dS2 dF alpha
     active-fields-csv heavy-fields-csv input-tree.wxf input-loop.wxf
     transition-tree.wxf transition-full.wxf cg-registry.wxf
     [eft1-running-insertion.wl]

   The optional final argument is the full-flavor EFT1 leading-log insertion
   produced by EFT1WilsonFlavorMatcheteExporter.py.  It is consumed ONLY by
   fixed-order tree propagation [C], never by the new one-loop threshold
   matching [B].

   The previous EFT is imported only after the SM, representations, fields,
   couplings, and CG tensors have been freshly registered for THIS threshold.
*)

ClearAll["Global`*"];

scriptDirectory = DirectoryName @ ExpandFileName[$InputFileName];

args = Rest[$ScriptCommandLine];
If[Length[args] < 14,
  Print[
    "ERROR: RunThresholdStage.wl received too few arguments: ",
    Length[args],
    " (expected 14). Raw $ScriptCommandLine = ",
    InputForm[$ScriptCommandLine]
  ];
  Exit[2]
];

resultPath = ExpandFileName @ args[[1]];
eftOrder = ToExpression @ args[[2]];
loopOrder = ToExpression @ args[[3]];
dS1 = ToExpression @ args[[4]];
dS2 = ToExpression @ args[[5]];
dF = ToExpression @ args[[6]];
alphaToken = args[[7]];
alpha = Which[
  StringStartsQ[alphaToken, "m"],
    -ToExpression @ StringDrop[alphaToken, 1],
  StringStartsQ[alphaToken, "p"],
    ToExpression @ StringDrop[alphaToken, 1],
  True,
    ToExpression @ alphaToken
];
activeFields = Select[StringSplit[args[[8]], ","], StringLength[#] > 0 &];
heavyFields = Select[StringSplit[args[[9]], ","], StringLength[#] > 0 &];
treePath = ExpandFileName @ args[[10]];
loopPath = ExpandFileName @ args[[11]];
transitionTreePath = ExpandFileName @ args[[12]];
transitionFullPath = ExpandFileName @ args[[13]];
cgPath = ExpandFileName @ args[[14]];

runningInsertionPath = If[
  Length[args] >= 15 && StringLength[args[[15]]] > 0,
  ExpandFileName @ args[[15]],
  None
];

Print["[fresh kernel] Loading Matchete..."];
matcheteLoaded = UsingFrontEnd[Needs["Matchete`"]; True];
If[!TrueQ[matcheteLoaded],
  Print["ERROR: Matchete failed to load in fresh threshold kernel."];
  Exit[3]
];

modelCatalogFile =
  FileNameJoin[{scriptDirectory, "T3ModelCatalog.wl"}];
builderFile =
  FileNameJoin[{scriptDirectory, "LagrangianBuilder.wl"}];
matchingFile =
  FileNameJoin[{scriptDirectory, "RunMatching.wl"}];

Get[modelCatalogFile];
Get[builderFile];
Get[matchingFile];

model = T3ModelFromDimensions[dS1, dS2, dF, alpha];
If[model === $Failed,
  Print["ERROR: could not reconstruct T3 model in fresh threshold kernel."];
  Exit[4]
];

Print[
  "[fresh kernel] Active T3 fields: ", activeFields,
  "; Heavy now: ", heavyFields
];

(* Completely new registry. *)
ResetAll[];
LSM = LoadModel["SM"];
If[LSM === $Failed,
  Print["ERROR: SM failed to load in fresh threshold kernel."];
  Exit[5]
];

If[SetT3HeavyFields[heavyFields] === $Failed,
  Print["ERROR: invalid heavy-field group."];
  Exit[6]
];

If[DefineT3FieldObjects[model, True, activeFields] === $Failed,
  Print["ERROR: active T3 field definition failed."];
  Exit[7]
];

If[Check[DefineT3Couplings[], $Failed] === $Failed,
  Print["ERROR: T3 coupling definition failed."];
  Exit[8]
];

cgStatus = DefineT3InvariantCGs[model];
If[cgStatus === $Failed,
  Print["ERROR: T3 CG definition failed."];
  Exit[9]
];

(* The intermediate EFT contains not only the three topology CGs above,
   but also every scalar-potential invariant tensor that was generated in
   the original UV kernel (for example T3HiggsPortal1CGInv1/Inv2,
   T3Scalar1SelfCGInv*, T3CrossScalarCGInv*, ...).

   WXF serialises the symbolic CG names, not Matchete's private tensor
   registry.  Therefore every one of those named CGs must be reconstructed
   in this fresh kernel before the imported EFT can be interpreted.  Calling
   the scalar-potential candidate builder here is intentional: we use its
   deterministic side effect of defining the same invariant CG families and
   their couplings; the candidate expressions themselves are discarded. *)
Print["[fresh kernel] Reconstructing scalar-potential invariant registry..."];
scalarRegistryCandidates = BuildGeneralScalarPotentialCandidates[model];

If[!ListQ[scalarRegistryCandidates],
  Print["ERROR: scalar-potential invariant registry reconstruction failed."];
  Exit[91]
];

Print[
  "[fresh kernel] Scalar-potential registry reconstructed: ",
  Length[scalarRegistryCandidates],
  " invariant candidate(s)."
];

Print["[fresh kernel] Restoring serialized intermediate-EFT CG registry..."];
serializedCGRegistry = Quiet @ Check[Import[cgPath, "WXF"], $Failed];

If[serializedCGRegistry === $Failed || !ListQ[serializedCGRegistry],
  Print["ERROR: failed to import intermediate-EFT CG registry."];
  Exit[92]
];

restoredCGCount = 0;

Scan[
  Function[item,
    Module[{name, reps, tensor, existing},
      name = item["Name"];
      reps = item["Reps"];
      tensor = item["Tensor"];

      existing = Quiet @ Check[
        Matchete`PackageScope`$CGproperties[name, Indices],
        $Failed
      ];

      If[
        existing === $Failed || MissingQ[existing],
        If[
          Quiet @ Check[
            DefineCG[name, reps, tensor];
            True,
            False
          ],
          restoredCGCount++
        ]
      ]
    ]
  ],
  serializedCGRegistry
];

Print[
  "[fresh kernel] Restored ",
  restoredCGCount,
  " previously missing CG definition(s)."
];

PreserveIntegratedT3MassCouplings[activeFields];

(* EOMSimplify at the previous threshold introduces effective couplings for
   lower-power / quadratic terms.  In the T3 intermediate EFT the scalar
   quadratic coefficients CNewScalar12 and CNewScalar22 multiply
       Bar[NewScalar1] NewScalar1
   and
       Bar[NewScalar2] NewScalar2,
   respectively, so they are Hermitian (self-conjugate) couplings.

   Their symbolic names survive WXF serialization, but Matchete's private
   coupling metadata does not.  Re-register them in the fresh kernel before
   importing/matching the inherited one-loop EFT. *)
Print["[fresh kernel] Reconstructing inherited effective-coupling metadata..."];

If[MemberQ[activeFields, "S1"],
  Quiet @ Check[
    DefineCoupling[CNewScalar12, SelfConjugate -> True],
    Null
  ]
];

If[MemberQ[activeFields, "S2"],
  Quiet @ Check[
    DefineCoupling[CNewScalar22, SelfConjugate -> True],
    Null
  ]
];

Print[
  "[fresh kernel] Effective scalar couplings restored for: ",
  Intersection[activeFields, {"S1", "S2"}]
];

(* Optional full-flavor leading-log running generated between threshold 1
   and this threshold.

   This insertion is already O(hbar):
       Delta L_run = hbar Log[muLow/muHigh] beta^(1)[L^(0)].

   It MUST therefore enter only [C], where O(hbar) terms are propagated by
   tree matching.  Feeding it to [B] would one-loop match an already one-loop
   quantity and create an O(hbar^2) contribution. *)
eft1WilsonRunningHeavyInsertion = 0;
eft1WilsonRunningWeinbergCoefficient = 0;
eft1WilsonRunningWeinbergDefinition = {};
eft1WilsonRunningLogReplacement = {};
eft1WilsonRunningDirectWeinberg = 0;
eft1WilsonRunningDirectWeinbergAdded = False;
eft1WilsonRunningDirectWeinbergEqualScaleVanishes = Missing["NotChecked"];
eft1WilsonRunningLoaded = False;

If[runningInsertionPath =!= None,
  If[!FileExistsQ[runningInsertionPath],
    Print[
      "ERROR: requested EFT1 running insertion does not exist: ",
      runningInsertionPath
    ];
    Exit[93]
  ];

  Print[
    "[fresh kernel] Loading full-flavor EFT1 running insertion: ",
    runningInsertionPath
  ];

  (* Do not wrap Get[...] in Check here.  The generated insertion can cause
     harmless Matchete/Wolfram messages while defining indexed expressions;
     Check would treat any such message as a hard load failure even when the
     file was read successfully.  Load quietly, then validate the definitions
     we actually require. *)
  runningLoadResult = Quiet[Get[runningInsertionPath]];

  If[!ValueQ[EFT1WilsonRunningHeavyInsertion],
    Print[
      "ERROR: EFT1 running insertion was read, but ",
      "EFT1WilsonRunningHeavyInsertion is not defined."
    ];
    Print["  Get result: ", InputForm[runningLoadResult]];
    Exit[94]
  ];

  eft1WilsonRunningHeavyInsertion =
    Expand[EFT1WilsonRunningHeavyInsertion];

  If[!FreeQ[eft1WilsonRunningHeavyInsertion, Get | Import],
    Print[
      "ERROR: running insertion contains unresolved file-loading constructs."
    ];
    Exit[95]
  ];

  If[ValueQ[EFT1WilsonRunningLogReplacement],
    eft1WilsonRunningLogReplacement =
      EFT1WilsonRunningLogReplacement;
  ];

  (* The pure-SM Weinberg running term does not need threshold-2 matching.
     Keep its coefficient separate so the parent/final-EFT layer can carry it
     directly once the Weinberg operator is assembled there. *)
  If[DownValues[EFT1WilsonRunningWeinbergCoefficient] =!= {},
    eft1WilsonRunningWeinbergDefinition =
      DownValues[EFT1WilsonRunningWeinbergCoefficient];

    eft1WilsonRunningWeinbergCoefficient =
      EFT1WilsonRunningWeinbergCoefficient[
        Index[EFT1FlavorP, Flavor],
        Index[EFT1FlavorQ, Flavor]
      ];

    If[ListQ[eft1WilsonRunningLogReplacement] ||
       Head[eft1WilsonRunningLogReplacement] === Rule,
      eft1WilsonRunningWeinbergCoefficient =
        Expand[
          eft1WilsonRunningWeinbergCoefficient /.
            eft1WilsonRunningLogReplacement
        ];
    ];
  ];

  eft1WilsonRunningLoaded = True;

  Print[
    "[fresh kernel] EFT1 running heavy insertion loaded. Terms: ",
    Lookup[
      ExpressionDiagnostics[eft1WilsonRunningHeavyInsertion],
      "TermCount",
      "?"
    ]
  ];
];

currentTree = Quiet @ Check[Import[treePath, "WXF"], $Failed];
currentLoop = Quiet @ Check[Import[loopPath, "WXF"], $Failed];
transitionTree =
  Quiet @ Check[Import[transitionTreePath, "WXF"], $Failed];
transitionFull =
  Quiet @ Check[Import[transitionFullPath, "WXF"], $Failed];

If[
  currentTree === $Failed ||
  currentLoop === $Failed ||
  transitionTree === $Failed ||
  transitionFull === $Failed,
  Print["ERROR: failed to import intermediate EFT."];
  Exit[10]
];

(* Rehydrate only after the fresh field/coupling/CG registry exists. *)
Print[
  "[fresh kernel] Rehydrating canonical transition EFTs under fresh registry..."
];

RehydrateMatcheteExpression[expr_] := Quiet @ Check[
  ToExpression[
    ToString[InputForm[Unevaluated[expr]]],
    InputForm
  ],
  $Failed
];

transitionTreeRehydrated =
  RehydrateMatcheteExpression[transitionTree];
transitionFullRehydrated =
  RehydrateMatcheteExpression[transitionFull];

If[
  transitionTreeRehydrated === $Failed ||
  transitionFullRehydrated === $Failed,
  Print["ERROR: failed to rehydrate canonical transition EFTs."];
  Exit[101]
];

transitionTree = transitionTreeRehydrated;
transitionFull = transitionFullRehydrated;

Print[
  "[fresh kernel] Transition rehydration complete. Tree bytes: ",
  ByteCount[transitionTree],
  "; full O(hbar) bytes: ",
  ByteCount[transitionFull]
];

(* Upstream pole provenance.  Inspect the inherited stage-1 O(hbar)
   transition before threshold-2 matching.  This answers whether stage 2 is
   generating a UV pole or merely propagating one already present in EFT1. *)
ClearAll[InheritedHasUVPole];

InheritedHasUVPole[expr_] := Module[{s = ToString[InputForm[expr]]},
  Or[
    StringContainsQ[s, "\\[Epsilon]"],
    StringContainsQ[s, "\\[CurlyEpsilon]"],
    StringContainsQ[s, "ϵ"],
    StringContainsQ[s, "ε"],
    StringContainsQ[s, "1/Epsilon"],
    StringContainsQ[s, "1/CurlyEpsilon"],
    StringContainsQ[s, "Power[Epsilon, -1]"],
    StringContainsQ[s, "Power[CurlyEpsilon, -1]"]
  ]
];

stage1InheritedLoopRaw = Expand[
  hbar Coefficient[
    Expand[transitionFull - transitionTree],
    hbar,
    1
  ]
];

stage1InheritedLoopHasUVPole =
  InheritedHasUVPole[stage1InheritedLoopRaw];

stage1InheritedLoopTerms = If[
  Head[Expand[stage1InheritedLoopRaw]] === Plus,
  List @@ Expand[stage1InheritedLoopRaw],
  {Expand[stage1InheritedLoopRaw]}
];

stage1InheritedPoleTerms = Select[
  stage1InheritedLoopTerms,
  InheritedHasUVPole
];

stage1InheritedPolePath = FileNameJoin[
  {DirectoryName[resultPath], "stage1_inherited_uv_pole_terms.txt"}
];

Export[
  stage1InheritedPolePath,
  ToString[InputForm[stage1InheritedPoleTerms]],
  "String"
];

Print[
  "[fresh kernel] Stage-1 inherited O(hbar) transition contains UV pole: ",
  stage1InheritedLoopHasUVPole,
  "; pole-term count: ",
  Length[stage1InheritedPoleTerms]
];

(* A raw inherited loop expression may still contain unevaluated Matchete
   loop functions.  Evaluate those functions BEFORE threshold-2 matching to
   distinguish:
     (i) a UV pole already encoded in the stage-1 one-loop EFT, from
     (ii) a pole introduced only by threshold-2 tree matching/canonicalisation.
*)
stage1InheritedLoopEvaluated = Quiet @ Check[
  EvaluateLoopFunctions[stage1InheritedLoopRaw],
  $Failed
];

stage1InheritedEvaluatedHasUVPole = If[
  stage1InheritedLoopEvaluated === $Failed,
  False,
  InheritedHasUVPole[stage1InheritedLoopEvaluated]
];

stage1InheritedEvaluatedTerms = If[
  stage1InheritedLoopEvaluated === $Failed,
  {},
  If[
    Head[Expand[stage1InheritedLoopEvaluated]] === Plus,
    List @@ Expand[stage1InheritedLoopEvaluated],
    {Expand[stage1InheritedLoopEvaluated]}
  ]
];

stage1InheritedEvaluatedPoleTerms = Select[
  stage1InheritedEvaluatedTerms,
  InheritedHasUVPole
];

stage1InheritedEvaluatedPolePath = FileNameJoin[
  {
    DirectoryName[resultPath],
    "stage1_inherited_evaluated_uv_pole_terms.txt"
  }
];

Export[
  stage1InheritedEvaluatedPolePath,
  ToString[InputForm[stage1InheritedEvaluatedPoleTerms]],
  "String"
];

Print[
  "[fresh kernel] Stage-1 inherited O(hbar) after EvaluateLoopFunctions ",
  "contains UV pole: ",
  stage1InheritedEvaluatedHasUVPole,
  "; pole-term count: ",
  Length[stage1InheritedEvaluatedPoleTerms]
];

(* Effective-coupling provenance.
   The inherited transition itself is finite, but RunT3Matching's RawEFT
   already contains the pole.  Test whether expanding reconstructed Matchete
   effective couplings BEFORE Match exposes that pole. *)
stage1TransitionExpandedEffective = Quiet @ Check[
  ReplaceEffectiveCouplings[transitionFullCanonical],
  $Failed
];

stage1TransitionExpandedEffectiveHasUVPole = If[
  stage1TransitionExpandedEffective === $Failed,
  False,
  InheritedHasUVPole[stage1TransitionExpandedEffective]
];

stage1ExpandedEffectiveLoop = If[
  stage1TransitionExpandedEffective === $Failed,
  $Failed,
  Expand[
    hbar Coefficient[
      Expand[stage1TransitionExpandedEffective - transitionTree],
      hbar,
      1
    ]
  ]
];

stage1ExpandedEffectiveLoopHasUVPole = If[
  stage1ExpandedEffectiveLoop === $Failed,
  False,
  InheritedHasUVPole[stage1ExpandedEffectiveLoop]
];

stage1ExpandedEffectivePoleTerms = If[
  stage1ExpandedEffectiveLoop === $Failed,
  {},
  Select[
    If[
      Head[Expand[stage1ExpandedEffectiveLoop]] === Plus,
      List @@ Expand[stage1ExpandedEffectiveLoop],
      {Expand[stage1ExpandedEffectiveLoop]}
    ],
    InheritedHasUVPole
  ]
];

stage1ExpandedEffectivePolePath = FileNameJoin[
  {
    DirectoryName[resultPath],
    "stage1_effective_coupling_expanded_uv_pole_terms.txt"
  }
];

Export[
  stage1ExpandedEffectivePolePath,
  ToString[InputForm[stage1ExpandedEffectivePoleTerms]],
  "String"
];

Print[
  "[fresh kernel] ReplaceEffectiveCouplings before [C]: success=",
  !TrueQ[stage1TransitionExpandedEffective === $Failed],
  "; full transition pole=",
  stage1TransitionExpandedEffectiveHasUVPole,
  "; O(hbar) pole=",
  stage1ExpandedEffectiveLoopHasUVPole,
  "; pole-term count=",
  Length[stage1ExpandedEffectivePoleTerms]
];

(* Diagnose which quadratic kinetic sector changes at O(hbar).  This is
   deliberately diagnostic only: do not remove or rewrite any term yet. *)
Print["[fresh kernel] Diagnosing transition kinetic sectors..."];

transitionTreeExpanded = Quiet @ Check[
  ReplaceEffectiveCouplings[transitionTree, Superleading -> True],
  $Failed
];
transitionFullExpanded = Quiet @ Check[
  ReplaceEffectiveCouplings[transitionFull, Superleading -> True],
  $Failed
];

SafeOperatorClass[lag_, fields_List, derivatives_Integer] := Quiet @ Check[
  SelectOperatorClass[lag, fields, derivatives],
  $Failed
];

NonzeroClassDifferenceQ[a_, b_] := Which[
  a === $Failed || b === $Failed, Missing["Unavailable"],
  True, !TrueQ[Expand[a - b] === 0]
];

If[
  transitionTreeExpanded =!= $Failed &&
  transitionFullExpanded =!= $Failed,

  gaugeDiagnostics = <|
    "U1Y(BB)" -> NonzeroClassDifferenceQ[
      SafeOperatorClass[transitionFullExpanded, {B, B}, 0],
      SafeOperatorClass[transitionTreeExpanded, {B, B}, 0]
    ],
    "SU2L(WW)" -> NonzeroClassDifferenceQ[
      SafeOperatorClass[transitionFullExpanded, {W, W}, 0],
      SafeOperatorClass[transitionTreeExpanded, {W, W}, 0]
    ],
    "SU3c(GG)" -> NonzeroClassDifferenceQ[
      SafeOperatorClass[transitionFullExpanded, {G, G}, 0],
      SafeOperatorClass[transitionTreeExpanded, {G, G}, 0]
    ]
  |>;

  scalar1KineticDiagnostics = DeleteDuplicates @ {
    NonzeroClassDifferenceQ[
      SafeOperatorClass[
        transitionFullExpanded,
        {NewScalar1, Bar[NewScalar1]},
        2
      ],
      SafeOperatorClass[
        transitionTreeExpanded,
        {NewScalar1, Bar[NewScalar1]},
        2
      ]
    ],
    NonzeroClassDifferenceQ[
      SafeOperatorClass[
        transitionFullExpanded,
        {NewScalar1, NewScalar1},
        2
      ],
      SafeOperatorClass[
        transitionTreeExpanded,
        {NewScalar1, NewScalar1},
        2
      ]
    ]
  };

  scalar2KineticDiagnostics = DeleteDuplicates @ {
    NonzeroClassDifferenceQ[
      SafeOperatorClass[
        transitionFullExpanded,
        {NewScalar2, Bar[NewScalar2]},
        2
      ],
      SafeOperatorClass[
        transitionTreeExpanded,
        {NewScalar2, Bar[NewScalar2]},
        2
      ]
    ],
    NonzeroClassDifferenceQ[
      SafeOperatorClass[
        transitionFullExpanded,
        {NewScalar2, NewScalar2},
        2
      ],
      SafeOperatorClass[
        transitionTreeExpanded,
        {NewScalar2, NewScalar2},
        2
      ]
    ]
  };

  Print[
    "[fresh kernel] O(hbar) gauge-kinetic differences: ",
    InputForm[gaugeDiagnostics]
  ];
  Print[
    "[fresh kernel] O(hbar) S1 two-derivative differences: ",
    InputForm[scalar1KineticDiagnostics]
  ];
  Print[
    "[fresh kernel] O(hbar) S2 two-derivative differences: ",
    InputForm[scalar2KineticDiagnostics]
  ],

  Print[
    "[fresh kernel] WARNING: could not expand effective couplings for ",
    "kinetic-sector diagnostics."
  ]
];


(* Gauge-coupling threshold canonicalisation.

   Matchete normalises gauge theories with gauge couplings multiplying the
   gauge kinetic operators rather than sitting inside the covariant
   derivative.  Therefore an O(hbar) threshold correction to F_{\mu\nu}^2
   is a renormalisable gauge-coupling matching correction.

   Before a later Match, absorb that correction into the definition of the
   low-energy gauge couplings by restoring the tree-level gauge-kinetic
   form.  The correction is NOT discarded: it is retained explicitly in
   gaugeKineticThresholdCorrection as the threshold matching relation.

   At the one-loop order used here, replacing a high-energy gauge coupling
   by its threshold-shifted low-energy value inside quantities already of
   O(hbar) would change only O(hbar^2) terms. *)

GaugeClass[lag_, fields_List] := Quiet @ Check[
  SelectOperatorClass[lag, fields, 0],
  0
];

transitionGaugeTree =
  GaugeClass[transitionTreeExpanded, {B, B}] +
  GaugeClass[transitionTreeExpanded, {W, W}] +
  GaugeClass[transitionTreeExpanded, {G, G}];

transitionGaugeFull =
  GaugeClass[transitionFullExpanded, {B, B}] +
  GaugeClass[transitionFullExpanded, {W, W}] +
  GaugeClass[transitionFullExpanded, {G, G}];

gaugeKineticThresholdCorrection = Expand[
  transitionGaugeFull - transitionGaugeTree
];

Print[
  "[fresh kernel] Gauge-coupling threshold shift present: ",
  !TrueQ[gaugeKineticThresholdCorrection === 0]
];

transitionFullCanonical = Expand[
  transitionFullExpanded - gaugeKineticThresholdCorrection
];

Print[
  "[fresh kernel] Absorbing O(hbar) gauge-kinetic threshold correction ",
  "into low-energy gauge couplings before [C]."
];

(* Do not call CheckLagrangian directly on the imported serialized EFT.
   For some representations Matchete's checker can Abort[] the kernel rather
   than returning a recoverable diagnostic.  Match itself performs the same
   consistency checks and gives us a usable failure message. *)

(* ---------------------------------------------------------------------- *)
(* Weinberg-coefficient provenance diagnostics.

   These helpers are diagnostic only.  They do not modify the matching
   result.  They let us isolate the C5 generated by [A], [B], inherited [C],
   and running [C], and in particular identify which perturbative source
   carries any explicit 1/epsilon pole. *)

ClearAll[
  C5TermCount,
  C5HasUVPole,
  C5PieceDiagnostic,
  ExportC5Piece
];

C5TermCount[expr_] := Module[{expanded = Expand[expr]},
  Which[
    TrueQ[expanded === 0], 0,
    Head[expanded] === Plus, Length[List @@ expanded],
    True, 1
  ]
];

C5HasUVPole[expr_] := Module[{text},
  text = ToString[InputForm[expr]];
  Or[
    StringContainsQ[text, "\\[Epsilon]"],
    StringContainsQ[text, "\\[CurlyEpsilon]"],
    StringContainsQ[text, "ϵ"],
    StringContainsQ[text, "ε"],
    StringContainsQ[text, "1/Epsilon"],
    StringContainsQ[text, "1/CurlyEpsilon"],
    StringContainsQ[text, "Power[Epsilon, -1]"],
    StringContainsQ[text, "Power[CurlyEpsilon, -1]"]
  ]
];

C5PieceDiagnostic[label_String, expr_] := Module[
  {extraction, coefficient, status, present},

  extraction = Quiet @ Check[
    ExtractWeinbergCoefficient[expr],
    <|"Status" -> "Failed", "Present" -> False|>
  ];

  status = If[
    AssociationQ[extraction],
    Lookup[extraction, "Status", "Failed"],
    "Failed"
  ];

  present = TrueQ[
    AssociationQ[extraction] &&
    Lookup[extraction, "Present", False]
  ];

  coefficient = If[
    AssociationQ[extraction] && status === "Success",
    Lookup[extraction, "Coefficient", 0],
    0
  ];

  <|
    "Label" -> label,
    "ExtractionStatus" -> status,
    "WeinbergPresent" -> present,
    "ContainsHbar" -> !FreeQ[coefficient, hbar],
    "ContainsUVPole" -> C5HasUVPole[coefficient],
    "CoefficientTermCount" -> C5TermCount[coefficient],
    "CoefficientInputForm" -> ToString[InputForm[coefficient]]
  |>
];

ExportC5Piece[path_, diagnostic_Association] := Export[
  path,
  Lookup[diagnostic, "CoefficientInputForm", "0"],
  "String"
];

(* A. Tree match of L^(0). *)
Print["[fresh kernel] [A] Tree matching L^(0)..."];
treeResult = RunT3Matching[transitionTree, eftOrder, 0];
If[
  !AssociationQ[treeResult] ||
    Lookup[treeResult, "Status", ""] =!= "Success",
  Print["ERROR: fresh-kernel tree matching failed."];
  Exit[11]
];

treeExplicit = CompactExplicitEFT[treeResult];
If[treeExplicit === $Failed,
  Print["ERROR: failed to construct explicit tree EFT."];
  Exit[12]
];

(* B. New one-loop threshold contribution M^(1)[L^(0)]. *)
localLoopCorrection = 0;
oneLoopResult = treeResult;
bSanityPath = Missing["NotRun"];

If[loopOrder >= 1,
  Print["[fresh kernel] [B] One-loop matching L^(0)..."];
  oneLoopResult = RunT3Matching[transitionTree, eftOrder, 1];

  If[
    !AssociationQ[oneLoopResult] ||
      Lookup[oneLoopResult, "Status", ""] =!= "Success",
    Print["ERROR: fresh-kernel one-loop matching failed."];
    Exit[13]
  ];

  fullFromTreeExplicit = CompactExplicitEFT[oneLoopResult];
  If[fullFromTreeExplicit === $Failed,
    Print["ERROR: failed to construct explicit one-loop EFT."];
    Exit[14]
  ];

  localLoopCorrection =
    Expand[fullFromTreeExplicit - treeExplicit];

  (* [B] sanity diagnostic.
     After integrating out F first, the tree EFT should already contain the
     LL S1 S2 operator generated by fermion exchange, while lambdaT3 is also
     present.  Their scalar one-loop matching is the expected T3 source of
     the Weinberg operator.  Record whether these ingredients are actually
     present in the stage-2 tree input and whether [B] generates C5. *)
  bTreeInputText = ToString[InputForm[transitionTree]];
  bRawText = ToString[
    InputForm[Lookup[oneLoopResult, "RawEFT", 0]]
  ];

  bTreeIngredientDiagnostic = <|
    "HasY1" -> StringContainsQ[bTreeInputText, "y1"],
    "HasY2" -> StringContainsQ[bTreeInputText, "y2"],
    "HasLambdaT3" -> StringContainsQ[bTreeInputText, "lambdaT3"],
    "HasS1" -> Or[
      StringContainsQ[bTreeInputText, "NewScalar1"],
      StringContainsQ[bTreeInputText, "S1"]
    ],
    "HasS2" -> Or[
      StringContainsQ[bTreeInputText, "NewScalar2"],
      StringContainsQ[bTreeInputText, "S2"]
    ]
  |>;

  bRawIngredientDiagnostic = <|
    "HasY1" -> StringContainsQ[bRawText, "y1"],
    "HasY2" -> StringContainsQ[bRawText, "y2"],
    "HasLambdaT3" -> StringContainsQ[bRawText, "lambdaT3"]
  |>;

  bC5Diagnostic = C5PieceDiagnostic[
    "B_one_loop_raw",
    Lookup[oneLoopResult, "RawEFT", 0]
  ];

  Print["[fresh kernel] [B] expected-topology sanity check:"];
  Print[
    "  tree input: y1=",
    Lookup[bTreeIngredientDiagnostic, "HasY1", False],
    ", y2=",
    Lookup[bTreeIngredientDiagnostic, "HasY2", False],
    ", lambdaT3=",
    Lookup[bTreeIngredientDiagnostic, "HasLambdaT3", False],
    ", S1=",
    Lookup[bTreeIngredientDiagnostic, "HasS1", False],
    ", S2=",
    Lookup[bTreeIngredientDiagnostic, "HasS2", False]
  ];
  Print[
    "  [B] RawEFT: y1=",
    Lookup[bRawIngredientDiagnostic, "HasY1", False],
    ", y2=",
    Lookup[bRawIngredientDiagnostic, "HasY2", False],
    ", lambdaT3=",
    Lookup[bRawIngredientDiagnostic, "HasLambdaT3", False],
    "; Weinberg present=",
    Lookup[bC5Diagnostic, "WeinbergPresent", False]
  ];

  bSanityPath = FileNameJoin[
    {DirectoryName[resultPath], "stage2_B_sanity.json"}
  ];

  Export[
    bSanityPath,
    <|
      "TreeInput" -> bTreeIngredientDiagnostic,
      "BRawEFT" -> bRawIngredientDiagnostic,
      "BC5" -> bC5Diagnostic
    |>,
    "JSON"
  ];
];

(* C. Inherited one-loop contribution M^(0)[L^(1)].
   Match the already-canonical full EFT only at tree order, then keep the
   coefficient linear in Matchete's native hbar. *)
inheritedLoopCorrection = 0;

If[
  currentLoop =!= 0 || !TrueQ[eft1WilsonRunningHeavyInsertion === 0],

  Print[
    "[fresh kernel] [C] Tree propagation of inherited/running O(hbar) EFT..."
  ];

  (* The full EFT was canonicalised in the parent kernel while S1/S2 were
     still Light.  Do NOT rebuild it here as currentTree + currentLoop:
     that sum contains the pre-field-redefinition one-loop kinetic terms
     that caused CheckLagrangian::CanonicallyNormalized. *)
  (* transitionFull is the EOMSimplify-canonical representation from the
     previous kernel, with effective couplings intentionally NOT expanded. *)
  fullInput = Expand[
    transitionFullCanonical + eft1WilsonRunningHeavyInsertion
  ];

  If[eft1WilsonRunningLoaded,
    Print[
      "[fresh kernel] Added EFT1 leading-log heavy insertion to [C] only."
    ]
  ];

  (* This probe has expanded effective couplings already.  The explicit
     hbar bookkeeping is therefore visible before/after tree matching. *)
  probeTreeResult = RunT3Matching[fullInput, eftOrder, 0];

  If[
    !AssociationQ[probeTreeResult] ||
      Lookup[probeTreeResult, "Status", ""] =!= "Success",
    Print["ERROR: inherited-loop tree propagation failed."];
    Exit[16]
  ];

  probeExplicit = CompactExplicitEFT[probeTreeResult];

  (* The dummy EFT1RunLog coupling exists only to keep CheckLagrangian happy.
     Once matching/simplification is complete, restore the actual threshold
     logarithm. *)
  If[Head[eft1WilsonRunningLogReplacement] === Rule,
    probeExplicit =
      Expand[probeExplicit /. eft1WilsonRunningLogReplacement];

    Print[
      "[fresh kernel] Restored physical EFT1 running logarithm after [C]."
    ];
  ];

  If[probeExplicit === $Failed,
    Print["ERROR: inherited-loop explicit EFT failed."];
    Exit[17]
  ];

  If[FreeQ[probeExplicit, hbar],
    Print[
      "ERROR: inherited-loop result still contains no explicit hbar after ",
      "effective-coupling expansion; fixed-order projection is unsafe."
    ];
    Exit[15]
  ];

  Print[
    "[fresh kernel] Explicit hbar restored after effective-coupling expansion."
  ];

  inheritedLoopCorrection = Expand[
    hbar Coefficient[
      Expand[probeExplicit - treeExplicit],
      hbar,
      1
    ]
  ];

  Print[
    "[fresh kernel] Inherited one-loop terms: ",
    Lookup[
      ExpressionDiagnostics[inheritedLoopCorrection],
      "TermCount",
      "?"
    ]
  ];
];

(* For provenance only, split [C] into:
     [C_inherited] = M2^(0)[L1_threshold^(1)]
     [C_running]   = M2^(0)[L1_run^(1)]

   The authoritative [C] above is still the actual result used downstream.
   When a running insertion is present, perform one extra tree match with the
   inherited transition alone and obtain the running piece by subtraction.
   This keeps the diagnostic aligned with the exact same Matchete
   canonicalisation path used by the authoritative combined [C]. *)
inheritedThresholdLoopCorrection = inheritedLoopCorrection;
runningHeavyLoopCorrection = 0;

If[eft1WilsonRunningLoaded,

  Print[
    "[fresh kernel] Provenance diagnostic: isolating inherited [C] from ",
    "running [C]..."
  ];

  inheritedOnlyResult = RunT3Matching[
    transitionFullCanonical,
    eftOrder,
    0
  ];

  If[
    !AssociationQ[inheritedOnlyResult] ||
      Lookup[inheritedOnlyResult, "Status", ""] =!= "Success",
    Print[
      "ERROR: provenance diagnostic failed while matching inherited-only [C]."
    ];
    Exit[22]
  ];

  (* Stage-by-stage provenance inside RunT3Matching.
     RunMatching.wl retains every intermediate representation under:
       RawEFT, GreenEFT, EOMEFT, LoopEFT, MatchedEFT.
     Inspect each one before CompactExplicitEFT changes anything further. *)
  c5InheritedPipelineDiagnostics = Association @ Map[
    Function[key,
      key -> If[
        KeyExistsQ[inheritedOnlyResult, key],
        C5PieceDiagnostic[
          "C_inherited_" <> key,
          inheritedOnlyResult[key]
        ],
        <|
          "Label" -> ("C_inherited_" <> key),
          "ExtractionStatus" -> "Missing",
          "WeinbergPresent" -> False,
          "ContainsHbar" -> False,
          "ContainsUVPole" -> False,
          "CoefficientTermCount" -> 0,
          "CoefficientInputForm" -> "0"
        |>
      ]
    ],
    {"RawEFT", "GreenEFT", "EOMEFT", "LoopEFT", "MatchedEFT"}
  ];

  Print["[fresh kernel] [C inherited] stage-by-stage C5 provenance:"];
  Do[
    diag = c5InheritedPipelineDiagnostics[key];
    Print[
      "  ", key,
      ": present=", Lookup[diag, "WeinbergPresent", False],
      "; pole=", Lookup[diag, "ContainsUVPole", False],
      "; hbar=", Lookup[diag, "ContainsHbar", False],
      "; terms=", Lookup[diag, "CoefficientTermCount", 0]
    ];
    If[
      TrueQ[Lookup[diag, "WeinbergPresent", False]],
      Print[
        "    coefficient=",
        Lookup[diag, "CoefficientInputForm", "0"]
      ]
    ],
    {key, {"RawEFT", "GreenEFT", "EOMEFT", "LoopEFT", "MatchedEFT"}}
  ];

  c5InheritedPipelinePath = FileNameJoin[
    {
      DirectoryName[resultPath],
      "c5_stage2_C_inherited_pipeline.json"
    }
  ];

  Export[
    c5InheritedPipelinePath,
    c5InheritedPipelineDiagnostics,
    "JSON"
  ];

  inheritedOnlyExplicit = CompactExplicitEFT[inheritedOnlyResult];

  If[inheritedOnlyExplicit === $Failed,
    Print[
      "ERROR: provenance diagnostic failed to construct inherited-only EFT."
    ];
    Exit[23]
  ];

  If[FreeQ[inheritedOnlyExplicit, hbar],
    Print[
      "ERROR: provenance inherited-only result contains no explicit hbar; ",
      "cannot isolate [C] pieces safely."
    ];
    Exit[24]
  ];

  inheritedThresholdLoopCorrection = Expand[
    hbar Coefficient[
      Expand[inheritedOnlyExplicit - treeExplicit],
      hbar,
      1
    ]
  ];

  runningHeavyLoopCorrection = Expand[
    inheritedLoopCorrection - inheritedThresholdLoopCorrection
  ];

  Print[
    "[fresh kernel] Provenance split complete: inherited [C] terms=",
    Lookup[
      ExpressionDiagnostics[inheritedThresholdLoopCorrection],
      "TermCount",
      "?"
    ],
    "; running [C] terms=",
    Lookup[
      ExpressionDiagnostics[runningHeavyLoopCorrection],
      "TermCount",
      "?"
    ]
  ];

  (* Additional diagnostic:
     match an input in which effective couplings have already been expanded.
     If the RawEFT still contains the same 1/epsilon term, then the pole is
     generated by Match itself from the explicit finite inherited EFT rather
     than by hidden effective-coupling definitions. *)
  If[
    stage1TransitionExpandedEffective =!= $Failed,

    Print[
      "[fresh kernel] Provenance diagnostic: matching inherited [C] with ",
      "effective couplings pre-expanded..."
    ];

    inheritedExpandedInputResult = RunT3Matching[
      stage1TransitionExpandedEffective,
      eftOrder,
      0
    ];

    If[
      AssociationQ[inheritedExpandedInputResult] &&
        Lookup[inheritedExpandedInputResult, "Status", ""] === "Success",

      inheritedExpandedRawDiagnostic = C5PieceDiagnostic[
        "C_inherited_expanded_input_RawEFT",
        Lookup[inheritedExpandedInputResult, "RawEFT", 0]
      ];

      Print[
        "[fresh kernel] Expanded-input RawEFT C5: present=",
        Lookup[
          inheritedExpandedRawDiagnostic,
          "WeinbergPresent",
          False
        ],
        "; pole=",
        Lookup[
          inheritedExpandedRawDiagnostic,
          "ContainsUVPole",
          False
        ],
        "; terms=",
        Lookup[
          inheritedExpandedRawDiagnostic,
          "CoefficientTermCount",
          0
        ]
      ];

      Print[
        "    coefficient=",
        Lookup[
          inheritedExpandedRawDiagnostic,
          "CoefficientInputForm",
          "0"
        ]
      ],

      inheritedExpandedRawDiagnostic = <|
        "ExtractionStatus" -> "MatchFailed",
        "WeinbergPresent" -> False,
        "ContainsUVPole" -> False,
        "CoefficientTermCount" -> 0,
        "CoefficientInputForm" -> "0"
      |>;

      Print[
        "[fresh kernel] Expanded-input inherited [C] match failed; ",
        "diagnostic inconclusive."
      ]
    ],

    inheritedExpandedRawDiagnostic = <|
      "ExtractionStatus" -> "ExpansionFailed",
      "WeinbergPresent" -> False,
      "ContainsUVPole" -> False,
      "CoefficientTermCount" -> 0,
      "CoefficientInputForm" -> "0"
    |>
  ];
];

(* D. Pure-SM Weinberg component generated directly by EFT1 running.

   This piece contains no S1/S2 fields, so threshold-2 matching is the
   identity on it.  Importantly, the final Matchete Weinberg sector used by
   ExtractWeinbergCoefficient has already contracted/stripped the explicit
   Flavor indices from the operator structure.  Therefore trying to rebuild a
   full-flavor LLHH operator from that final Matchete expression is not
   well-defined: the flavor information belongs to the coefficient, not to
   the stripped operator template.

   Carry the full-flavor coefficient separately and combine it with the
   threshold-generated C5 at the coefficient/report layer.  This preserves
   the p,q flavor dependence exactly and avoids inventing a Matchete operator
   convention that is no longer present in the compact final EFT. *)

nextTree = treeExplicit;
nextLoop = Expand[localLoopCorrection + inheritedLoopCorrection];
nextCompact = Expand[nextTree + nextLoop];

Print[
  "[fresh kernel] Direct Weinberg definition available: ",
  eft1WilsonRunningWeinbergDefinition =!= {}
];

If[
  eft1WilsonRunningLoaded &&
  eft1WilsonRunningWeinbergDefinition =!= {},

  Print[
    "[fresh kernel] Carrying direct EFT1-generated Weinberg running " ,
    "coefficient separately from threshold matching..."
  ];

  (* eft1WilsonRunningWeinbergCoefficient was sampled at formal Flavor
     indices when the insertion file was loaded.  hbar is restored here
     because the exported beta/transport coefficient is the O(hbar^0)
     coefficient multiplying the one-loop correction. *)
  eft1WilsonRunningDirectWeinberg = Expand[
    hbar * eft1WilsonRunningWeinbergCoefficient
  ];

  eft1WilsonRunningDirectWeinbergEqualScaleVanishes = TrueQ[
    Quiet @ Check[
      Simplify[
        Expand[
          eft1WilsonRunningDirectWeinberg /.
            MS -> Coupling[MF, {}, 0]
        ]
      ] === 0,
      False
    ]
  ];

  If[!TrueQ[eft1WilsonRunningDirectWeinbergEqualScaleVanishes],
    Print[
      "ERROR: direct EFT1 Weinberg running coefficient does not vanish " ,
      "at equal scales."
    ];
    Exit[20]
  ];

  eft1WilsonRunningDirectWeinbergAdded = True;

  Print[
    "[fresh kernel] Direct EFT1 Weinberg running coefficient carried " ,
    "separately to final C5."
  ];
  Print[
    "[fresh kernel] Direct Weinberg equal-scale check: ",
    eft1WilsonRunningDirectWeinbergEqualScaleVanishes
  ];
];

(* Extract the threshold-generated Weinberg coefficient from the
   AUTHORITATIVE resumed threshold-2 result.  This includes:
     [B] new one-loop matching at the S1/S2 threshold,
     [C] tree propagation of inherited stage-1 O(hbar) terms, and
     [C] tree propagation of the EFT1 heavy-field running insertion.

   The pure-SM Weinberg running coefficient remains separate because it does
   not pass through threshold-2 matching. *)
authoritativeWeinbergExtraction = Quiet @ Check[
  ExtractWeinbergCoefficient[nextCompact],
  <|"Status" -> "Failed", "Present" -> False|>
];

authoritativeThresholdC5 = If[
  AssociationQ[authoritativeWeinbergExtraction] &&
  Lookup[authoritativeWeinbergExtraction, "Status", "Failed"] === "Success",
  Lookup[authoritativeWeinbergExtraction, "Coefficient", 0],
  0
];

authoritativeC5Dir = DirectoryName[resultPath];
authoritativeThresholdC5Path = FileNameJoin[
  {authoritativeC5Dir, "c5_threshold_with_eft1_running.txt"}
];
authoritativeDirectC5Path = FileNameJoin[
  {authoritativeC5Dir, "c5_direct_eft1_running.txt"}
];

If[
  authoritativeThresholdC5 =!= 0,
  Export[
    authoritativeThresholdC5Path,
    ToString[InputForm[authoritativeThresholdC5]],
    "String"
  ];
  Print[
    "[fresh kernel] Authoritative resumed threshold C5 exported: ",
    authoritativeThresholdC5Path
  ],
  Print[
    "ERROR: could not extract authoritative Weinberg coefficient from ",
    "resumed threshold-2 EFT."
  ];
  Exit[21]
];

Export[
  authoritativeDirectC5Path,
  ToString[InputForm[eft1WilsonRunningDirectWeinberg]],
  "String"
];

Print[
  "[fresh kernel] Direct EFT1 Weinberg running coefficient exported: ",
  authoritativeDirectC5Path
];

(* ---------------------------------------------------------------------- *)
(* Pole/RGE consistency diagnostic.

   The inherited stage-1 hard contribution contains a 1/epsilon pole, while
   the separately constructed EFT1 running gives the finite logarithm.
   In the one-generation reduction expected for this T3-B benchmark, the
   coefficient of Log[MS/MF] should be twice the residue of the hard-region
   1/epsilon pole.  That is exactly the relation needed for the usual
   hard/soft cancellation:
       r [1/eps + log(mu^2/MF^2) + finite]
     - r [1/eps + log(mu^2/MS^2)]
       = finite + 2 r log(MS/MF).

   This is a diagnostic only: it does not subtract the pole or alter the
   authoritative coefficient. *)

ClearAll[CollapseT3FlavorIndices];

CollapseT3FlavorIndices[expr_] := Expand[
  expr /.
    Index[_, Flavor] -> Index[T3Flavor1, Flavor] /.
    Index[_, NFlavor] -> Index[T3NFlavor1, NFlavor]
];

c5HardPoleResidue = Quiet @ Check[
  Coefficient[
    Expand[authoritativeThresholdC5],
    1/\[Epsilon]
  ],
  $Failed
];

c5DirectLogCoefficient = Quiet @ Check[
  Coefficient[
    Expand[eft1WilsonRunningDirectWeinberg],
    Log[MS/Coupling[MF, {}, 0]]
  ],
  $Failed
];

c5HardPoleResidueOneGen = If[
  c5HardPoleResidue === $Failed,
  $Failed,
  CollapseT3FlavorIndices[c5HardPoleResidue]
];

c5DirectLogCoefficientOneGen = If[
  c5DirectLogCoefficient === $Failed,
  $Failed,
  CollapseT3FlavorIndices[c5DirectLogCoefficient]
];

c5PoleRGEConsistency = If[
  MemberQ[
    {
      c5HardPoleResidueOneGen,
      c5DirectLogCoefficientOneGen
    },
    $Failed
  ],
  False,
  TrueQ[
    Quiet @ Check[
      Simplify[
        Expand[
          c5DirectLogCoefficientOneGen -
            2 c5HardPoleResidueOneGen
        ]
      ] === 0,
      False
    ]
  ]
];

c5PoleRGERatio = If[
  MemberQ[
    {
      c5HardPoleResidueOneGen,
      c5DirectLogCoefficientOneGen
    },
    $Failed | 0
  ],
  Missing["Undefined"],
  Quiet @ Check[
    FullSimplify[
      c5DirectLogCoefficientOneGen /
        c5HardPoleResidueOneGen
    ],
    Missing["SimplifyFailed"]
  ]
];

c5PoleRGEPath = FileNameJoin[
  {authoritativeC5Dir, "c5_pole_rge_consistency.json"}
];

c5PoleRGEDiagnostic = <|
  "Status" -> "Success",
  "HardPoleResidueInputForm" ->
    ToString[InputForm[c5HardPoleResidue]],
  "DirectRunningLogCoefficientInputForm" ->
    ToString[InputForm[c5DirectLogCoefficient]],
  "HardPoleResidueOneGenerationInputForm" ->
    ToString[InputForm[c5HardPoleResidueOneGen]],
  "DirectRunningLogCoefficientOneGenerationInputForm" ->
    ToString[InputForm[c5DirectLogCoefficientOneGen]],
  "DirectLogToPoleRatioOneGenerationInputForm" ->
    ToString[InputForm[c5PoleRGERatio]],
  "DirectLogEqualsTwicePoleResidueOneGeneration" ->
    c5PoleRGEConsistency
|>;

Export[c5PoleRGEPath, c5PoleRGEDiagnostic, "JSON"];

Print["[fresh kernel] C5 pole/RGE consistency:"];
Print[
  "  direct-log coefficient / hard-pole residue = ",
  c5PoleRGERatio
];
Print[
  "  direct-log coefficient == 2 * hard-pole residue: ",
  c5PoleRGEConsistency
];
Print[
  "  diagnostic exported: ",
  c5PoleRGEPath
];

(* Export perturbative-source C5 diagnostics.  These files answer the
   renormalisation question without changing the final result:

     A              : tree threshold matching
     B_one_loop     : new one-loop threshold matching
     C_inherited    : tree propagation of inherited stage-1 O(hbar)
     C_running      : tree propagation of EFT1 leading-log heavy running

   The direct pure-SM Weinberg running term is already exported separately
   above because it bypasses threshold-2 matching. *)
c5ProvenanceDir = DirectoryName[resultPath];

c5Stage2APath = FileNameJoin[
  {c5ProvenanceDir, "c5_stage2_A.txt"}
];
c5Stage2BPath = FileNameJoin[
  {c5ProvenanceDir, "c5_stage2_B_one_loop.txt"}
];
c5Stage2CInheritedPath = FileNameJoin[
  {c5ProvenanceDir, "c5_stage2_C_inherited.txt"}
];
c5Stage2CRunningPath = FileNameJoin[
  {c5ProvenanceDir, "c5_stage2_C_running.txt"}
];
c5Stage2ProvenancePath = FileNameJoin[
  {c5ProvenanceDir, "c5_stage2_provenance.json"}
];

c5Stage2ADiagnostic =
  C5PieceDiagnostic["A_tree", treeExplicit];

c5Stage2BDiagnostic =
  C5PieceDiagnostic["B_one_loop", localLoopCorrection];

c5Stage2CInheritedDiagnostic =
  C5PieceDiagnostic[
    "C_inherited",
    inheritedThresholdLoopCorrection
  ];

c5Stage2CRunningDiagnostic =
  C5PieceDiagnostic[
    "C_running",
    runningHeavyLoopCorrection
  ];

ExportC5Piece[c5Stage2APath, c5Stage2ADiagnostic];
ExportC5Piece[c5Stage2BPath, c5Stage2BDiagnostic];
ExportC5Piece[
  c5Stage2CInheritedPath,
  c5Stage2CInheritedDiagnostic
];
ExportC5Piece[
  c5Stage2CRunningPath,
  c5Stage2CRunningDiagnostic
];

c5Stage2Provenance = <|
  "Status" -> "Success",
  "Stage1InheritedTransition" -> <|
    "ContainsUVPole" -> stage1InheritedLoopHasUVPole,
    "PoleTermCount" -> Length[stage1InheritedPoleTerms],
    "PoleTermsPath" -> stage1InheritedPolePath,
    "EvaluateLoopFunctionsSucceeded" ->
      !TrueQ[stage1InheritedLoopEvaluated === $Failed],
    "ContainsUVPoleAfterEvaluateLoopFunctions" ->
      stage1InheritedEvaluatedHasUVPole,
    "EvaluatedPoleTermCount" ->
      Length[stage1InheritedEvaluatedPoleTerms],
    "EvaluatedPoleTermsPath" ->
      stage1InheritedEvaluatedPolePath,
    "ReplaceEffectiveCouplingsBeforeC" -> <|
      "Succeeded" ->
        !TrueQ[stage1TransitionExpandedEffective === $Failed],
      "FullTransitionContainsUVPole" ->
        stage1TransitionExpandedEffectiveHasUVPole,
      "LoopPieceContainsUVPole" ->
        stage1ExpandedEffectiveLoopHasUVPole,
      "PoleTermCount" ->
        Length[stage1ExpandedEffectivePoleTerms],
      "PoleTermsPath" ->
        stage1ExpandedEffectivePolePath
    |>
  |>,
  "A" -> c5Stage2ADiagnostic,
  "B" -> c5Stage2BDiagnostic,
  "CInherited" -> c5Stage2CInheritedDiagnostic,
  "CRunning" -> c5Stage2CRunningDiagnostic,
  "DirectWeinbergRunning" -> <|
    "Present" -> !TrueQ[eft1WilsonRunningDirectWeinberg === 0],
    "ContainsHbar" -> !FreeQ[eft1WilsonRunningDirectWeinberg, hbar],
    "ContainsUVPole" -> C5HasUVPole[
      eft1WilsonRunningDirectWeinberg
    ],
    "CoefficientTermCount" -> C5TermCount[
      eft1WilsonRunningDirectWeinberg
    ],
    "CoefficientInputForm" ->
      ToString[InputForm[eft1WilsonRunningDirectWeinberg]]
  |>,
  "AuthoritativeThresholdC5" -> <|
    "ContainsHbar" -> !FreeQ[authoritativeThresholdC5, hbar],
    "ContainsUVPole" -> C5HasUVPole[authoritativeThresholdC5],
    "CoefficientTermCount" -> C5TermCount[
      authoritativeThresholdC5
    ],
    "CoefficientInputForm" ->
      ToString[InputForm[authoritativeThresholdC5]]
  |>
|>;

Export[c5Stage2ProvenancePath, c5Stage2Provenance, "JSON"];

Print["[fresh kernel] C5 provenance diagnostics:"];
Print[
  "  [A] present=", Lookup[c5Stage2ADiagnostic, "WeinbergPresent", False],
  "; pole=", Lookup[c5Stage2ADiagnostic, "ContainsUVPole", False],
  "; terms=", Lookup[c5Stage2ADiagnostic, "CoefficientTermCount", 0]
];
Print[
  "  [B] present=", Lookup[c5Stage2BDiagnostic, "WeinbergPresent", False],
  "; pole=", Lookup[c5Stage2BDiagnostic, "ContainsUVPole", False],
  "; terms=", Lookup[c5Stage2BDiagnostic, "CoefficientTermCount", 0]
];
Print[
  "  [C inherited] present=",
  Lookup[c5Stage2CInheritedDiagnostic, "WeinbergPresent", False],
  "; pole=",
  Lookup[c5Stage2CInheritedDiagnostic, "ContainsUVPole", False],
  "; terms=",
  Lookup[c5Stage2CInheritedDiagnostic, "CoefficientTermCount", 0]
];
Print[
  "  [C running] present=",
  Lookup[c5Stage2CRunningDiagnostic, "WeinbergPresent", False],
  "; pole=",
  Lookup[c5Stage2CRunningDiagnostic, "ContainsUVPole", False],
  "; terms=",
  Lookup[c5Stage2CRunningDiagnostic, "CoefficientTermCount", 0]
];
Print[
  "[fresh kernel] C5 provenance JSON exported: ",
  c5Stage2ProvenancePath
];

nextReport = SafeStage[EvaluateLoopFunctions[nextCompact]];
If[MemberQ[{$Failed, $Aborted}, nextReport],
  nextReport = nextCompact
];

(* Prepare the canonical representations to be consumed by a still-later
   threshold.  The tree transition comes directly from the tree match's
   EOMEFT.  The full transition is canonicalised from the consistently
   truncated one-loop result while the surviving fields are still Light in
   this stage's output. *)
nextTransitionTree = CanonicalTransitionEFT[treeResult];

nextTransitionFull = PrepareSequentialMatchInput[nextCompact];

If[
  nextTransitionTree === $Failed ||
  nextTransitionFull === $Failed,
  Print["ERROR: failed to prepare canonical transition EFT for next stage."];
  Exit[18]
];

result = <|
  "Status" -> "Success",
  "TreeLagrangian" -> nextTree,
  "OneLoopCorrection" -> nextLoop,
  "GaugeKineticThresholdCorrection" -> gaugeKineticThresholdCorrection,
  "EFT1WilsonRunningLoaded" -> eft1WilsonRunningLoaded,
  "EFT1WilsonRunningHeavyInsertion" -> eft1WilsonRunningHeavyInsertion,
  "EFT1WilsonRunningWeinbergCoefficient" ->
    eft1WilsonRunningWeinbergCoefficient,
  "EFT1WilsonRunningWeinbergDefinition" ->
    eft1WilsonRunningWeinbergDefinition,
  "EFT1WilsonRunningLogReplacement" ->
    eft1WilsonRunningLogReplacement,
  "EFT1WilsonRunningDirectWeinberg" ->
    eft1WilsonRunningDirectWeinberg,
  "EFT1WilsonRunningDirectWeinbergCoefficientContribution" ->
    eft1WilsonRunningDirectWeinberg,
  "EFT1WilsonRunningDirectWeinbergAdded" ->
    eft1WilsonRunningDirectWeinbergAdded,
  "EFT1WilsonRunningDirectWeinbergCarriedSeparately" ->
    eft1WilsonRunningDirectWeinbergAdded,
  "EFT1WilsonRunningDirectWeinbergEqualScaleVanishes" ->
    eft1WilsonRunningDirectWeinbergEqualScaleVanishes,
  "AuthoritativeWeinbergExtraction" ->
    authoritativeWeinbergExtraction,
  "AuthoritativeThresholdC5" ->
    authoritativeThresholdC5,
  "AuthoritativeThresholdC5Path" ->
    authoritativeThresholdC5Path,
  "AuthoritativeDirectWeinbergC5Path" ->
    authoritativeDirectC5Path,
  "C5Stage2Provenance" -> c5Stage2Provenance,
  "C5Stage2ProvenancePath" -> c5Stage2ProvenancePath,
  "C5Stage2CInheritedPipeline" -> c5InheritedPipelineDiagnostics,
  "C5Stage2CInheritedPipelinePath" -> c5InheritedPipelinePath,
  "C5Stage2CInheritedExpandedInputRaw" ->
    inheritedExpandedRawDiagnostic,
  "Stage2BSanityPath" -> bSanityPath,
  "C5PoleRGEConsistency" -> c5PoleRGEDiagnostic,
  "C5PoleRGEConsistencyPath" -> c5PoleRGEPath,
  "C5Stage2APath" -> c5Stage2APath,
  "C5Stage2BPath" -> c5Stage2BPath,
  "C5Stage2CInheritedPath" -> c5Stage2CInheritedPath,
  "C5Stage2CRunningPath" -> c5Stage2CRunningPath,
  "TransitionTreeLagrangian" -> nextTransitionTree,
  "TransitionFullLagrangian" -> nextTransitionFull,
  "CompactLagrangian" -> nextCompact,
  "Lagrangian" -> nextReport,
  "FreshKernel" -> True,
  "ActiveFieldsAtStart" -> activeFields,
  "HeavyFieldsAtStage" -> heavyFields
|>;

Export[resultPath, result, "WXF"];

Print["[fresh kernel] Threshold stage completed successfully."];
Exit[0];

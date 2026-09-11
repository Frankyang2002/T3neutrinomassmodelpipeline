(* RunThresholdStage.wl

   Execute ONE sequential T3 threshold in a completely fresh Matchete kernel.

   CLI:
     result.wxf EFTOrder LoopOrder dS1 dS2 dF alpha
     active-fields-csv heavy-fields-csv input-tree.wxf input-loop.wxf

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
];

(* C. Inherited one-loop contribution M^(0)[L^(1)].
   Match the already-canonical full EFT only at tree order, then keep the
   coefficient linear in Matchete's native hbar. *)
inheritedLoopCorrection = 0;

If[currentLoop =!= 0,
  Print["[fresh kernel] [C] Tree propagation of inherited O(hbar) EFT..."];

  (* The full EFT was canonicalised in the parent kernel while S1/S2 were
     still Light.  Do NOT rebuild it here as currentTree + currentLoop:
     that sum contains the pre-field-redefinition one-loop kinetic terms
     that caused CheckLagrangian::CanonicallyNormalized. *)
  (* transitionFull is the EOMSimplify-canonical representation from the
     previous kernel, with effective couplings intentionally NOT expanded. *)
  fullInput = transitionFullCanonical;

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

nextTree = treeExplicit;
nextLoop = Expand[localLoopCorrection + inheritedLoopCorrection];
nextCompact = Expand[nextTree + nextLoop];

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

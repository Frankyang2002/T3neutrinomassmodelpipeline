(*
  RunMatching.wl
  Matchete UV -> EFT matching plus Weinberg-operator extraction.
  Intermediate EFTs are retained for diagnostics.

  All this file does is 
  1. We match each stage
*)

ClearAll[
  SafeStage,
  RunStages,
  MatchingStages,
  RunT3Matching,
  RunSequentialT3Matching,
  ExportIntermediateEFTForRGE,
  RunSMBaselineMatching,
  BuildMatchedEFTDifference
];


(* ---------------------------------------------------------------------- *)
(* Matching pipeline                                                       *)
(* ---------------------------------------------------------------------- *)

(* Wrapper to make the code deal with errors.

   Important for sequential matching:
   Matchete can emit diagnostic/warning messages when a theory contains both
   Heavy and Light BSM fields.  Mathematica's Check[expr, $Failed] treats
   *any* message as a failure, even when expr itself returns a perfectly valid
   EFT.  That made the first partial-threshold Match look like a hard failure.

   We therefore distinguish:
     - an actual returned $Failed,
     - an abort,
     - ordinary messages/warnings.

   Messages remain visible in the Wolfram log, but do not by themselves turn
   a successful Matchete result into $Failed. *)
SafeStage[operation_] := CheckAbort[
  operation,
  $Aborted
];

(* Run ordered matching stages while retaining every intermediate result.
 We are running many functions in a row, this a helper function for function running*)
RunStages[input_, stages_List, prefix_Association : <||>] := Module[
  {current = input, results = <||>, stageName, function, key, failureStatus, value},

  Do[
    {stageName, function, key, failureStatus} = specification;
    Print["Running ", stageName, "..."];

    (* Specification looks like {  humanReadableName,  functionToRun,  resultKey,  errorStatus}
    Example is {"Match",  Match[#, EFTOrder -> eftOrder, LoopOrder -> loopOrder] &, "RawEFT",  "MatchFailed"}*)
    value = SafeStage[function[current]];
    results[key] = value; (* We keep intermediate representations, not just leaving the last one *)

    Print[
      "  ", stageName, " result head: ",
      Quiet @ Check[Head[value], "Unknown"]
    ];

    (* Detect genuine returned failures.  A warning message alone is not a
       failure; it is still present in stdout/stderr for inspection. *)
    If[MemberQ[{$Failed, $Aborted}, value],
      Print["ERROR: ", stageName, " failed or aborted."];
      Return[Join[<|"Status" -> failureStatus|>, prefix, results]]
    ];

    (* Output to next stage *)
    current = value,
    {specification, stages}
  ];

  Join[<|"Status" -> "Success"|>, prefix, results]
];

(* We get actual Matchete workflow 
1. First we match 
2. then we green simplify to get green basis, getting rid of redundancies coming from integration by parts and identities etc 
3. EOM simplify remove operator redundancies from EOM 
4. Evaluate loop functions 
5. replace effective couplings with our original couplings*)
MatchingStages[eftOrder_Integer, loopOrder_Integer] := {
  {
    "Match",
    Match[#, EFTOrder -> eftOrder, LoopOrder -> loopOrder] &,
    "RawEFT",
    "MatchFailed"
  },
  {"GreensSimplify", GreensSimplify, "GreenEFT", "GreensSimplifyFailed"},
  {"EOMSimplify", EOMSimplify, "EOMEFT", "EOMSimplifyFailed"},
  {
    "EvaluateLoopFunctions",
    EvaluateLoopFunctions,
    "LoopEFT",
    "LoopEvaluationFailed"
  },
  {
    "ReplaceEffectiveCouplings",
    ReplaceEffectiveCouplings,
    "MatchedEFT",
    "EffectiveCouplingReplacementFailed"
  }
};

(* Convenient Wrapper to get EFT order and loop order, just does the matching *)
RunT3Matching[LUV_, eftOrder_Integer : 5, loopOrder_Integer : 1] := Module[
  {metadata},

  metadata = <|"EFTOrder" -> eftOrder, "LoopOrder" -> loopOrder|>;

  Print["\nRunning Matchete matching..."];
  Print["EFT order: ", eftOrder, "; loop order: ", loopOrder];

  With[{result = RunStages[LUV, MatchingStages[eftOrder, loopOrder], metadata]},
    If[Lookup[result, "Status", ""] === "Success",
      Print["Matching pipeline completed successfully."]
    ];
    result
  ]
];



(* ---------------------------------------------------------------------- *)
(* Intermediate-EFT diagnostics                                            *)
(* ---------------------------------------------------------------------- *)

ExpressionDiagnostics[expr_] := Module[
  {expanded, terms, held, byteCount, leafCount, fields, couplings},

  expanded = Quiet @ Check[Expand[expr], expr];
  terms = If[Head[expanded] === Plus, List @@ expanded, {expanded}];
  held = HoldComplete[expr];

  fields = DeleteDuplicates @ Cases[
    held,
    Field[name_, ___] :> ToString[Unevaluated[name], InputForm],
    Infinity
  ];

  couplings = DeleteDuplicates @ Cases[
    held,
    Coupling[name_, ___] :> ToString[Unevaluated[name], InputForm],
    Infinity
  ];

  <|
    "Head" -> ToString[Head[expr], InputForm],
    "TermCount" -> Length[terms],
    "LeafCount" -> Quiet @ Check[LeafCount[expr], Missing["Unavailable"]],
    "ByteCount" -> Quiet @ Check[ByteCount[expr], Missing["Unavailable"]],
    "Fields" -> fields,
    "Couplings" -> couplings
  |>
];

ExportExpressionDiagnostics[
  outputDirectory_String,
  modelLabel_String,
  level_Integer,
  result_Association
] := Module[
  {keys, payload, path},

  keys = {"RawEFT", "GreenEFT", "EOMEFT", "LoopEFT", "MatchedEFT"};

  payload = Association @ Map[
    Function[key,
      key -> If[
        KeyExistsQ[result, key],
        ExpressionDiagnostics[result[key]],
        <|"Missing" -> True|>
      ]
    ],
    keys
  ];

  path = FileNameJoin @ {
    outputDirectory,
    "debug",
    "sequential_stage_" <> ToString[level] <> "_expressions.json"
  };

  Quiet @ CreateDirectory[DirectoryName[path], CreateIntermediateDirectories -> True];
  Export[path, payload, "RawJSON"];
  path
];


(* ---------------------------------------------------------------------- *)
(* Intermediate-EFT export for later RGE construction                      *)
(* ---------------------------------------------------------------------- *)

ExportIntermediateEFTForRGE[
  outputDirectory_String,
  level_Integer,
  integratedFields_List,
  activeFields_List,
  treeEFT_,
  oneLoopCorrection_,
  compactEFT_
] := Module[
  {
    dataDirectory,
    prefix,
    treePath,
    loopPath,
    compactPath,
    metadataPath,
    wilsonSeedPath,
    quarticSeedPath,
    payload,
    expandedTree,
    treeTerms,
    wilsonTerms,
    scalarQuarticTerms,
    normalizeCGName,
    relevantCGNames,
    quarticCGNames,
    cgRegistry,
    quarticCGRegistry,
    wilsonSeed,
    quarticSeed
  },

  dataDirectory = FileNameJoin[{outputDirectory, "data"}];
  Quiet @ CreateDirectory[
    dataDirectory,
    CreateIntermediateDirectories -> True
  ];

  prefix = "eft" <> ToString[level] <> "_after_" <>
    StringRiffle[integratedFields, "_"];

  treePath = FileNameJoin[{
    dataDirectory,
    prefix <> "_tree_inputform.txt"
  }];

  loopPath = FileNameJoin[{
    dataDirectory,
    prefix <> "_one_loop_inputform.txt"
  }];

  compactPath = FileNameJoin[{
    dataDirectory,
    prefix <> "_compact_inputform.txt"
  }];

  metadataPath = FileNameJoin[{
    dataDirectory,
    prefix <> "_rge_export.json"
  }];

  wilsonSeedPath = FileNameJoin[{
    dataDirectory,
    prefix <> "_wilson_seed.json"
  }];

  quarticSeedPath = FileNameJoin[{
    dataDirectory,
    prefix <> "_scalar_quartic_seed.json"
  }];

  (* InputForm is intentionally used here rather than TeX.  The next RGE
     stage needs the exact Matchete symbolic structure: Field, CG, NCM,
     projectors, hbar, effective couplings, and all gauge-index data. *)
  Export[
    treePath,
    ToString[InputForm[treeEFT]],
    "Text"
  ];

  Export[
    loopPath,
    ToString[InputForm[oneLoopCorrection]],
    "Text"
  ];

  Export[
    compactPath,
    ToString[InputForm[compactEFT]],
    "Text"
  ];

  (* Export the exact tree-level Majorana-type psi^2 phi^2 sector together
     with every CG tensor it references.  This is the bridge to the Python
     MasterWeinbergRGE tensor representation: Python should not infer SU(2)
     contractions from field dimensions alone. *)
  expandedTree = Expand[treeEFT];
  treeTerms = If[
    Head[expandedTree] === Plus,
    List @@ expandedTree,
    {expandedTree}
  ];

  (* Use the serialized InputForm for classification rather than matching
     Matchete's private-context symbols directly.  Objects such as Field,
     GammaCC and Proj can live in Matchete private contexts, so literal
     Mathematica patterns may silently miss valid operators. *)
  wilsonTerms = Select[
    treeTerms,
    Function[term,
      Module[{termText, scalarCount},
        termText = ToString[InputForm[term]];
        scalarCount =
          StringCount[termText, "Field[NewScalar1"] +
          StringCount[termText, "Field[NewScalar2"];

        StringContainsQ[termText, "GammaCC"] &&
        scalarCount === 2
      ]
    ]
  ];

  normalizeCGName[name_] := Replace[name, Bar[x_] :> x];

  relevantCGNames = DeleteDuplicates @ Cases[
    HoldComplete @@ wilsonTerms,
    CG[name_, ___] :> normalizeCGName[name],
    Infinity
  ];

  cgRegistry = DeleteCases[
    Map[
      Function[name,
        Module[{reps, tensor},
          reps = Quiet @ Check[
            Matchete`PackageScope`$CGproperties[name, Indices],
            $Failed
          ];
          tensor = Quiet @ Check[
            Matchete`CGManipulations`PackagePrivate`$CGtensors[name],
            $Failed
          ];

          If[
            reps === $Failed || tensor === $Failed ||
              MissingQ[reps] || MissingQ[tensor],
            Nothing,
            <|
              "Name" -> ToString[Unevaluated[name], InputForm],
              "RepsInputForm" -> ToString[InputForm[reps]],
              "TensorInputForm" -> ToString[InputForm[tensor]]
            |>
          ]
        ]
      ],
      relevantCGNames
    ],
    Nothing
  ];

  wilsonSeed = <|
    "Level" -> level,
    "IntegratedFields" -> integratedFields,
    "ActiveHeavyFields" -> activeFields,
    "TreeWilsonTermCount" -> Length[wilsonTerms],
    "TreeWilsonTerms" -> MapIndexed[
      Function[{term, position},
        <|
          "Index" -> First[position],
          "Chirality" -> Module[{termText = ToString[InputForm[term]]},
            Which[
              StringContainsQ[termText, "Proj[-1]"], "PL",
              StringContainsQ[termText, "Proj[1]"], "PR",
              True, "Unknown"
            ]
          ],
          "Scalar1Count" ->
            StringCount[ToString[InputForm[term]], "Field[NewScalar1"],
          "Scalar2Count" ->
            StringCount[ToString[InputForm[term]], "Field[NewScalar2"],
          "BarredScalar1Count" ->
            StringCount[
              ToString[InputForm[term]],
              "Bar[Field[NewScalar1"
            ],
          "BarredScalar2Count" ->
            StringCount[
              ToString[InputForm[term]],
              "Bar[Field[NewScalar2"
            ],
          "CGNames" -> (
            ToString[Unevaluated[#], InputForm] & /@
              DeleteDuplicates @ Cases[
                HoldComplete[term],
                CG[name_, ___] :> normalizeCGName[name],
                Infinity
              ]
          ),
          "TermInputForm" -> ToString[InputForm[term]]
        |>
      ],
      wilsonTerms
    ],
    "CGRegistry" -> cgRegistry
  |>;

  Export[wilsonSeedPath, wilsonSeed, "RawJSON"];

  Print[
    "  Intermediate EFT Wilson seed: ",
    Length[wilsonTerms],
    " tree-level psi^2 phi^2 term(s), ",
    Length[cgRegistry],
    " CG definition(s) -> ",
    wilsonSeedPath
  ];

  (* Export the complete renormalisable scalar-quartic sector in exact
     InputForm.  The Python RGE layer will convert these complex-field
     monomials into the fully symmetric real-scalar lambda_abcd tensor using

         V4 = (1/4!) lambda_abcd phi_a phi_b phi_c phi_d.

     Select structurally from the canonical tree EFT: exactly four scalar
     field factors, no fermion fields, and no derivatives.  This keeps all
     Higgs, S1, S2, portal, cross-scalar, and T3 loop-closing quartics without
     hard-coding a representation-dependent list of couplings. *)
  scalarQuarticTerms = Select[
    treeTerms,
    Function[term,
      Module[{termText, scalarCount},
        termText = ToString[InputForm[term]];
        scalarCount =
          StringCount[termText, "Field[H, Scalar"] +
          StringCount[termText, "Field[NewScalar1, Scalar"] +
          StringCount[termText, "Field[NewScalar2, Scalar"];

        scalarCount === 4 &&
        !StringContainsQ[termText, ", Fermion,"] &&
        !StringContainsQ[termText, "FieldStrength["] &&
        !StringContainsQ[termText, "GammaM["] &&
        !StringContainsQ[termText, "CovariantD["]
      ]
    ]
  ];

  quarticCGNames = DeleteDuplicates @ Cases[
    HoldComplete @@ scalarQuarticTerms,
    CG[name_, ___] :> normalizeCGName[name],
    Infinity
  ];

  quarticCGRegistry = DeleteCases[
    Map[
      Function[name,
        Module[{reps, tensor},
          reps = Quiet @ Check[
            Matchete`PackageScope`$CGproperties[name, Indices],
            $Failed
          ];
          tensor = Quiet @ Check[
            Matchete`CGManipulations`PackagePrivate`$CGtensors[name],
            $Failed
          ];

          If[
            reps === $Failed || tensor === $Failed ||
              MissingQ[reps] || MissingQ[tensor],
            Nothing,
            <|
              "Name" -> ToString[Unevaluated[name], InputForm],
              "RepsInputForm" -> ToString[InputForm[reps]],
              "TensorInputForm" -> ToString[InputForm[tensor]]
            |>
          ]
        ]
      ],
      quarticCGNames
    ],
    Nothing
  ];

  quarticSeed = <|
    "Level" -> level,
    "IntegratedFields" -> integratedFields,
    "ActiveHeavyFields" -> activeFields,
    "ScalarQuarticTermCount" -> Length[scalarQuarticTerms],
    "ScalarQuarticTerms" -> MapIndexed[
      Function[{term, position},
        Module[{termText, couplingNames},
          termText = ToString[InputForm[term]];
          couplingNames = DeleteDuplicates @ Cases[
            HoldComplete[term],
            coupling_[___] /;
              Head[Unevaluated[coupling]] === Symbol &&
              StringMatchQ[
                SymbolName[Unevaluated[coupling]],
                "lambda" ~~ ___ | "\[Lambda]"
              ] :>
                SymbolName[Unevaluated[coupling]],
            Infinity
          ];

          <|
            "Index" -> First[position],
            "Couplings" -> couplingNames,
            "CGNames" -> (
              ToString[Unevaluated[#], InputForm] & /@
                DeleteDuplicates @ Cases[
                  HoldComplete[term],
                  CG[name_, ___] :> normalizeCGName[name],
                  Infinity
                ]
            ),
            "HiggsFieldCount" ->
              StringCount[termText, "Field[H, Scalar"],
            "Scalar1FieldCount" ->
              StringCount[termText, "Field[NewScalar1, Scalar"],
            "Scalar2FieldCount" ->
              StringCount[termText, "Field[NewScalar2, Scalar"],
            "TermInputForm" -> termText
          |>
        ]
      ],
      scalarQuarticTerms
    ],
    "CGRegistry" -> quarticCGRegistry
  |>;

  Export[quarticSeedPath, quarticSeed, "RawJSON"];

  Print[
    "  Intermediate EFT scalar-quartic seed: ",
    Length[scalarQuarticTerms],
    " tree-level scalar quartic term(s), ",
    Length[quarticCGRegistry],
    " CG definition(s) -> ",
    quarticSeedPath
  ];

  payload = <|
    "Level" -> level,
    "IntegratedFields" -> integratedFields,
    "ActiveHeavyFields" -> activeFields,
    "TreeFile" -> FileNameTake[treePath],
    "OneLoopFile" -> FileNameTake[loopPath],
    "CompactFile" -> FileNameTake[compactPath],
    "WilsonSeedFile" -> FileNameTake[wilsonSeedPath],
    "ScalarQuarticSeedFile" -> FileNameTake[quarticSeedPath],
    "TreeDiagnostics" -> ExpressionDiagnostics[treeEFT],
    "OneLoopDiagnostics" -> ExpressionDiagnostics[oneLoopCorrection],
    "CompactDiagnostics" -> ExpressionDiagnostics[compactEFT]
  |>;

  Export[metadataPath, payload, "RawJSON"];

  Print[
    "  Intermediate EFT RGE export: ",
    metadataPath
  ];

  <|
    "TreePath" -> treePath,
    "OneLoopPath" -> loopPath,
    "CompactPath" -> compactPath,
    "MetadataPath" -> metadataPath,
    "WilsonSeedPath" -> wilsonSeedPath,
    "ScalarQuarticSeedPath" -> quarticSeedPath
  |>
];


(* ---------------------------------------------------------------------- *)
(* Fixed-order sequential matching helpers                                 *)
(* ---------------------------------------------------------------------- *)

(* A private bookkeeping symbol used only to propagate an already-generated
   one-loop correction through the tree-level matching at the next threshold.
   We always extract only the coefficient linear in this symbol. *)
ClearAll[T3SequentialLoopBookkeeper];

EnsureSequentialLoopBookkeeper[] := Module[{result},
  result = Quiet @ Check[
    DefineCoupling[
      T3SequentialLoopBookkeeper,
      SelfConjugate -> True
    ],
    $Failed
  ];

  Print[
    "  Sequential loop bookkeeper registration: ",
    If[result === $Failed, "FAILED", "OK"]
  ];

  If[result === $Failed, $Failed, True]
];

PrepareSequentialMatchInput[lag_] := Module[{green, eom},
  Print["  Canonicalising EFT before the next threshold..."];

  green = SafeStage[GreensSimplify[lag]];
  If[MemberQ[{$Failed, $Aborted}, green], Return[$Failed]];

  eom = SafeStage[EOMSimplify[green]];
  If[MemberQ[{$Failed, $Aborted}, eom], Return[$Failed]];

  eom
];

CompactExplicitEFT[result_Association] := Module[{eom, explicit},
  eom = Lookup[result, "EOMEFT", $Failed];
  If[eom === $Failed, Return[$Failed]];

  (* Replace the effective couplings so loop-order bookkeeping remains
     explicit, but deliberately do NOT evaluate loop functions here. *)
  explicit = SafeStage[ReplaceEffectiveCouplings[eom]];
  explicit
];

CanonicalTransitionEFT[result_Association] := Module[{eom},
  eom = Lookup[result, "EOMEFT", $Failed];
  If[eom === $Failed, Return[$Failed]];

  (* IMPORTANT:
     This is the representation carried to a later threshold.
     Do not call ReplaceEffectiveCouplings here.  EOMSimplify has already
     performed the field redefinitions needed to keep the kinetic terms
     canonical.  Expanding the effective couplings can re-expose the
     wave-function / renormalizable threshold shifts and make the next
     Match reject the Lagrangian as non-canonical. *)
  eom
];


LinearBookkeepingCoefficient[expr_, n_Integer] := Module[{expanded},
  expanded = Quiet @ Check[Expand[expr], expr];
  Coefficient[expanded, T3SequentialLoopBookkeeper[], n]
];

(* ---------------------------------------------------------------------- *)
(* Sequential threshold matching                                           *)
(* ---------------------------------------------------------------------- *)


(* ---------------------------------------------------------------------- *)
(* Fresh-kernel threshold transition                                      *)
(* ---------------------------------------------------------------------- *)


PrintFreshKernelSummary[stdout_String, stderr_String] := Module[
  {lines, important, compact},

  lines = Select[
    StringSplit[stdout <> "\n" <> stderr, {"\r\n", "\n", "\r"}],
    StringLength[StringTrim[#]] > 0 &
  ];

  important = Select[
    lines,
    Function[line,
      Or[
        StringStartsQ[StringTrim[line], "[fresh kernel]"],
        StringStartsQ[StringTrim[line], "ERROR:"],
        StringStartsQ[StringTrim[line], "WARNING:"],
        StringStartsQ[StringTrim[line], "Running Matchete matching"],
        StringStartsQ[StringTrim[line], "EFT order:"],
        StringStartsQ[StringTrim[line], "Running Match"],
        StringContainsQ[line, "Match result head:"],
        StringContainsQ[line, "Matching pipeline completed successfully"],
        StringContainsQ[line, "CheckLagrangian::"],
        StringContainsQ[line, "FindDummyIndices::"],
        StringContainsQ[line, "UndefinedObject"],
        StringTrim[line] === "$Aborted"
      ]
    ]
  ];

  compact = DeleteDuplicates @ Map[
    Function[line,
      Which[
        StringContainsQ[line, "CheckLagrangian::Hermiticity"],
          "  [fresh kernel] Matchete error: CheckLagrangian::Hermiticity",
        StringContainsQ[line, "CheckLagrangian::CanonicallyNormalized"],
          "  [fresh kernel] Matchete error: CheckLagrangian::CanonicallyNormalized",
        StringContainsQ[line, "FindDummyIndices::trippleindex"],
          "  [fresh kernel] Matchete error: FindDummyIndices::trippleindex",
        StringLength[line] > 220,
          StringTake[line, 220] <> " ... [full message in stdout.log]",
        True,
          line
      ]
    ],
    important
  ];

  Scan[Print, compact];

  If[Length[compact] < Length[lines],
    Print[
      "  [fresh kernel] Suppressed ",
      Length[lines] - Length[compact],
      " verbose line(s); full child output is saved to stdout.log."
    ]
  ];
];

RunFreshThresholdStage[
  currentTree_,
  currentLoop_,
  transitionTree_,
  transitionFull_,
  model_Association,
  group_List,
  active_List,
  level_Integer,
  eftOrder_Integer,
  loopOrder_Integer,
  outputDirectory_String
] := Module[
  {
    stageScript, transitionDir, treePath, loopPath,
    transitionTreePath, transitionFullPath, cgPath,
    resultPath, logPath, command, process, result, dS1, dS2, dF, alpha,
    cgNames, cgRegistry, normalizeCGName
  },

  stageScript = FileNameJoin @ {
    DirectoryName @ ExpandFileName[$InputFileName],
    "RunThresholdStage.wl"
  };

  transitionDir = FileNameJoin @ {
    outputDirectory,
    "debug",
    "fresh_kernel_stage_" <> ToString[level]
  };
  Quiet @ CreateDirectory[
    transitionDir,
    CreateIntermediateDirectories -> True
  ];

  treePath = FileNameJoin[{transitionDir, "input_tree.wxf"}];
  loopPath = FileNameJoin[{transitionDir, "input_loop.wxf"}];
  transitionTreePath =
    FileNameJoin[{transitionDir, "input_transition_tree.wxf"}];
  transitionFullPath =
    FileNameJoin[{transitionDir, "input_transition_full.wxf"}];
  cgPath = FileNameJoin[{transitionDir, "input_cg_registry.wxf"}];
  resultPath = FileNameJoin[{transitionDir, "result.wxf"}];
  logPath = FileNameJoin[{transitionDir, "stdout.log"}];

  (* The explicit tree/loop split is retained for reporting and diagnostics.
     The separate transition pair is what the next Match actually consumes:
       transitionTree = canonical EOM EFT at O(hbar^0)
       transitionFull = canonical EOM EFT through O(hbar)
     Crucially, effective couplings remain unexpanded in these transition
     expressions. *)
  Export[treePath, currentTree, "WXF"];
  Export[loopPath, currentLoop, "WXF"];
  Export[transitionTreePath, transitionTree, "WXF"];
  Export[transitionFullPath, transitionFull, "WXF"];

  Print[
    "  Exporting canonical transition EFTs without expanding effective ",
    "couplings."
  ];

  normalizeCGName[name_] := Replace[name, Bar[x_] :> x];

  cgNames = DeleteDuplicates @ Cases[
    HoldComplete[
      currentTree,
      currentLoop,
      transitionTree,
      transitionFull
    ],
    CG[name_, ___] :> normalizeCGName[name],
    Infinity
  ];

  cgRegistry = DeleteCases[
    Map[
      Function[name,
        Module[{reps, tensor},
          reps = Quiet @ Check[
            Matchete`PackageScope`$CGproperties[name, Indices],
            $Failed
          ];
          tensor = Quiet @ Check[
            Matchete`CGManipulations`PackagePrivate`$CGtensors[name],
            $Failed
          ];
          If[
            reps === $Failed || tensor === $Failed ||
              MissingQ[reps] || MissingQ[tensor],
            Nothing,
            <|"Name" -> name, "Reps" -> reps, "Tensor" -> tensor|>
          ]
        ]
      ],
      cgNames
    ],
    Nothing
  ];

  Export[cgPath, cgRegistry, "WXF"];

  Print[
    "  Exported fresh-kernel CG registry: ",
    Length[cgRegistry],
    " definition(s)."
  ];

  dS1 = model["Scalar1", "SU2"];
  dS2 = model["Scalar2", "SU2"];
  dF = model["Fermion", "SU2"];
  alpha = model["Alpha"];

  Print[
    "\n  Launching fresh Wolfram kernel for threshold stage ",
    level,
    " (", StringRiffle[group, ", "], ")..."
  ];

  command = {
    "wolframscript",
    "-file",
    stageScript,
    resultPath,
    ToString[eftOrder],
    ToString[loopOrder],
    ToString[dS1],
    ToString[dS2],
    ToString[dF],
    If[alpha < 0, "m" <> ToString[Abs[alpha]], "p" <> ToString[alpha]],
    StringRiffle[active, ","],
    StringRiffle[group, ","],
    treePath,
    loopPath,
    transitionTreePath,
    transitionFullPath,
    cgPath
  };

  Print[
    "  Fresh-kernel CLI payload count (after -file script): ",
    Length[command] - 3
  ];

  process = RunProcess[command];

  Export[
    logPath,
    Lookup[process, "StandardOutput", ""] <>
      Lookup[process, "StandardError", ""],
    "Text"
  ];

  PrintFreshKernelSummary[
    Lookup[process, "StandardOutput", ""],
    Lookup[process, "StandardError", ""]
  ];

  If[
    Lookup[process, "ExitCode", 1] =!= 0 ||
      !FileExistsQ[resultPath],
    Print[
      "ERROR: fresh-kernel threshold stage ", level,
      " failed. See ", logPath
    ];
    Return[$Failed]
  ];

  result = Quiet @ Check[Import[resultPath, "WXF"], $Failed];

  If[
    !AssociationQ[result] ||
      Lookup[result, "Status", ""] =!= "Success",
    Print[
      "ERROR: fresh-kernel threshold stage ", level,
      " returned an invalid result. See ", logPath
    ];
    Return[$Failed]
  ];

  result
];

RunSequentialT3Matching[
  LUV_,
  model_Association,
  thresholdPlan_List,
  eftOrder_Integer : 5,
  loopOrder_Integer : 1,
  outputDirectory_String : "",
  modelLabel_String : "T3"
] := Module[
  {
    stageResults = {},
    active = {"F", "S1", "S2"},
    currentTree = LUV,
    currentLoop = 0,
    currentTransitionTree = LUV,
    currentTransitionFull = LUV,
    group,
    level,
    treeInput,
    treeResult,
    oneLoopResult,
    treeExplicit,
    fullFromTreeExplicit,
    localLoopCorrection,
    inheritedLoopCorrection,
    probeInput,
    probeCanonical,
    probeTreeResult,
    probeExplicit,
    probeTreeCheck,
    nextTree,
    nextLoop,
    nextFullCompact,
    nextFullReport,
    nextTransitionTree,
    nextTransitionFull,
    diagnosticPath,
    freshStage
  },

  If[thresholdPlan === {},
    Return[<|"Status" -> "EmptyThresholdPlan", "Stages" -> {}|>]
  ];

  Do[
    group = thresholdPlan[[level]];

    Print[
      "\n=== Threshold stage ", level, ": integrating ",
      StringRiffle[group, ", "], " ==="
    ];

    (* Stage 1 is matched in the UV kernel that built the model.  Every later
       threshold is deliberately delegated to a brand-new Wolfram kernel.
       This prevents Matchete's global field registry from carrying the
       previous stage's Light/Heavy assignment into the next EFT. *)
    If[level > 1,
      freshStage = RunFreshThresholdStage[
        currentTree,
        currentLoop,
        currentTransitionTree,
        currentTransitionFull,
        model,
        group,
        active,
        level,
        eftOrder,
        loopOrder,
        outputDirectory
      ];

      If[freshStage === $Failed,
        Return[<|
          "Status" -> "FreshKernelThresholdFailed",
          "FailedLevel" -> level,
          "Stages" -> stageResults
        |>]
      ];

      active = Complement[active, group];

      AppendTo[
        stageResults,
        <|
          "Level" -> level,
          "IntegratedFields" -> group,
          "ActiveHeavyFields" -> active,
          "TreeLagrangian" -> freshStage["TreeLagrangian"],
          "OneLoopCorrection" -> freshStage["OneLoopCorrection"],
          "CompactLagrangian" -> freshStage["CompactLagrangian"],
          "Lagrangian" -> freshStage["Lagrangian"]
        |>
      ];

      currentTree = freshStage["TreeLagrangian"];
      currentLoop = freshStage["OneLoopCorrection"];
      currentTransitionTree = freshStage["TransitionTreeLagrangian"];
      currentTransitionFull = freshStage["TransitionFullLagrangian"];

      Continue[];
    ];

    (* -------------------------------------------------------------- *)
    (* Prepare the incoming EFT while the PREVIOUS stage registry is   *)
    (* still active.  At this point the fields that survive the last   *)
    (* threshold are still flagged Light, so EOMSimplify is allowed to *)
    (* perform the field redefinitions needed for canonical form.      *)
    (* -------------------------------------------------------------- *)

    treeInput = If[
      level === 1,
      currentTree,
      PrepareSequentialMatchInput[currentTree]
    ];

    If[treeInput === $Failed,
      Return[<|
        "Status" -> "TreeInputCanonicalisationFailed",
        "FailedLevel" -> level,
        "Stages" -> stageResults
      |>]
    ];

    (* Prepare the full EFT through one loop BEFORE changing the surviving
       fields from Light to Heavy.

       Matchete already carries loop order with its native symbol hbar
       (hbar = 1/(16 Pi^2)).  Using the actual one-loop EFT here is preferable
       to introducing an artificial epsilon: currentTree + currentLoop is a
       canonical EFT through O(hbar).  After the tree-level match at the next
       threshold we explicitly project the result back to O(hbar), removing
       all hbar^2 and higher terms generated by nonlinear tree matching. *)
    probeCanonical = $Failed;

    If[level > 1 && currentLoop =!= 0,
      Print[
        "  Preparing full inherited EFT through O(hbar) while surviving ",
        "fields are still Light..."
      ];

      probeInput = Expand[currentTree + currentLoop];

      If[FreeQ[probeInput, hbar],
        Print[
          "ERROR: inherited one-loop EFT contains no explicit Matchete hbar. ",
          "Cannot perform fixed-order sequential projection safely."
        ];
        Return[<|
          "Status" -> "MissingHbarBookkeeping",
          "FailedLevel" -> level,
          "Stages" -> stageResults
        |>]
      ];

      probeCanonical = PrepareSequentialMatchInput[probeInput];

      If[probeCanonical === $Failed,
        Return[<|
          "Status" -> "LoopProbeCanonicalisationFailed",
          "FailedLevel" -> level,
          "Stages" -> stageResults
        |>]
      ]
    ];

    If[level > 1,
      Print[
        "  Re-registering active T3 fields: ",
        active,
        "; heavy at this stage: ",
        group
      ];

      If[ReclassifyT3Fields[model, group, active] === $Failed,
        Return[<|
          "Status" -> "FieldReclassificationFailed",
          "FailedLevel" -> level,
          "Stages" -> stageResults
        |>]
      ];

      Print[
        "  Preserved integrated mass parameters as real couplings: ",
        Complement[{"F", "S1", "S2"}, active]
      ];

    ];

    (* -------------------------------------------------------------- *)
    (* 1. Tree-level image of the tree-level EFT                      *)
    (* -------------------------------------------------------------- *)

    Print["\n  [A] Tree matching of the tree-level EFT"];
    treeResult = RunT3Matching[treeInput, eftOrder, 0];

    If[
      !AssociationQ[treeResult] ||
        Lookup[treeResult, "Status", ""] =!= "Success",
      Return[<|
        "Status" -> "TreeMatchingFailed",
        "FailedLevel" -> level,
        "StageResult" -> treeResult,
        "Stages" -> stageResults
      |>]
    ];

    treeExplicit = CompactExplicitEFT[treeResult];
    nextTransitionTree = CanonicalTransitionEFT[treeResult];

    If[
      treeExplicit === $Failed || nextTransitionTree === $Failed,
      Return[<|
        "Status" -> "TreeExplicitEFTFailed",
        "FailedLevel" -> level,
        "Stages" -> stageResults
      |>]
    ];

    (* -------------------------------------------------------------- *)
    (* 2. One-loop matching acts ONLY on the tree-level EFT           *)
    (* -------------------------------------------------------------- *)

    If[loopOrder >= 1,
      Print["\n  [B] One-loop matching of the tree-level EFT"];
      oneLoopResult = RunT3Matching[treeInput, eftOrder, 1];

      If[
        !AssociationQ[oneLoopResult] ||
          Lookup[oneLoopResult, "Status", ""] =!= "Success",
        Return[<|
          "Status" -> "OneLoopMatchingFailed",
          "FailedLevel" -> level,
          "StageResult" -> oneLoopResult,
          "Stages" -> stageResults
        |>]
      ];

      fullFromTreeExplicit = CompactExplicitEFT[oneLoopResult];
      nextTransitionFull = CanonicalTransitionEFT[oneLoopResult];

      If[
        fullFromTreeExplicit === $Failed ||
          nextTransitionFull === $Failed,
        Return[<|
          "Status" -> "OneLoopExplicitEFTFailed",
          "FailedLevel" -> level,
          "Stages" -> stageResults
        |>]
      ];

      (* M^(1)[L^(0)] = M_{tree+1loop}[L^(0)] - M_tree[L^(0)] *)
      localLoopCorrection = Expand[
        fullFromTreeExplicit - treeExplicit
      ],
      localLoopCorrection = 0;
      oneLoopResult = treeResult;
      nextTransitionFull = nextTransitionTree
    ];

    (* -------------------------------------------------------------- *)
    (* 3. Propagate the loop correction inherited from earlier stages *)
    (*    through TREE matching only.                                 *)
    (*                                                                *)
    (*    M0[L0 + eps L1] = M0[L0] + eps dM0[L1] + O(eps^2).         *)
    (*    We explicitly extract the eps coefficient, so no two-loop   *)
    (*    contamination from (L1)^2 is retained.                      *)
    (* -------------------------------------------------------------- *)

    inheritedLoopCorrection = 0;

    If[level > 1 && currentLoop =!= 0,
      Print[
        "\n  [C] Tree propagation of the inherited one-loop EFT "
        "(native hbar projection)"
      ];

      (* probeCanonical was prepared above while the surviving fields were
         still Light.  Perform only a tree-level Match here. *)
      probeTreeResult = RunT3Matching[probeCanonical, eftOrder, 0];

      If[
        !AssociationQ[probeTreeResult] ||
          Lookup[probeTreeResult, "Status", ""] =!= "Success",
        Return[<|
          "Status" -> "InheritedLoopTreeMatchingFailed",
          "FailedLevel" -> level,
          "StageResult" -> probeTreeResult,
          "Stages" -> stageResults
        |>]
      ];

      probeExplicit = CompactExplicitEFT[probeTreeResult];

      If[probeExplicit === $Failed,
        Return[<|
          "Status" -> "InheritedLoopExplicitEFTFailed",
          "FailedLevel" -> level,
          "Stages" -> stageResults
        |>]
      ];

      (* M0[L0 + hbar L1] - M0[L0] can contain hbar^2, hbar^3, ... because
         tree matching is nonlinear in EFT couplings.  Keep exactly the
         coefficient linear in Matchete's native loop-counting parameter. *)
      probeTreeCheck = Expand[probeExplicit /. hbar -> 0];

      inheritedLoopCorrection = Expand[
        hbar Coefficient[
          Expand[probeExplicit - treeExplicit],
          hbar,
          1
        ]
      ];

      Print[
        "    inherited one-loop terms: ",
        Lookup[
          ExpressionDiagnostics[inheritedLoopCorrection],
          "TermCount",
          "?"
        ]
      ];

      Print[
        "    hbar^0 tree-check byte difference: ",
        Quiet @ Check[
          ByteCount @ Expand[probeTreeCheck - treeExplicit],
          "Unavailable"
        ]
      ];

      Print[
        "    discarded higher-loop terms present: ",
        !FreeQ[
          Expand[probeExplicit - treeExplicit -
            inheritedLoopCorrection],
          hbar
        ]
      ];
    ];

    (* -------------------------------------------------------------- *)
    (* 4. Assemble the EFT truncated consistently through one loop.   *)
    (* -------------------------------------------------------------- *)

    nextTree = treeExplicit;
    nextLoop = Expand[
      localLoopCorrection + inheritedLoopCorrection
    ];
    nextFullCompact = Expand[nextTree + nextLoop];

    (* Evaluate loop functions only for the report/final exported form.
       The compact tree/loop split is what is carried to another threshold. *)
    nextFullReport = SafeStage[
      EvaluateLoopFunctions[nextFullCompact]
    ];

    If[MemberQ[{$Failed, $Aborted}, nextFullReport],
      nextFullReport = nextFullCompact
    ];

    active = Complement[active, group];

    AppendTo[
      stageResults,
      <|
        "Level" -> level,
        "IntegratedFields" -> group,
        "ActiveHeavyFields" -> active,
        "TreeMatching" -> treeResult,
        "OneLoopMatching" -> oneLoopResult,
        "TreeLagrangian" -> nextTree,
        "OneLoopCorrection" -> nextLoop,
        "TransitionTreeLagrangian" -> nextTransitionTree,
        "TransitionFullLagrangian" -> nextTransitionFull,
        "CompactLagrangian" -> nextFullCompact,
        "Lagrangian" -> nextFullReport
      |>
    ];

    If[outputDirectory =!= "",
      diagnosticPath = FileNameJoin @ {
        outputDirectory,
        "debug",
        "perturbative_stage_" <> ToString[level] <> ".json"
      };

      Quiet @ CreateDirectory[
        DirectoryName[diagnosticPath],
        CreateIntermediateDirectories -> True
      ];

      Export[
        diagnosticPath,
        <|
          "Level" -> level,
          "IntegratedFields" -> group,
          "ActiveHeavyFields" -> active,
          "Tree" -> ExpressionDiagnostics[nextTree],
          "LocalOneLoop" -> ExpressionDiagnostics[localLoopCorrection],
          "InheritedOneLoop" -> ExpressionDiagnostics[inheritedLoopCorrection],
          "TotalOneLoop" -> ExpressionDiagnostics[nextLoop],
          "FullCompact" -> ExpressionDiagnostics[nextFullCompact]
        |>,
        "RawJSON"
      ];

      Print["  Perturbative-stage diagnostics: ", diagnosticPath];

      ExportIntermediateEFTForRGE[
        outputDirectory,
        level,
        group,
        active,
        nextTree,
        nextLoop,
        nextFullCompact
      ]
    ];

    If[level < Length[thresholdPlan],
      Print[
        "  Carrying separate tree and one-loop EFT pieces to threshold stage ",
        level + 1,
        "."
      ]
    ];

    currentTree = nextTree;
    currentLoop = nextLoop;
    currentTransitionTree = nextTransitionTree;
    currentTransitionFull = nextTransitionFull,
    {level, Length[thresholdPlan]}
  ];

  <|
    "Status" -> "Success",
    "Stages" -> stageResults,
    "MatchedEFT" -> Last[stageResults]["Lagrangian"],
    "CompactMatchedEFT" -> Last[stageResults]["CompactLagrangian"],
    "TreeEFT" -> Last[stageResults]["TreeLagrangian"],
    "OneLoopCorrection" -> Last[stageResults]["OneLoopCorrection"],
    "EFTOrder" -> eftOrder,
    "LoopOrder" -> loopOrder,
    "PerturbativeTruncation" -> "O(1-loop)"
  |>
];

(* Pure SM has no heavy field: if Match fails, canonicalise LSM directly and use standard forms for it. 
 We are comparing SM and full lagrangian so we match Pure SM as well*)
RunSMBaselineMatching[LSM_, eftOrder_Integer : 5, loopOrder_Integer : 1] := Module[
  {attempt, metadata, stages, result},

  attempt = RunT3Matching[LSM, eftOrder, loopOrder];

  If[AssociationQ[attempt] && Lookup[attempt, "Status", ""] === "Success",
    Return[Append[attempt, "BaselineMode" -> "FullMatchPipeline"]]
  ];

  If[!(AssociationQ[attempt] && Lookup[attempt, "Status", ""] === "MatchFailed"),
    Return[attempt]
  ];

  Print["Pure SM has no matchable heavy field; canonicalising LSM directly."];

  (* Keep the historical failure-status strings used by downstream diagnostics. *)
  stages = {
    {"GreensSimplifySM", GreensSimplify, "GreenEFT", "SMGreensSimplifyFailed"},
    {"EOMSimplifySM", EOMSimplify, "EOMEFT", "SMEOMSimplifyFailed"},
    {
      "EvaluateLoopFunctionsSM",
      EvaluateLoopFunctions,
      "LoopEFT",
      "SMLoopEvaluationFailed"
    },
    {
      "ReplaceEffectiveCouplingsSM",
      ReplaceEffectiveCouplings,
      "MatchedEFT",
      "SMEffectiveCouplingReplacementFailed"
    }
  };

  metadata = <|
    "BaselineMode" -> "DirectCanonicalisation",
    "EFTOrder" -> eftOrder,
    "LoopOrder" -> loopOrder
  |>;

  RunStages[LSM, stages, metadata]
];

(* Gives us BSM contribution to EFT by doing FullEFT-SMEFT=BSMEFT *)
BuildMatchedEFTDifference[fullEFT_, smEFT_] := Module[
  {
    input,
    greenDifference,
    eomDifference,
    canonicalInput,
    loopDifference,
    bsmEFT
  },

  input = Expand[fullEFT - smEFT];

  Print["Running GreensSimplifyDifference..."];
  greenDifference = SafeStage[GreensSimplify[input]];
  If[MemberQ[{$Failed, $Aborted}, greenDifference],
    Print["ERROR: GreensSimplifyDifference failed or aborted."];
    Return[<|
      "Status" -> "DifferenceGreensSimplifyFailed",
      "InputDifference" -> input,
      "GreenDifference" -> greenDifference
    |>]
  ];

  Print["Running EOMSimplifyDifference..."];
  eomDifference = SafeStage[EOMSimplify[greenDifference]];

  canonicalInput = If[
    MemberQ[{$Failed, $Aborted}, eomDifference],
    Print[
      "EOMSimplifyDifference could not canonicalise the difference-only EFT; ",
      "continuing from GreenDifference."
    ];
    greenDifference,
    eomDifference
  ];

  Print["Running EvaluateLoopFunctionsDifference..."];
  loopDifference = SafeStage[EvaluateLoopFunctions[canonicalInput]];
  If[MemberQ[{$Failed, $Aborted}, loopDifference],
    Print["ERROR: EvaluateLoopFunctionsDifference failed or aborted."];
    Return[<|
      "Status" -> "DifferenceLoopEvaluationFailed",
      "InputDifference" -> input,
      "GreenDifference" -> greenDifference,
      "EOMDifference" -> eomDifference,
      "LoopDifference" -> loopDifference
    |>]
  ];

  Print["Running ReplaceEffectiveCouplingsDifference..."];
  bsmEFT = SafeStage[ReplaceEffectiveCouplings[loopDifference]];
  If[MemberQ[{$Failed, $Aborted}, bsmEFT],
    Print["ERROR: ReplaceEffectiveCouplingsDifference failed or aborted."];
    Return[<|
      "Status" -> "DifferenceEffectiveCouplingReplacementFailed",
      "InputDifference" -> input,
      "GreenDifference" -> greenDifference,
      "EOMDifference" -> eomDifference,
      "LoopDifference" -> loopDifference,
      "BSMEFT" -> bsmEFT
    |>]
  ];

  Print["Matched EFT difference canonicalised successfully."];

  <|
    "Status" -> "Success",
    "InputDifference" -> input,
    "GreenDifference" -> greenDifference,
    "EOMDifference" -> eomDifference,
    "EOMFallbackUsed" -> MemberQ[{$Failed, $Aborted}, eomDifference],
    "LoopDifference" -> loopDifference,
    "BSMEFT" -> bsmEFT
  |>
];


(* ---------------------------------------------------------------------- *)
(* Weinberg-operator extraction                                            *)
(* ---------------------------------------------------------------------- *)

ClearAll[
  InternalHeadName,
  InternalSymbolName,
  MatcheteFieldName,
  CountMatcheteField,
  ContainsNamedSymbolQ,
  ContainsProjectorQ,
  BarredMatcheteFieldQ,
  WeinbergLHHTermQ,
  ExtractWeinbergTerms,
  StripWeinbergOperatorStructure,
  CompactWeinbergCoefficient,
  ExtractWeinbergCoefficient
];

InternalHeadName[x_] := Quiet@Check[
  SymbolName[Unevaluated[Head[x]]],
  ToString[Unevaluated[Head[x]], InputForm]
];

InternalSymbolName[x_Symbol] := SymbolName[Unevaluated[x]];
InternalSymbolName[x_] := ToString[Unevaluated[x], InputForm];

MatcheteFieldName[field_] := Module[{args},
  If[InternalHeadName[Unevaluated[field]] =!= "Field", Return[""]];
  args = List @@ Unevaluated[field];
  If[args === {}, "", InternalSymbolName[args[[1]]]]
];

CountMatcheteField[expr_, name_String] := Count[
  Unevaluated[expr],
  object_ /; InternalHeadName[Unevaluated[object]] === "Field" &&
    MatcheteFieldName[Unevaluated[object]] === name,
  Infinity
];

ContainsNamedSymbolQ[expr_, name_String] := !FreeQ[
  Unevaluated[expr],
  symbol_Symbol /; InternalSymbolName[Unevaluated[symbol]] === name,
  Infinity
];

ContainsProjectorQ[expr_, chirality_Integer] := !FreeQ[
  Unevaluated[expr],
  object_ /; InternalHeadName[Unevaluated[object]] === "Proj" &&
    Quiet@Check[(List @@ Unevaluated[object]) === {chirality}, False],
  Infinity
];

BarredMatcheteFieldQ[object_, name_String] := Module[{args, inner},
  If[InternalHeadName[Unevaluated[object]] =!= "Bar", Return[False]];
  args = List @@ Unevaluated[object];
  If[Length[args] =!= 1, Return[False]];

  inner = args[[1]];
  InternalHeadName[Unevaluated[inner]] === "Field" &&
    MatcheteFieldName[Unevaluated[inner]] === name
];

(* GammaCC is the marker for the Weinberg operator. *)
WeinbergLHHTermQ[term_] := !FreeQ[Unevaluated[term], GammaCC];

(*
1. We find all terms with GammaCC and thats it
*)
ExtractWeinbergTerms[eft_] := Module[{expanded, terms},
  expanded = Expand[eft];

  (* Get Weinberg Terms based on the prefactor with Times[] with GammaCC*)
  terms = DeleteDuplicates @ Cases[
    expanded,
    term_Times /; !FreeQ[term, GammaCC],
    Infinity
  ];

  (* If Times[] doesnt exist we look for NCM (non commutative multiplication chain) instead which should have GammaCC *)
  If[terms === {} && !FreeQ[expanded, GammaCC],
    terms = DeleteDuplicates @ Cases[
      expanded,
      chain_NCM /; !FreeQ[chain, GammaCC],
      Infinity
    ]
  ];

  terms
];

(* Strip only the universal LLHH structure; retain the physical prefactor. *)
StripWeinbergOperatorStructure[term_] := Module[{stripped},
  stripped = term /. {
    chain_NCM /; !FreeQ[chain, GammaCC] :> 1,
    Field[H, Scalar, inds_, derivs_] :> 1,
    CG[___] :> 1
  };

  Quiet@Check[Simplify[Expand[stripped]], stripped]
];

(* 
Combine all non conjugated term prefactors into a compact C5 expression.
So A*LLHH + B*LLHH + C*LLHH -> A+B+C
 *)
CompactWeinbergCoefficient[terms_List] := Module[{pieces, combined},
  If[terms === {}, Return[Missing["NoHolomorphicTerms"]]];

  pieces = StripWeinbergOperatorStructure /@ terms;
  combined = Total[pieces];

  Quiet@Check[
    FactorTerms[Cancel[Together[combined]]],
    Simplify[combined]
  ]
];

(*
1. Check for existence of coefficient
2. Use ExtractWeinbergTerms to get our Weinberg terms
3. Sort our terms into holomorphic and conjugate terms
4. We use CompactWeinbergCoefficient to get our coefficient
5. We return a report on this
*)
ExtractWeinbergCoefficient[eft_] := Module[
  {
    present,
    terms,
    holomorphicTerms,
    conjugateTerms,
    sector,
    holomorphicSector,
    coefficient,
    status
  },

  (* Check if weinberg charged conjugated gamma is present in eft *)
  present = !FreeQ[eft, GammaCC];

  If[!TrueQ[present],
    Return[<|
      "Status" -> "NotFound",
      "Present" -> False,
      "TermCount" -> 0,
      "Terms" -> {},
      "Sector" -> 0,
      "Coefficient" -> 0
    |>]
  ];

  (* Gets all Weinberg terms *)
  terms = ExtractWeinbergTerms[eft];

  If[terms === {},
    Return[<|
      "Status" -> "PresentIsolationPending",
      "Present" -> True,
      "TermCount" -> 0,
      "Terms" -> {},
      "Sector" -> 0,
      "Coefficient" -> Missing["PendingCanonicalisation"]
    |>]
  ];

  sector = Simplify[Expand[Total[terms]]];

  (* C5 uses the P_L LLHH sector, we look at holomorphic terms;
   keep P_R only as the Hermitian-conjugate audit. *)
  holomorphicTerms = Select[terms, !FreeQ[#, Proj[-1]] &];
  conjugateTerms = Select[terms, !FreeQ[#, Proj[1]] &];

  If[holomorphicTerms === {},
    Return[<|
      "Status" -> "PresentHolomorphicIsolationPending",
      "Present" -> True,
      "TermCount" -> Length[terms],
      "HolomorphicTermCount" -> 0,
      "ConjugateTermCount" -> Length[conjugateTerms],
      "Terms" -> terms,
      "HolomorphicTerms" -> {},
      "ConjugateTerms" -> conjugateTerms,
      "Sector" -> sector,
      "HolomorphicSector" -> 0,
      "Coefficient" -> Missing["PendingCanonicalisation"]
    |>]
  ];

  holomorphicSector = Simplify[Expand[Total[holomorphicTerms]]];
  (* Obtain Weinberg Coefficient itself from our terms *)
  coefficient = CompactWeinbergCoefficient[holomorphicTerms];
  status = If[MissingQ[coefficient], "CoefficientPending", "Success"];

  <|
    "Status" -> status,
    "Present" -> True,
    "TermCount" -> Length[terms],
    "HolomorphicTermCount" -> Length[holomorphicTerms],
    "ConjugateTermCount" -> Length[conjugateTerms],
    "Terms" -> terms,
    "HolomorphicTerms" -> holomorphicTerms,
    "ConjugateTerms" -> conjugateTerms,
    "Sector" -> sector,
    "HolomorphicSector" -> holomorphicSector,
    "Coefficient" -> coefficient
  |>
];

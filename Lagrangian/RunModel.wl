(* RunModel.wl

   Run one T3 model from the command line, match it onto the EFT, extract C5,
   and write a compact machine-readable summary.

   CLI modes:
     CLASS: output-dir EFT-order loop-order T3-class alpha eg: (output 5 1 B -1) < If we use the classes
     DIMS : output-dir EFT-order loop-order DIMS dS1 dS2 dF alpha eg: (output 5 1 DIMS 3 5 4 0) < if we use non classified

    All this file does is
    1. Parser Translate inputs like m and p into -2 and +2 etc
    2. Parser Get the argument inputs of this file into a dictionary
    3. Get our helper functions form other WL files
    4. We check validity of model
    5. We use BuildT3Lagrangian
*)

ClearAll["Global`*"];
scriptDirectory = DirectoryName @ ExpandFileName[$InputFileName];


(* ------------------------------------------------------------------------- *)
(* Command-line input                                                         *)
(* ------------------------------------------------------------------------- *)

(* m2 -> -2, p2 -> +2 etc, translation device *)
ParseSignedCLI[s_String] := Which[
  StringMatchQ[s, "m" ~~ DigitCharacter ..], -ToExpression @ StringDrop[s, 1],
  StringMatchQ[s, "p" ~~ DigitCharacter ..],  ToExpression @ StringDrop[s, 1],
  True, ToExpression[s]
];


ParseThresholdPlanCLI[args_List] := Module[{token, payload, groups},
  token = SelectFirst[
    args,
    StringStartsQ[ToUpperCase[#], "THRESHOLDS="] &,
    Missing["NotFound"]
  ];

  If[MissingQ[token],
    Return[{{"F", "S1", "S2"}}]
  ];

  payload = StringDrop[token, StringLength["THRESHOLDS="]];
  groups = StringSplit[payload, ";"];

  Select[
    StringSplit[#, ","] & /@ groups,
    Length[#] > 0 &
  ]
];

ParseCLI[args_List] := Module[
  {output, eftOrder, loopOrder, mode, alpha, dS1, dS2, dF,
   debugReports, exportRGETensors, thresholdPlan},

  output = If[Length[args] >= 1, ExpandFileName @ args[[1]],
    FileNameJoin @ {scriptDirectory, "output", "T3"}];
  eftOrder = If[Length[args] >= 2, ToExpression @ args[[2]], 5];
  loopOrder = If[Length[args] >= 3, ToExpression @ args[[3]], 1];
  mode = If[Length[args] >= 4, ToUpperCase @ args[[4]], "B"];

  If[mode === "DIMS",
    If[Length[args] < 8, Return[$Failed]]; (* Needs enough dimensions for all fields *)
    {dS1, dS2, dF} = ToExpression /@ args[[5 ;; 7]]; 
    alpha = ParseSignedCLI @ args[[8]],
    alpha = If[Length[args] >= 5, ParseSignedCLI @ args[[5]], -1]
  ];

  debugReports = MemberQ[ToUpperCase /@ args, "DEBUG"];
  exportRGETensors = MemberQ[ToUpperCase /@ args, "RGETENSORS"];
  thresholdPlan = ParseThresholdPlanCLI[args];

  <|
    "Output" -> output, "EFTOrder" -> eftOrder, "LoopOrder" -> loopOrder,
    "Mode" -> mode, "Alpha" -> alpha,
    "Dimensions" -> If[mode === "DIMS", {dS1, dS2, dF}, None],
    "DebugReports" -> debugReports,
    "ExportRGETensors" -> exportRGETensors,
    "ThresholdPlan" -> thresholdPlan
  |>
];


(* ------------------------------------------------------------------------- *)
(* Output helpers                                                             *)
(* ------------------------------------------------------------------------- *)

(* Recieves Weinberg Coefficient information and Writes to files *)
ExportWeinbergArtifacts[data_Association, output_String, debugReports_] := Module[
  {
    terms,
    holomorphicTerms,
    conjugateTerms,
    coefficient,
    sectorTeX,
    coefficientTeX
  },

  terms = Lookup[data, "Terms", {}];
  holomorphicTerms = Lookup[data, "HolomorphicTerms", {}];
  conjugateTerms = Lookup[data, "ConjugateTerms", {}];
  coefficient = Lookup[data, "Coefficient", Missing["PendingCanonicalisation"]];
  sectorTeX = Lookup[data, "SectorTeX", <||>];
  coefficientTeX = Lookup[data, "CoefficientTeX", <||>];

  (* Always keep the isolated Weinberg terms as audit artifacts.  These are
     much easier to inspect than the full matched EFT and are needed to check
     flavor symmetrisation explicitly. *)
  If[Length[holomorphicTerms] > 0,
    Export[
      FileNameJoin @ {output, "weinberg_holomorphic_terms.txt"},
      StringRiffle[
        MapIndexed[
          "--- Holomorphic Weinberg term " <> ToString[First[#2]] <> " ---\n" <>
            ToString[#1, InputForm] &,
          holomorphicTerms
        ],
        "\n\n"
      ],
      "Text"
    ];
  ];

  If[Length[conjugateTerms] > 0,
    Export[
      FileNameJoin @ {output, "weinberg_conjugate_terms.txt"},
      StringRiffle[
        MapIndexed[
          "--- Hermitian-conjugate Weinberg term " <>
            ToString[First[#2]] <> " ---\n" <>
            ToString[#1, InputForm] &,
          conjugateTerms
        ],
        "\n\n"
      ],
      "Text"
    ];
  ];

  If[TrueQ[debugReports] && Length[terms] > 0,
    Export[
      FileNameJoin @ {output, "c5_raw.txt"},
      StringRiffle[
        MapIndexed[
          "--- Weinberg term " <> ToString[First[#2]] <> " ---\n" <>
            ToString[#1, InputForm] &,
          terms
        ],
        "\n\n"
      ],
      "Text"
    ];
    Export[FileNameJoin @ {output, "c5_raw.tex"}, Lookup[sectorTeX, "LaTeX", ""], "Text"];
  ];

  If[!MissingQ[coefficient],
    Export[
      FileNameJoin @ {output, "c5_coefficient.txt"},
      ToString[coefficient, InputForm],
      "Text"
    ];
    If[
      TrueQ[debugReports] &&
        TrueQ @ Lookup[coefficientTeX, "Success", False],
      Export[
        FileNameJoin @ {output, "c5_coefficient.tex"},
        Lookup[coefficientTeX, "LaTeX", ""],
        "Text"
      ]
    ];
  ];
];

(* Build our outputs *)
BuildSummary[model_, alpha_, build_, matching_, smMatching_, difference_, data_,
  free_, interaction_, uv_, bsm_, eft_, bsmEft_, debugReports_] := Module[
  {coefficient, canonicalCoefficient, coefficientTeX, canonicalCoefficientTeX, loopFunctionTeX, sectorTeX, terms},

  coefficient = Lookup[data, "Coefficient", Missing["PendingCanonicalisation"]];
  canonicalCoefficient = Lookup[
    data,
    "CanonicalCoefficient",
    Missing["PendingCanonicalisation"]
  ];
  coefficientTeX = Lookup[data, "CoefficientTeX", <||>];
  canonicalCoefficientTeX = Lookup[data, "CanonicalCoefficientTeX", <||>];
  loopFunctionTeX = Lookup[data, "LoopFunctionTeX", <||>];
  sectorTeX = Lookup[data, "SectorTeX", <||>];
  terms = Lookup[data, "Terms", {}];

  <|
    "ModelClass" -> model["Class"], "Alpha" -> alpha,
    "Scalar1SU2" -> model["Scalar1", "SU2"],
    "Scalar1Hypercharge" -> ToString[model["Scalar1", "Y"], InputForm],
    "Scalar2SU2" -> model["Scalar2", "SU2"],
    "Scalar2Hypercharge" -> ToString[model["Scalar2", "Y"], InputForm],
    "FermionSU2" -> model["Fermion", "SU2"],
    "FermionHypercharge" -> ToString[model["Fermion", "Y"], InputForm],

    "BuildStatus" -> build["Status"],
    "MatchingStatus" -> matching["Status"],
    "AcceptedInteractions" -> build["AllowedInteractions"],
    "RejectedInteractions" -> build["RejectedInteractions"],
    "T3IngredientsPresent" -> build["T3IngredientsPresent"],

    "WeinbergOperatorPresent" -> TrueQ @ Lookup[data, "Present", False],
    "WeinbergExtractionStatus" -> Lookup[data, "Status", "Unknown"],
    "WeinbergTermCount" -> Lookup[data, "TermCount", 0],
    "WeinbergHolomorphicTermCount" -> Lookup[data, "HolomorphicTermCount", 0],
    "WeinbergConjugateTermCount" -> Lookup[data, "ConjugateTermCount", 0],
    "WeinbergHolomorphicTermsFile" -> If[
      Length @ Lookup[data, "HolomorphicTerms", {}] > 0,
      "weinberg_holomorphic_terms.txt",
      ""
    ],
    "WeinbergConjugateTermsFile" -> If[
      Length @ Lookup[data, "ConjugateTerms", {}] > 0,
      "weinberg_conjugate_terms.txt",
      ""
    ],
    "WeinbergRawFile" -> If[
      TrueQ[debugReports] && Length[terms] > 0, "c5_raw.txt", ""
    ],
    "WeinbergRawLaTeXFile" -> If[
      TrueQ[debugReports] && Length[terms] > 0, "c5_raw.tex", ""
    ],
    "WeinbergCoefficientFile" -> If[!MissingQ[coefficient], "c5_coefficient.txt", ""],
    "WeinbergCoefficientLaTeXFile" -> If[
      TrueQ[debugReports] &&
        TrueQ @ Lookup[coefficientTeX, "Success", False],
      "c5_coefficient.tex",
      ""
    ],
    "WeinbergCoefficientInputForm" -> If[
      MissingQ[coefficient],
      "Pending exact Matchete operator canonicalisation",
      ToString[coefficient, InputForm]
    ],
    "WeinbergCoefficientLaTeX" -> Lookup[coefficientTeX, "LaTeX", ""],
    "WeinbergCanonicalCoefficientInputForm" -> If[
      MissingQ[canonicalCoefficient],
      "",
      ToString[canonicalCoefficient, InputForm]
    ],
    "WeinbergCanonicalCoefficientLaTeX" -> Lookup[
      canonicalCoefficientTeX,
      "LaTeX",
      ""
    ],
    "WeinbergCanonicalCoefficientConversionSuccess" -> TrueQ @ Lookup[
      canonicalCoefficientTeX,
      "Success",
      False
    ],
    "WeinbergLoopFunctionLaTeX" -> Lookup[loopFunctionTeX, "LaTeX", ""],
    "WeinbergLoopFunctionConversionSuccess" -> TrueQ @ Lookup[
      loopFunctionTeX,
      "Success",
      False
    ],
    "WeinbergOperatorConvention" -> Lookup[data, "OperatorConvention", ""],
    "WeinbergSectorLaTeX" -> Lookup[sectorTeX, "LaTeX", ""],
    "WeinbergSectorConversionSuccess" -> TrueQ @ Lookup[sectorTeX, "Success", False],
    "WeinbergCoefficientConversionSuccess" -> TrueQ @ Lookup[coefficientTeX, "Success", False],

    "UVConversionSuccess" -> TrueQ[uv["Success"]],
    "BSMUVConversionSuccess" -> TrueQ[bsm["Success"]],
    "EFTConversionSuccess" -> TrueQ[eft["Success"]],
    "BSMEFTConversionSuccess" -> TrueQ[bsmEft["Success"]],
    "SMMatchingStatus" -> Lookup[smMatching, "Status", "Unknown"],
    "SMBaselineMode" -> Lookup[smMatching, "BaselineMode", "Unknown"],
    "BSMDifferenceStatus" -> Lookup[difference, "Status", "Unknown"],
    "BSMDifferenceEOMFallbackUsed" -> TrueQ @ Lookup[
      difference, "EOMFallbackUsed", False
    ],
    "FreeLagrangianConversionSuccess" -> TrueQ[free["Success"]],
    "InteractionLagrangianConversionSuccess" -> TrueQ[interaction["Success"]],
    "FreeLagrangianLaTeX" -> free["LaTeX"],
    "InteractionLagrangianLaTeX" -> interaction["LaTeX"],
    "UVLagrangianLaTeX" -> uv["LaTeX"],
    "BSMUVLagrangianLaTeX" -> bsm["LaTeX"],
    "EFTLagrangianLaTeX" -> eft["LaTeX"],
    "BSMEFTLagrangianLaTeX" -> bsmEft["LaTeX"]
  |>
];


(* ------------------------------------------------------------------------- *)
(* Main pipeline                                                              *)
(* ------------------------------------------------------------------------- *)

(* We run parser *)
config = ParseCLI @ Rest[$ScriptCommandLine];
If[config === $Failed,
  Print["ERROR: DIMS mode requires dS1 dS2 dF alpha."];
  Exit[2]
];

outputDirectory = config["Output"];
If[!DirectoryQ[outputDirectory],
  CreateDirectory[outputDirectory, CreateIntermediateDirectories -> True]
];

Print["Loading Matchete..."];
matcheteLoaded = UsingFrontEnd[Needs["Matchete`"]; True];
If[!TrueQ[matcheteLoaded], Print["ERROR: Matchete failed to load."]; Exit[6]];
Print["Matchete loaded successfully."];

lagrangianDir = DirectoryName[$InputFileName];
projectRoot = DirectoryName[lagrangianDir];

(* Get modules and helpers from these files *)
physicsLaTeXFile =
    FileNameJoin[{lagrangianDir, "PhysicsLaTeX.wl"}];

modelCatalogFile =
    FileNameJoin[{lagrangianDir, "T3ModelCatalog.wl"}];

lagrangianBuilderFile =
    FileNameJoin[{lagrangianDir, "LagrangianBuilder.wl"}];

matchingFile =
    FileNameJoin[{lagrangianDir, "RunMatching.wl"}];

rgeTensorExporterFile = FileNameJoin[
  {projectRoot, "RGE", "general", "wolfram", "T3RGETensorExport.wl"}
];

Get[physicsLaTeXFile];
Get[modelCatalogFile];
Get[lagrangianBuilderFile];
Get[matchingFile];
If[TrueQ[config["ExportRGETensors"]], Get[rgeTensorExporterFile]];

(* Get model from SU2 dimensions *)
model = If[
  config["Mode"] === "DIMS",
  T3ModelFromDimensions[Sequence @@ config["Dimensions"], config["Alpha"]],
  T3ModelFromClass[config["Mode"], config["Alpha"]]
];
If[model === $Failed,
  If[config["Mode"] === "DIMS",
    Print["ERROR: requested SU(2) dimensions do not form a T3 topology."],
    Print["ERROR: unknown T3 class."]
  ];
  Exit[2]
];

(* Print models *)
Print["Model: ", model["Class"], ", alpha = ", config["Alpha"]];
Print["Scalar1 (d,Y) = ", {model["Scalar1", "SU2"], model["Scalar1", "Y"]}];
Print["Scalar2 (d,Y) = ", {model["Scalar2", "SU2"], model["Scalar2", "Y"]}];
Print["Fermion (d,Y) = ", {model["Fermion", "SU2"], model["Fermion", "Y"]}];

(* Configure the UV model so only the first requested threshold group is
   marked Heavy. The other T3 fields remain dynamical Light fields. *)
If[SetT3HeavyFields[First[config["ThresholdPlan"]]] === $Failed,
  Print["ERROR: invalid first threshold group."];
  Exit[9]
];
Print["Threshold plan: ", config["ThresholdPlan"]];
Print["Threshold group count: ", Length[config["ThresholdPlan"]]];

Export[
  FileNameJoin @ {outputDirectory, "threshold_plan_debug.json"},
  <|
    "ParsedThresholdPlan" -> config["ThresholdPlan"],
    "ParsedThresholdCount" -> Length[config["ThresholdPlan"]]
  |>,
  "RawJSON"
];

(* We build Lagrangian *)
build = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
If[!AssociationQ[build] || build === $Aborted || build === $Failed,
  Print["ERROR: build failed."]; Exit[10]
];
If[build["Status"] =!= "Success",
  Print["ERROR: full UV validation failed."]; Exit[11]
];
Print["Accepted interactions: ", build["AllowedInteractions"]];
Print["T3 ingredients present: ", build["T3IngredientsPresent"]];


(* All of this stuff is to help RGE tensors where we expose the indices of our terms 
May be removed
rgeTensorExportStatus = "NotRequested";
rgeTensorExchangeFile = "";
rgeTensorQuarticCount = 0;
rgeTensorYukawaCount = 0;

(* 
1. Check for RGETensors and store things into directory
2. Get results from ExportT3RGETensors and put in the path  *)
If[TrueQ[config["ExportRGETensors"]],
  dataDirectory = FileNameJoin[{outputDirectory, "data"}];
  If[!DirectoryQ[dataDirectory],
    CreateDirectory[dataDirectory, CreateIntermediateDirectories -> True]
  ];

  rgeTensorExchangePath = FileNameJoin[
    {dataDirectory, "t3_rge_tensor_exchange.json"}
  ];
  
  rgeTensorExportResult = Quiet @ Check[
    ExportT3RGETensors[model, rgeTensorExchangePath],
    $Failed
  ];

  If[rgeTensorExportResult === $Failed || !FileExistsQ[rgeTensorExchangePath],
    rgeTensorExportStatus = "Failed";
    Print["ERROR: T3 RGE tensor export failed."],
    rgeTensorExchange = Import[rgeTensorExchangePath, "RawJSON"];
    rgeTensorExportStatus = "Success";
    rgeTensorExchangeFile = "data/t3_rge_tensor_exchange.json";
    rgeTensorQuarticCount = Length @ Lookup[
      rgeTensorExchange, "quartic_components", {}
    ];
    rgeTensorYukawaCount = Length @ Lookup[
      rgeTensorExchange, "raw_yukawa_components", {}
    ];
    Print[
      "T3 RGE tensor export: Success (",
      rgeTensorQuarticCount, " quartic; ",
      rgeTensorYukawaCount, " Yukawa components)."
    ];
  ];
];
*)

(* We match our Lagrangian*)
matching = CheckAbort[
  RunSequentialT3Matching[
    build["LUV"],
    model,
    config["ThresholdPlan"],
    config["EFTOrder"],
    config["LoopOrder"],
    outputDirectory,
    ToString[model["Class"]]
  ],
  $Aborted
];
If[!AssociationQ[matching] || Lookup[matching, "Status", ""] =!= "Success",
  Print["ERROR: matching failed."]; Exit[12]
];

Print[
  "Sequential matching completed with ",
  Length @ Lookup[matching, "Stages", {}],
  " stage(s)."
];

Export[
  FileNameJoin @ {outputDirectory, "sequential_matching_debug.json"},
  <|
    "Status" -> Lookup[matching, "Status", "Unknown"],
    "RequestedThresholdCount" -> Length[config["ThresholdPlan"]],
    "MatchedStageCount" -> Length @ Lookup[matching, "Stages", {}],
    "StageLabels" -> Map[
      Function[stage,
        "EFT_" <> ToString[stage["Level"]] <> "_after_" <>
          StringRiffle[stage["IntegratedFields"], "_"]
      ],
      Lookup[matching, "Stages", {}]
    ]
  |>,
  "RawJSON"
];

(* Get the MatchedEFT, if that fails, get loopEFT *)
matchedEFT = Lookup[matching, "MatchedEFT", Lookup[matching, "LoopEFT", 0]];

(* We make SM Lagrangian. The builder returns LUV = LSM + LBSM, we get LSM by removing LBSM. *)
smLagrangian = Expand[build["LUV"] - build["LBSM"]];

(* We match just the SM to compare and to obtain just the BSM matched later *)
Print["\nStarting pure-SM baseline matching..."];
smMatching = CheckAbort[
  RunSMBaselineMatching[
    smLagrangian,
    config["EFTOrder"],
    config["LoopOrder"]
  ],
  $Aborted
];

If[
  !AssociationQ[smMatching] ||
    Lookup[smMatching, "Status", ""] =!= "Success",
  Print["ERROR: pure-SM baseline matching failed."];
  Exit[13]
];

(* get smEFT *)
smMatchedEFT = Lookup[
  smMatching,
  "MatchedEFT",
  Lookup[smMatching, "LoopEFT", smLagrangian]
];


(* Build compact, JSON-safe descriptions of every real sequential EFT stage.
   For comparison tables we subtract the pure-SM baseline directly from each
   already-canonical stage Lagrangian. *)
sequentialStageData = Map[
  Function[stage,
    Module[
      {level, integrated, active, stageLag, stageBSM, stageTeX, stageBSMTeX, label},

      level = stage["Level"];
      integrated = stage["IntegratedFields"];
      active = stage["ActiveHeavyFields"];
      stageLag = stage["Lagrangian"];
      stageBSM = Expand[stageLag - smMatchedEFT];

      label = "EFT_" <> ToString[level] <> "_after_" <>
        StringRiffle[integrated, "_"];

      stageTeX = ExpressionToLaTeX[stageLag];
      stageBSMTeX = ExpressionToLaTeX[stageBSM];

      <|
        "Level" -> level,
        "Label" -> label,
        "IntegratedFields" -> integrated,
        "ActiveHeavyFields" -> active,
        "EFTConversionSuccess" -> TrueQ[stageTeX["Success"]],
        "BSMEFTConversionSuccess" -> TrueQ[stageBSMTeX["Success"]],
        "EFTLagrangianLaTeX" -> stageTeX["LaTeX"],
        "BSMEFTLagrangianLaTeX" -> stageBSMTeX["LaTeX"]
      |>
    ]
  ],
  Lookup[matching, "Stages", {}]
];

(* get BSM EFT from subtraction *)
difference = CheckAbort[
  BuildMatchedEFTDifference[matchedEFT, smMatchedEFT],
  $Aborted
];

If[
  !AssociationQ[difference] ||
    Lookup[difference, "Status", ""] =!= "Success",
  Print["ERROR: BSM EFT difference construction failed."];
  Exit[14]
];

bsmMatchedEFT = Lookup[difference, "BSMEFT", 0];

(* Get latex of everything *)
freeTeX = ExpressionToLaTeX[build["LFree"]];
interactionTeX = ExpressionToLaTeX[build["LInt"]];
uvTeX = ExpressionToLaTeX[build["LUV"]];
bsmTeX = ExpressionToLaTeX[build["LBSM"]];
eftTeX = ExpressionToLaTeX[matchedEFT];
bsmEftTeX = ExpressionToLaTeX[bsmMatchedEFT];

(* Get Weinberg Coefficient from matchedEFT *)
weinbergData = ExtractWeinbergCoefficient[matchedEFT];
weinbergCoefficient = Lookup[
  weinbergData,
  "Coefficient",
  Missing["PendingCanonicalisation"]
];
weinbergCanonicalCoefficient = Lookup[
  weinbergData,
  "CanonicalCoefficient",
  Missing["PendingCanonicalisation"]
];

weinbergData = Join[
  weinbergData,
  <|
    "SectorTeX" -> ExpressionToLaTeX @ Lookup[weinbergData, "Sector", 0],
    "CoefficientTeX" -> If[
      MissingQ[weinbergCoefficient],
      <|"Success" -> False, "LaTeX" -> ""|>,
      ExpressionToLaTeX[weinbergCoefficient]
    ],
    "CanonicalCoefficientTeX" -> If[
      MissingQ[weinbergCanonicalCoefficient],
      <|"Success" -> False, "LaTeX" -> ""|>,
      ExpressionToLaTeX[weinbergCanonicalCoefficient]
    ],
    "LoopFunctionTeX" -> ExpressionToLaTeX[
      T3LoopI[
        Coupling[MF, {}, 0]^2,
        Coupling[MS1, {}, 0]^2,
        Coupling[MS2, {}, 0]^2
      ]
    ]
  |>
];


(* Write the Weinberg files *)
ExportWeinbergArtifacts[
  weinbergData,
  outputDirectory,
  config["DebugReports"]
];
summary = BuildSummary[
  model,
  config["Alpha"],
  build,
  matching,
  smMatching,
  difference,
  weinbergData,
  freeTeX,
  interactionTeX,
  uvTeX,
  bsmTeX,
  eftTeX,
  bsmEftTeX,
  config["DebugReports"]
];
(*
  summary = Join[
    summary,
    <|
      "T3RGETensorExportStatus" -> rgeTensorExportStatus,
      "T3RGETensorExchangeFile" -> rgeTensorExchangeFile,
      "T3RGEQuarticComponentCount" -> rgeTensorQuarticCount,
      "T3RGEYukawaComponentCount" -> rgeTensorYukawaCount
    |>
  ];
*)
If[
  Length[sequentialStageData] =!= Length[config["ThresholdPlan"]],
  Print[
    "ERROR: stage-report export count mismatch: ",
    Length[sequentialStageData],
    " vs ",
    Length[config["ThresholdPlan"]]
  ];
];

summary = Join[
  summary,
  <|
    "ThresholdPlan" -> config["ThresholdPlan"],
    "SequentialMatchingStatus" -> Lookup[matching, "Status", "Unknown"],
    "EFTStages" -> sequentialStageData
  |>
];

Export[FileNameJoin @ {outputDirectory, "comparison_summary.json"}, summary, "RawJSON"];

Print["Weinberg operator present: ", summary["WeinbergOperatorPresent"]];
If[summary["WeinbergOperatorPresent"],
  If[summary["WeinbergTermCount"] > 0,
    Print["Weinberg contributions isolated: ", summary["WeinbergTermCount"]],
    Print["Weinberg present; contribution isolation pending"]
  ];
  If[summary["WeinbergExtractionStatus"] === "Success",
    Print[
      "C5 extraction: Success (holomorphic terms: ",
      summary["WeinbergHolomorphicTermCount"],
      "; HC terms: ", summary["WeinbergConjugateTermCount"], ")"
    ],
    Print["C5 canonicalisation: pending"]
  ];
];

Print["Pure-SM matching status: ", summary["SMMatchingStatus"]];
Print["BSM EFT difference status: ", summary["BSMDifferenceStatus"]];
Print["BSM EFT conversion success: ", summary["BSMEFTConversionSuccess"]];

Exit[0];

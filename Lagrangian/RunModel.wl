(* RunModel.wl

   Run one T3 model from the command line, match it onto the EFT, extract C5,
   and write a compact machine-readable summary.

   CLI modes:
     CLASS: output-dir EFT-order loop-order T3-class alpha eg: (output 5 1 B -1) < If we use the classes
     DIMS : output-dir EFT-order loop-order DIMS dS1 dS2 dF alpha eg: (output 5 1 DIMS 3 5 4 0) < if we use non classified
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

ParseCLI[args_List] := Module[
  {output, eftOrder, loopOrder, mode, alpha, dS1, dS2, dF,
   debugReports, debugPosition},

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

  debugPosition = If[mode === "DIMS", 9, 6];
  debugReports = Length[args] >= debugPosition &&
    ToUpperCase @ args[[debugPosition]] === "DEBUG";

  <|
    "Output" -> output, "EFTOrder" -> eftOrder, "LoopOrder" -> loopOrder,
    "Mode" -> mode, "Alpha" -> alpha,
    "Dimensions" -> If[mode === "DIMS", {dS1, dS2, dF}, None],
    "DebugReports" -> debugReports
  |>
];


(* ------------------------------------------------------------------------- *)
(* Output helpers                                                             *)
(* ------------------------------------------------------------------------- *)

(* Recieves Weinberg Coefficient information and Writes to files *)
ExportWeinbergArtifacts[data_Association, output_String, debugReports_] := Module[
  {terms, coefficient, sectorTeX, coefficientTeX},

  terms = Lookup[data, "Terms", {}];
  coefficient = Lookup[data, "Coefficient", Missing["PendingCanonicalisation"]];
  sectorTeX = Lookup[data, "SectorTeX", <||>];
  coefficientTeX = Lookup[data, "CoefficientTeX", <||>];

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
  {coefficient, coefficientTeX, sectorTeX, terms},

  coefficient = Lookup[data, "Coefficient", Missing["PendingCanonicalisation"]];
  coefficientTeX = Lookup[data, "CoefficientTeX", <||>];
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

physicsLaTeXFile =
    FileNameJoin[{lagrangianDir, "PhysicsLaTeX.wl"}];

modelCatalogFile =
    FileNameJoin[{lagrangianDir, "T3ModelCatalog.wl"}];

lagrangianBuilderFile =
    FileNameJoin[{lagrangianDir, "LagrangianBuilder.wl"}];

matchingFile =
    FileNameJoin[{lagrangianDir, "RunMatching.wl"}];

Get[physicsLaTeXFile];
Get[modelCatalogFile];
Get[lagrangianBuilderFile];
Get[matchingFile];

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

(* We match *)
matching = CheckAbort[
  RunT3Matching[build["LUV"], config["EFTOrder"], config["LoopOrder"]],
  $Aborted
];
If[!AssociationQ[matching] || Lookup[matching, "Status", ""] =!= "Success",
  Print["ERROR: matching failed."]; Exit[12]
];

matchedEFT = Lookup[matching, "MatchedEFT", Lookup[matching, "LoopEFT", 0]];

(* The builder returns LUV = LSM + LBSM, we get LSM be removing LBSM. *)
smLagrangian = Expand[build["LUV"] - build["LBSM"]];

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

freeTeX = ExpressionToLaTeX[build["LFree"]];
interactionTeX = ExpressionToLaTeX[build["LInt"]];
uvTeX = ExpressionToLaTeX[build["LUV"]];
bsmTeX = ExpressionToLaTeX[build["LBSM"]];
eftTeX = ExpressionToLaTeX[matchedEFT];
bsmEftTeX = ExpressionToLaTeX[bsmMatchedEFT];

(* Get Weinberg Coefficient from matchedEFT *)
weinbergData = ExtractWeinbergCoefficient[matchedEFT];
weinbergCoefficient = Lookup[weinbergData, "Coefficient", Missing["PendingCanonicalisation"]];
weinbergData = Join[
  weinbergData,
  <|
    "SectorTeX" -> ExpressionToLaTeX @ Lookup[weinbergData, "Sector", 0],
    "CoefficientTeX" -> If[
      MissingQ[weinbergCoefficient],
      <|"Success" -> False, "LaTeX" -> ""|>,
      ExpressionToLaTeX[weinbergCoefficient]
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

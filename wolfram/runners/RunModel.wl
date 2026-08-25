(* RunModel.wl
   Run one T3 model from the command line, match it onto the EFT, extract C5,
   and write a compact machine-readable summary.

   CLI modes:
     CLASS: output-dir EFT-order loop-order T3-class alpha
     DIMS : output-dir EFT-order loop-order DIMS dS1 dS2 dF alpha
*)

ClearAll["Global`*"];
scriptDirectory = DirectoryName @ ExpandFileName[$InputFileName];


(* ------------------------------------------------------------------------- *)
(* Command-line input                                                         *)
(* ------------------------------------------------------------------------- *)

ParseSignedCLI[s_String] := Which[
  StringMatchQ[s, "m" ~~ DigitCharacter ..], -ToExpression @ StringDrop[s, 1],
  StringMatchQ[s, "p" ~~ DigitCharacter ..],  ToExpression @ StringDrop[s, 1],
  True, ToExpression[s]
];

ParseCLI[args_List] := Module[
  {output, eftOrder, loopOrder, mode, alpha, dS1, dS2, dF},

  output = If[Length[args] >= 1, ExpandFileName @ args[[1]],
    FileNameJoin @ {scriptDirectory, "output", "T3"}];
  eftOrder = If[Length[args] >= 2, ToExpression @ args[[2]], 5];
  loopOrder = If[Length[args] >= 3, ToExpression @ args[[3]], 1];
  mode = If[Length[args] >= 4, ToUpperCase @ args[[4]], "B"];

  If[mode === "DIMS",
    If[Length[args] < 8, Return[$Failed]];
    {dS1, dS2, dF} = ToExpression /@ args[[5 ;; 7]];
    alpha = ParseSignedCLI @ args[[8]],
    alpha = If[Length[args] >= 5, ParseSignedCLI @ args[[5]], -1]
  ];

  <|
    "Output" -> output, "EFTOrder" -> eftOrder, "LoopOrder" -> loopOrder,
    "Mode" -> mode, "Alpha" -> alpha,
    "Dimensions" -> If[mode === "DIMS", {dS1, dS2, dF}, None]
  |>
];


(* ------------------------------------------------------------------------- *)
(* Output helpers                                                             *)
(* ------------------------------------------------------------------------- *)

ExportWeinbergArtifacts[data_Association, output_String] := Module[
  {terms, coefficient, sectorTeX, coefficientTeX},

  terms = Lookup[data, "Terms", {}];
  coefficient = Lookup[data, "Coefficient", Missing["PendingCanonicalisation"]];
  sectorTeX = Lookup[data, "SectorTeX", <||>];
  coefficientTeX = Lookup[data, "CoefficientTeX", <||>];

  If[Length[terms] > 0,
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
    If[TrueQ @ Lookup[coefficientTeX, "Success", False],
      Export[
        FileNameJoin @ {output, "c5_coefficient.tex"},
        Lookup[coefficientTeX, "LaTeX", ""],
        "Text"
      ]
    ];
  ];
];

BuildSummary[model_, alpha_, build_, matching_, data_, uv_, bsm_, eft_] := Module[
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
    "WeinbergRawFile" -> If[Length[terms] > 0, "c5_raw.txt", ""],
    "WeinbergRawLaTeXFile" -> If[Length[terms] > 0, "c5_raw.tex", ""],
    "WeinbergCoefficientFile" -> If[!MissingQ[coefficient], "c5_coefficient.txt", ""],
    "WeinbergCoefficientLaTeXFile" -> If[
      TrueQ @ Lookup[coefficientTeX, "Success", False], "c5_coefficient.tex", ""
    ],
    "WeinbergCoefficientInputForm" -> If[
      MissingQ[coefficient],
      "Pending exact Matchete operator canonicalisation",
      ToString[coefficient, InputForm]
    ],
    "WeinbergCoefficientLaTeX" -> Lookup[coefficientTeX, "LaTeX", ""],
    "WeinbergSectorConversionSuccess" -> TrueQ @ Lookup[sectorTeX, "Success", False],
    "WeinbergCoefficientConversionSuccess" -> TrueQ @ Lookup[coefficientTeX, "Success", False],

    "UVConversionSuccess" -> TrueQ[uv["Success"]],
    "BSMUVConversionSuccess" -> TrueQ[bsm["Success"]],
    "EFTConversionSuccess" -> TrueQ[eft["Success"]],
    "UVLagrangianLaTeX" -> uv["LaTeX"],
    "BSMUVLagrangianLaTeX" -> bsm["LaTeX"],
    "EFTLagrangianLaTeX" -> eft["LaTeX"]
  |>
];


(* ------------------------------------------------------------------------- *)
(* Main pipeline                                                              *)
(* ------------------------------------------------------------------------- *)

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

Get[FileNameJoin @ {scriptDirectory, "..", "core", "PhysicsLaTeX.wl"}];
Get[FileNameJoin @ {scriptDirectory, "..", "t3", "T3ModelCatalog.wl"}];
Get[FileNameJoin @ {scriptDirectory, "..", "t3", "LagrangianBuilder.wl"}];
Get[FileNameJoin @ {scriptDirectory, "..", "matching", "RunMatching.wl"}];

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

Print["Model: ", model["Class"], ", alpha = ", config["Alpha"]];
Print["Scalar1 (d,Y) = ", {model["Scalar1", "SU2"], model["Scalar1", "Y"]}];
Print["Scalar2 (d,Y) = ", {model["Scalar2", "SU2"], model["Scalar2", "Y"]}];
Print["Fermion (d,Y) = ", {model["Fermion", "SU2"], model["Fermion", "Y"]}];

(* The generic builder uses Matchete group-theory routines that can require
   front-end services. Keep that dependency scoped to model construction. *)
build = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
If[!AssociationQ[build] || build === $Aborted || build === $Failed,
  Print["ERROR: build failed."]; Exit[10]
];
If[build["Status"] =!= "Success",
  Print["ERROR: full UV validation failed."]; Exit[11]
];
Print["Accepted interactions: ", build["AllowedInteractions"]];
Print["T3 ingredients present: ", build["T3IngredientsPresent"]];

matching = CheckAbort[
  RunT3Matching[build["LUV"], config["EFTOrder"], config["LoopOrder"]],
  $Aborted
];
If[!AssociationQ[matching] || Lookup[matching, "Status", ""] =!= "Success",
  Print["ERROR: matching failed."]; Exit[12]
];

matchedEFT = Lookup[matching, "MatchedEFT", Lookup[matching, "LoopEFT", 0]];
uvTeX = ExpressionToLaTeX[build["LUV"]];
bsmTeX = ExpressionToLaTeX[build["LBSM"]];
eftTeX = ExpressionToLaTeX[matchedEFT];

(* Weinberg detection and coefficient extraction are kept together so the
   summary and exported diagnostics are derived from the same result. *)
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

ExportWeinbergArtifacts[weinbergData, outputDirectory];
summary = BuildSummary[
  model, config["Alpha"], build, matching, weinbergData, uvTeX, bsmTeX, eftTeX
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

Exit[0];

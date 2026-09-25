(* RunT3ScalarOnlyRGBeta.wl
   Renormalisable running between thresholds after the heavy fermion F is removed.

   The underlying RGBeta model-builder function names are historical and are
   retained here to avoid changing the verified Wolfram algebra in this
   source-cleanup batch.
*)

ClearAll["Global`*"];
Needs["RGBeta`"];

scriptDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{scriptDir, "T3RGBetaRunnerCommon.wl"}]];
Get[FileNameJoin[{scriptDir, "T3RGBetaModel.wl"}]];

config = T3RGBetaParseRunnerArguments[];

sharedMode = config["SharedMode"];
dS1 = config["dS1"];
dS2 = config["dS2"];
dF = config["dF"];
alpha = config["alpha"];
outputPath = config["OutputPath"];

build = CheckAbort[
    Quiet[
        If[
            sharedMode,
            T3RGBetaBuildSharedEFT1[dS1, dF],
            T3RGBetaBuildEFT1[dS1, dS2, dF, alpha]
        ]
    ],
    $Aborted
];

If[!AssociationQ[build],
    T3RGBetaWriteJSON[
        outputPath,
        <|
            "status" -> "Failed",
            "stage" -> "Build",
            "dimensions" -> {dS1, dS2, dF},
            "alpha" -> alpha
        |>
    ];
    Exit[1];
];

jsonMetadata = Join[
    T3RGBetaCommonMetadata[config, build],
    <|
        "IntegratedField" -> "F",
        "ActiveBSMFields" -> If[sharedMode, {"S"}, {"S1", "S2"}]
    |>
];

betaAssociation = CheckAbort[
    Quiet[
        If[
            sharedMode,
            T3RGBetaSharedEFT1OneLoopBetas[],
            T3RGBetaEFT1OneLoopBetas[]
        ]
    ],
    $Aborted
];

If[betaAssociation === $Aborted || !AssociationQ[betaAssociation],
    T3RGBetaWriteJSON[
        outputPath,
        <|
            "status" -> "Failed",
            "stage" -> "BetaGeneration",
            "metadata" -> jsonMetadata
        |>
    ];
    Exit[1];
];

stringBetas = T3RGBetaSerializeAssociation[betaAssociation];
reportBetaAssociation = T3RGBetaConventionalReportBetas[betaAssociation];
stringReportBetas = T3RGBetaSerializeAssociation[reportBetaAssociation];
latexReportBetas = T3RGBetaLaTeXAssociation[reportBetaAssociation];

result = <|
    "status" -> "Success",
    "metadata" -> Append[
        jsonMetadata,
        "ReportBetaConvention" ->
            "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"
    ],
    "betas" -> stringBetas,
    "report_betas" -> stringReportBetas,
    "report_beta_latex" -> latexReportBetas
|>;

T3RGBetaWriteJSON[outputPath, result];

Print["RGBeta T3 scalar-only renormalisable RGE export: SUCCESS"];
Print["Output: ", outputPath];

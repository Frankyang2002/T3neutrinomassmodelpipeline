(* 
Command line runner for full UV RGBeta
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

(* Build our UV model *)
build = CheckAbort[
    Quiet[
        If[
            sharedMode,
            T3RGBetaBuildShared[dS1, dF],
            T3RGBetaBuild[dS1, dS2, dF, alpha]
        ]
    ],
    $Aborted
];

(* If building fails we have this *)
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


(* Builds the jsonmeta data, as some things dont work with JSON *)
jsonMetadata = Join[
    T3RGBetaCommonMetadata[config, build],
    <|
        "YF" -> ToString[InputForm[build["YF"]]],
        "MajoranaFermion" -> TrueQ[build["MajoranaFermion"]]
    |>
];

(* We get our one loop beta function calculation here *)
betaAssociation = CheckAbort[
    Quiet[
        If[
            sharedMode,
            T3RGBetaSharedOneLoopBetas[],
            T3RGBetaOneLoopBetas[]
        ]
    ],
    $Aborted
];

(* Fail condition for rgbeta *)
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

(* Reverse RGBeta convention to our convention *)
reportBetaAssociation = T3RGBetaConventionalReportBetas[betaAssociation];

stringReportBetas = T3RGBetaSerializeAssociation[reportBetaAssociation];

(* TeXForm is retained as a first-pass rendering.  RGBeta has internal heads
   such as Matrix, Trans and Bar that TeXForm does not know how to typeset.
   Reports/RGEComparison.py performs the final physics-aware display cleanup
   while using these exact expressions. *)
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

Print["RGBeta T3 UV RGE export: SUCCESS"];
Print["Output: ", outputPath];

(* RunT3RGBeta.wl
   Command-line runner for the validated d<=3 T3 RGBeta model.

   Usage:
       wolframscript -file RunT3RGBeta.wl dS1 dS2 dF alpha output.json

   Exact symbolic quantities are exported as InputForm strings because JSON
   has no representation for Mathematica Rational or symbolic expressions.
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
            T3RGBetaBuildShared[dS1, dF],
            T3RGBetaBuild[dS1, dS2, dF, alpha]
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


(* JSON cannot encode exact Mathematica Rational objects.  Preserve the
   project's exact hypercharge convention as InputForm strings. *)
jsonMetadata = Join[
    T3RGBetaCommonMetadata[config, build],
    <|
        "YF" -> ToString[InputForm[build["YF"]]],
        "MajoranaFermion" -> TrueQ[build["MajoranaFermion"]]
    |>
];


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

(* RGBeta uses a special convention for gauge couplings:

       beta_g^RGBeta = d(g^2)/d ln(mu).

   For the human-facing report we instead use the conventional derivative

       16 pi^2 d g/d ln(mu).

   At one loop this is BetaTerm[g,1]/(2 g).  All non-gauge couplings already
   use the ordinary derivative convention, so their one-loop report term is
   simply BetaTerm[X,1].  The original BetaTerm output above is preserved
   unchanged in "betas" for machine use. *)
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

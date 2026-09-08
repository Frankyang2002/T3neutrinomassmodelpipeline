(* RunT3RGBeta.wl
   Command-line runner for the validated d<=3 T3 RGBeta model.

   Usage:
       wolframscript -file RunT3RGBeta.wl dS1 dS2 dF alpha output.json

   Exact symbolic quantities are exported as InputForm strings because JSON
   has no representation for Mathematica Rational or symbolic expressions.
*)

ClearAll["Global`*"];
Needs["RGBeta`"];

If[Length[$ScriptCommandLine] < 6,
    Print[
        "Usage: wolframscript -file RunT3RGBeta.wl ",
        "dS1 dS2 dF alpha output.json (negative integers use mN, e.g. m1)"
    ];
    Exit[2];
];

scriptDir = DirectoryName[$InputFileName];

Get[
    FileNameJoin[{
        scriptDir,
        "T3RGBetaModel.wl"
    }]
];

ParseIntegerToken[token_String] := If[
    StringStartsQ[token, "m"],
    -ToExpression[StringDrop[token, 1]],
    ToExpression[token]
];

dS1 = ParseIntegerToken[$ScriptCommandLine[[2]]];
dS2 = ParseIntegerToken[$ScriptCommandLine[[3]]];
dF = ParseIntegerToken[$ScriptCommandLine[[4]]];
alpha = ParseIntegerToken[$ScriptCommandLine[[5]]];
outputPath = $ScriptCommandLine[[6]];


WriteJSON[payload_Association] := Module[{json},
    json = ExportString[payload, "RawJSON"];
    If[!StringQ[json],
        Print["JSON serialization failed."];
        Exit[1];
    ];
    Export[outputPath, json, "Text"];
];


build = CheckAbort[
    Quiet[T3RGBetaBuild[dS1, dS2, dF, alpha]],
    $Aborted
];

If[!AssociationQ[build],
    WriteJSON[
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
jsonMetadata = <|
    "dS1" -> dS1,
    "dS2" -> dS2,
    "dF" -> dF,
    "alpha" -> alpha,
    "YS1" -> ToString[InputForm[build["YS1"]]],
    "YS2" -> ToString[InputForm[build["YS2"]]],
    "YF" -> ToString[InputForm[build["YF"]]],
    "MajoranaFermion" -> TrueQ[build["MajoranaFermion"]]
|>;


betaAssociation = CheckAbort[
    Quiet[T3RGBetaOneLoopBetas[]],
    $Aborted
];

If[betaAssociation === $Aborted || !AssociationQ[betaAssociation],
    WriteJSON[
        <|
            "status" -> "Failed",
            "stage" -> "BetaGeneration",
            "metadata" -> jsonMetadata
        |>
    ];
    Exit[1];
];


stringBetas = Association @ KeyValueMap[
    (#1 -> ToString[InputForm[#2]]) &,
    betaAssociation
];

(* RGBeta uses a special convention for gauge couplings:

       beta_g^RGBeta = d(g^2)/d ln(mu).

   For the human-facing report we instead use the conventional derivative

       16 pi^2 d g/d ln(mu).

   At one loop this is BetaTerm[g,1]/(2 g).  All non-gauge couplings already
   use the ordinary derivative convention, so their one-loop report term is
   simply BetaTerm[X,1].  The original BetaTerm output above is preserved
   unchanged in "betas" for machine use. *)
reportBetaAssociation = Association @ KeyValueMap[
    Function[{name, beta},
        name -> Switch[
            name,
            "gY", Cancel[beta/(2 gY)],
            "g2", Cancel[beta/(2 g2)],
            "g3", Cancel[beta/(2 g3)],
            _, beta
        ]
    ],
    betaAssociation
];

stringReportBetas = Association @ KeyValueMap[
    (#1 -> ToString[InputForm[#2]]) &,
    reportBetaAssociation
];

(* TeXForm is retained as a first-pass rendering.  RGBeta has internal heads
   such as Matrix, Trans and Bar that TeXForm does not know how to typeset.
   Reports/RGEComparison.py performs the final physics-aware display cleanup
   while using these exact expressions. *)
latexReportBetas = Association @ KeyValueMap[
    (#1 -> ToString[TeXForm[#2]]) &,
    reportBetaAssociation
];


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

WriteJSON[result];

Print["RGBeta T3 UV RGE export: SUCCESS"];
Print["Output: ", outputPath];

(* T3RGBetaRunnerCommon.wl
   Shared command-line, JSON, and report-beta utilities for the UV and EFT1
   RGBeta runners.

   Physics/model construction remains in T3RGBetaModel.wl and is selected
   explicitly by each stage-specific runner.
*)

T3RGBetaParseIntegerToken[token_String] := If[
    StringStartsQ[token, "m"],
    -ToExpression[StringDrop[token, 1]],
    ToExpression[token]
];


T3RGBetaParseRunnerArguments[] := Module[
    {sharedMode, dS, dS1, dS2, dF, alpha, outputPath},

    sharedMode = (
        Length[$ScriptCommandLine] >= 2 &&
        ToUpperCase[$ScriptCommandLine[[2]]] === "SHARED"
    );

    If[sharedMode,
        If[Length[$ScriptCommandLine] < 5,
            Print["Usage: ... SHARED dS dF output.json"];
            Exit[2];
        ],
        If[Length[$ScriptCommandLine] < 6,
            Print["Usage: ... dS1 dS2 dF alpha output.json"];
            Exit[2];
        ]
    ];

    If[sharedMode,
        dS = T3RGBetaParseIntegerToken[$ScriptCommandLine[[3]]];
        dF = T3RGBetaParseIntegerToken[$ScriptCommandLine[[4]]];
        dS1 = dS;
        dS2 = dS;
        alpha = -1;
        outputPath = $ScriptCommandLine[[5]],
        dS1 = T3RGBetaParseIntegerToken[$ScriptCommandLine[[2]]];
        dS2 = T3RGBetaParseIntegerToken[$ScriptCommandLine[[3]]];
        dF = T3RGBetaParseIntegerToken[$ScriptCommandLine[[4]]];
        alpha = T3RGBetaParseIntegerToken[$ScriptCommandLine[[5]]];
        outputPath = $ScriptCommandLine[[6]]
    ];

    <|
        "SharedMode" -> sharedMode,
        "dS1" -> dS1,
        "dS2" -> dS2,
        "dF" -> dF,
        "alpha" -> alpha,
        "OutputPath" -> outputPath
    |>
];


T3RGBetaWriteJSON[outputPath_String, payload_Association] := Module[{json},
    json = ExportString[payload, "RawJSON"];
    If[!StringQ[json],
        Print["JSON serialization failed."];
        Exit[1];
    ];
    Export[outputPath, json, "Text"];
];


T3RGBetaCommonMetadata[
    config_Association,
    build_Association
] := Module[
    {
        sharedMode = config["SharedMode"],
        dS1 = config["dS1"],
        dS2 = config["dS2"],
        dF = config["dF"],
        alpha = config["alpha"]
    },

    <|
        "SharedScalar" -> sharedMode,
        "dS" -> If[sharedMode, dS1, Null],
        "dS1" -> dS1,
        "dS2" -> dS2,
        "dF" -> dF,
        "alpha" -> alpha,
        "YS1" -> ToString[
            InputForm[If[sharedMode, -build["YS"], build["YS1"]]]
        ],
        "YS2" -> ToString[
            InputForm[If[sharedMode, build["YS"], build["YS2"]]]
        ],
        "YS" -> ToString[
            InputForm[If[sharedMode, build["YS"], Null]]
        ]
    |>
];


T3RGBetaSerializeAssociation[association_Association] :=
    Association @ KeyValueMap[
        (#1 -> ToString[InputForm[#2]]) &,
        association
    ];


T3RGBetaConventionalReportBetas[betaAssociation_Association] :=
    Association @ KeyValueMap[
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


T3RGBetaLaTeXAssociation[association_Association] :=
    Association @ KeyValueMap[
        (#1 -> ToString[TeXForm[#2]]) &,
        association
    ];

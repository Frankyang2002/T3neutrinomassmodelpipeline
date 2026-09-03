(* TestT3RGBetaModelBuilderDiagnostic_v2.wl
   Diagnose beta generation one coupling at a time.

   Important:
   RGBeta emits harmless messages in headless wolframscript.
   We therefore do NOT use Check[...] around successful model construction,
   because Check converts any message into $Failed.
*)

ClearAll["Global`*"];
Needs["RGBeta`"];

Get[
    FileNameJoin[{
        DirectoryName[$InputFileName],
        "..",
        "..",
        "RGE",
        "running",
        "wolfram",
        "T3RGBetaModel.wl"
    }]
];

models = {
    {"A", {1, 3, 2}, 0},
    {"B", {2, 2, 1}, -1},
    {"C", {2, 2, 3}, -1},
    {"D", {3, 1, 2}, -2},
    {"E", {3, 3, 2}, 0}
};

couplings = {
    {"gY", gY},
    {"g2", g2},
    {"g3", g3},
    {"yu", yu},
    {"yd", yd},
    {"ye", ye},
    {"y1", y1},
    {"y2", y2},
    {"MF", MF},
    {"mS1Sq", mS1Sq},
    {"mS2Sq", mS2Sq}
};

ClearAll[TryBeta];

TryBeta[name_String, coupling_] := Module[{result},
    Print["    ", name, " ..."];

    result = CheckAbort[
        Quiet[BetaTerm[coupling, 1]],
        $Aborted
    ];

    Which[
        result === $Aborted,
            Print["      ABORTED"],
        Head[result] === BetaTerm,
            Print["      UNEVALUATED"],
        True,
            Print["      SUCCESS"]
    ];

    result
];


Print[""];
Print["================ RGBeta T3 MODEL BUILDER DIAGNOSTIC V2 ================"];

Do[
    label = model[[1]];
    dims = model[[2]];
    alpha = model[[3]];

    (* Do not use Check here: ordinary RGBeta messages are not build failures. *)
    build = CheckAbort[
        Quiet[
            T3RGBetaBuild[
                dims[[1]],
                dims[[2]],
                dims[[3]],
                alpha
            ]
        ],
        $Aborted
    ];

    Print[""];
    Print["T3-", label, " build: ",
        Which[
            build === $Aborted, "ABORTED",
            build === $Failed, "FAILED",
            AssociationQ[build], "SUCCESS",
            True, "UNKNOWN"
        ]
    ];

    If[AssociationQ[build],
        Print["  metadata: ", build];
        Print["  beta functions:"];

        Do[
            result = TryBeta[item[[1]], item[[2]]];

            If[result === $Aborted, Break[]],
            {item, couplings}
        ];
    ],
    {model, models}
];

Print[""];
Print["RGBETA T3 MODEL BUILDER DIAGNOSTIC V2: COMPLETE"];

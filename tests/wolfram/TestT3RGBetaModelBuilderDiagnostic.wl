(* TestT3RGBetaModelBuilderDiagnostic.wl
   Diagnose beta generation one coupling at a time so an Abort does not hide
   which RGBeta sector caused the problem.
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
    result = Quiet @ CheckAbort[
        Check[BetaTerm[coupling, 1], $Failed],
        $Aborted
    ];

    Which[
        result === $Aborted,
            Print["      ABORTED"],
        result === $Failed,
            Print["      FAILED"],
        True,
            Print["      SUCCESS"]
    ];

    result
];


Print[""];
Print["================ RGBeta T3 MODEL BUILDER DIAGNOSTIC ================"];

Do[
    label = model[[1]];
    dims = model[[2]];
    alpha = model[[3]];

    build = Quiet @ CheckAbort[
        Check[
            T3RGBetaBuild[
                dims[[1]],
                dims[[2]],
                dims[[3]],
                alpha
            ],
            $Failed
        ],
        $Aborted
    ];

    Print[""];
    Print["T3-", label, " build: ",
        Which[
            build === $Aborted, "ABORTED",
            build === $Failed, "FAILED",
            True, "SUCCESS"
        ]
    ];

    If[AssociationQ[build],
        Print["  metadata: ", build];
        Print["  beta functions:"];

        Do[
            result = TryBeta[item[[1]], item[[2]]];

            (* Stop this model at the first Abort.  ResetModel in the next
               model iteration gives us a clean state. *)
            If[result === $Aborted, Break[]],
            {item, couplings}
        ];
    ],
    {model, models}
];

Print[""];
Print["RGBETA T3 MODEL BUILDER DIAGNOSTIC: COMPLETE"];

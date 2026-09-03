(* TestT3RGBetaFullModel.wl
   Final A-E smoke test for the complete current d<=3 UV RGBeta builder.
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
    {"gY", gY}, {"g2", g2}, {"g3", g3},
    {"yu", yu}, {"yd", yd}, {"ye", ye},
    {"y1", y1}, {"y2", y2},
    {"MF", MF}, {"mS1Sq", mS1Sq}, {"mS2Sq", mS2Sq},
    {"lambdaH", lambdaH},
    {"lambdaS1", lambdaS1}, {"lambdaS2", lambdaS2},
    {"lambdaH1", lambdaH1}, {"lambdaH2", lambdaH2},
    {"lambda12", lambda12}, {"lambdaT3", lambdaT3}
};

ClearAll[TryBeta];

TryBeta[name_String, coupling_] := Module[{result, residual},
    result = CheckAbort[
        Quiet[BetaTerm[coupling, 1]],
        $Aborted
    ];

    If[result === $Aborted,
        Print["    ", name, ": ABORTED"];
        Return[$Aborted];
    ];

    residual = DeleteDuplicates @ Cases[
        result,
        _del | _delS2 | _fStruct | _tGen | _eps,
        Infinity
    ];

    Print[
        "    ", name, ": ",
        If[residual === {}, "SUCCESS", "SUCCESS (residual tensors)"]
    ];

    If[residual =!= {},
        Print["      residual: ", residual]
    ];

    result
];


Print[""];
Print["================ FULL RGBeta T3 A-E MODEL ================"];

Do[
    label = model[[1]];
    dims = model[[2]];
    alpha = model[[3]];

    build = CheckAbort[
        Quiet[
            T3RGBetaBuild[
                dims[[1]], dims[[2]], dims[[3]], alpha
            ]
        ],
        $Aborted
    ];

    Print[""];
    Print[
        "T3-", label, " build: ",
        If[AssociationQ[build], "SUCCESS", "FAILED"]
    ];

    If[AssociationQ[build],
        Print["  metadata: ", build];
        Print["  beta functions:"];

        Do[
            TryBeta[item[[1]], item[[2]]],
            {item, couplings}
        ];
    ],
    {model, models}
];

Print[""];
Print["FULL RGBETA T3 A-E MODEL CHECK: COMPLETE"];

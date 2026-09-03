(* TestT3RGBetaRGClosedFull_v3.wl *)

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
    {"A", {1,3,2}, 0},
    {"B", {2,2,1}, -1},
    {"C", {2,2,3}, -1},
    {"D", {3,1,2}, -2},
    {"E", {3,3,2}, 0}
};

baseNames = {
    "gY", "g2", "g3",
    "yu", "yd", "ye", "y1", "y2",
    "MF", "mS1Sq", "mS2Sq",
    "lambdaH", "lambdaS1", "lambdaS2",
    "lambdaH1", "lambdaH2", "lambda12", "lambdaT3"
};

expectedExtras = <|
    "A" -> {
        "lambdaH2Adj",
        "lambdaS2Adj"
    },
    "B" -> {
        "lambdaH1Adj",
        "lambdaH2Adj",
        "lambda12Adj",
        "lambdaHHdagS2S2",
        "lambdaHHdagS1barS1bar",
        "lambdaS1bar2S2bar2",
        "lambdaS1barS2S2bar2",
        "lambdaS1S1bar2S2bar",
        "lambdaHHdagS1barS2barCross"
    },
    "C" -> {
        "lambdaH1Adj",
        "lambdaH2Adj",
        "lambda12Adj",
        "lambdaHHdagS2S2",
        "lambdaHHdagS1barS1bar",
        "lambdaS1bar2S2bar2",
        "lambdaS1barS2S2bar2",
        "lambdaS1S1bar2S2bar",
        "lambdaHHdagS1barS2barCross"
    },
    "D" -> {
        "lambdaH1Adj",
        "lambdaS1Adj"
    },
    "E" -> {
        "lambdaH1Adj",
        "lambdaH2Adj",
        "lambdaS1Adj",
        "lambdaS2Adj",
        "lambda12Adj",
        "lambda12Cross"
    }
|>;

HasInternalTensor[expr_] := Length @ Cases[
    expr,
    _RGBeta`PackageScope`QuarticTensors |
    _RGBeta`PackageScope`UpsilonQuarticTensors,
    Infinity
] > 0;

HasResidualGroupTensor[expr_] := Length @ Cases[
    expr,
    _del | _delS2 | _fStruct | _tGen | _eps,
    Infinity
] > 0;

Print[""];
Print["================ STRICT RG-CLOSED T3 A-E TEST v3 ================"];

Do[
    label = model[[1]];
    dims = model[[2]];
    alpha = model[[3]];

    Print[""];
    Print["T3-", label, ":"];

    build = CheckAbort[
        Quiet[T3RGBetaBuild[dims[[1]], dims[[2]], dims[[3]], alpha]],
        $Aborted
    ];

    If[!AssociationQ[build],
        Print["  build: FAILED"];
        Continue[];
    ];

    betas = CheckAbort[
        Quiet[T3RGBetaOneLoopBetas[]],
        $Aborted
    ];

    If[betas === $Aborted || !AssociationQ[betas],
        Print["  beta association: FAILED"];
        Continue[];
    ];

    expected = Join[baseNames, expectedExtras[label]];
    missing = Complement[expected, Keys[betas]];
    unexpected = Complement[Keys[betas], expected];

    Print["  beta count: ", Length[betas], " / expected ", Length[expected]];
    Print["  missing beta names: ", missing];
    Print["  unexpected beta names: ", unexpected];

    allGood = (missing === {} && unexpected === {});

    KeyValueMap[
        Function[{name, expr},
            status = Which[
                HasInternalTensor[expr],
                    "INCONCLUSIVE (internal RGBeta tensors)",
                HasResidualGroupTensor[expr],
                    "INCONCLUSIVE (residual group tensors)",
                True,
                    "SUCCESS"
            ];

            Print["    ", name, ": ", status];

            If[status =!= "SUCCESS",
                allGood = False;
                Print["      ", InputForm[expr]];
            ];
        ],
        betas
    ];

    Print[
        "  model result: ",
        If[allGood, "PASS", "CHECK REQUIRED"]
    ];

    ,
    {model, models}
];

Print[""];
Print["STRICT RG-CLOSED T3 A-E TEST v3: COMPLETE"];

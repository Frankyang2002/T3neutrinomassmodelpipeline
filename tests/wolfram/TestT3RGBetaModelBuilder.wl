(* TestT3RGBetaModelBuilder.wl *)

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

Print[""];
Print["================ RGBeta T3 MODEL BUILDER ================"];

Do[
    label = model[[1]];
    dims = model[[2]];
    alpha = model[[3]];

    result = Check[
        Quiet @ T3RGBetaBuild[
            dims[[1]],
            dims[[2]],
            dims[[3]],
            alpha
        ],
        $Failed
    ];

    Print[""];
    Print["T3-", label, " build: ",
        If[result === $Failed, "FAILED", "SUCCESS"]
    ];

    If[result =!= $Failed,
        Print["  metadata: ", result];

        betas = Check[
            Quiet @ T3RGBetaOneLoopBetas[],
            $Failed
        ];

        Print[
            "  one-loop beta generation: ",
            If[betas === $Failed, "FAILED", "SUCCESS"]
        ];

        If[betas =!= $Failed,
            Print[
                "  failed beta entries: ",
                Select[
                    Keys[betas],
                    betas[#] === $Failed &
                ]
            ];
        ];
    ],
    {model, models}
];

Print[""];
Print["RGBETA T3 MODEL BUILDER CHECK: COMPLETE"];

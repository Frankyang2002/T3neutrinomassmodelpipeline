(* TestRGBetaT3BasicQuartics.wl
   Test the norm-type scalar potential for the current d <= 3 T3 models.

   This deliberately excludes lambdaT3 at first.  We want to separate:
     1. ordinary self/portal quartics;
     2. the topology-defining HH S1 S2^\dagger invariant.

   Potential convention matches the current Matchete builder:

     V ⊃ lambdaH/2 (H†H)^2
        + lambdaS1/2 (S1†S1)^2
        + lambdaS2/2 (S2†S2)^2
        + lambdaH1 (H†H)(S1†S1)
        + lambdaH2 (H†H)(S2†S2)
        + lambda12 (S1†S1)(S2†S2).

   RGBeta's quartic field ordering is explicit, so each invariant below merely
   contracts each conjugate/non-conjugate pair in its own SU(2) representation.
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

ClearAll[PairInvariant, AddBasicQuartics, TryQuarticBeta];

PairInvariant[1, i_, j_] := 1;
PairInvariant[2, i_, j_] := del[SU2L @ fund, i, j];
PairInvariant[3, i_, j_] := del[SU2L @ S2, i, j];


AddBasicQuartics[dS1_Integer, dS2_Integer] := Module[{},
    SetReal[lambdaH, lambdaS1, lambdaS2, lambdaH1, lambdaH2, lambda12];

    AddQuartic[
        lambdaH,
        {Bar @ H, H, Bar @ H, H},
        GroupInvariant -> (
            del[SU2L @ fund, #1, #2] *
            del[SU2L @ fund, #3, #4] / 2 &
        )
    ];

    AddQuartic[
        lambdaS1,
        {Bar @ S1, S1, Bar @ S1, S1},
        GroupInvariant -> (
            PairInvariant[dS1, #1, #2] *
            PairInvariant[dS1, #3, #4] / 2 &
        )
    ];

    AddQuartic[
        lambdaS2,
        {Bar @ S2, S2, Bar @ S2, S2},
        GroupInvariant -> (
            PairInvariant[dS2, #1, #2] *
            PairInvariant[dS2, #3, #4] / 2 &
        )
    ];

    AddQuartic[
        lambdaH1,
        {Bar @ H, H, Bar @ S1, S1},
        GroupInvariant -> (
            del[SU2L @ fund, #1, #2] *
            PairInvariant[dS1, #3, #4] &
        )
    ];

    AddQuartic[
        lambdaH2,
        {Bar @ H, H, Bar @ S2, S2},
        GroupInvariant -> (
            del[SU2L @ fund, #1, #2] *
            PairInvariant[dS2, #3, #4] &
        )
    ];

    AddQuartic[
        lambda12,
        {Bar @ S1, S1, Bar @ S2, S2},
        GroupInvariant -> (
            PairInvariant[dS1, #1, #2] *
            PairInvariant[dS2, #3, #4] &
        )
    ];
];


TryQuarticBeta[name_String, coupling_] := Module[{result},
    Print["    ", name, " ..."];
    result = CheckAbort[
        Quiet[BetaTerm[coupling, 1]],
        $Aborted
    ];

    Print[
        "      ",
        Which[
            result === $Aborted, "ABORTED",
            Head[result] === BetaTerm, "UNEVALUATED",
            True, "SUCCESS"
        ]
    ];

    result
];


models = {
    {"A", {1, 3, 2}, 0},
    {"B", {2, 2, 1}, -1},
    {"C", {2, 2, 3}, -1},
    {"D", {3, 1, 2}, -2},
    {"E", {3, 3, 2}, 0}
};

quartics = {
    {"lambdaH", lambdaH},
    {"lambdaS1", lambdaS1},
    {"lambdaS2", lambdaS2},
    {"lambdaH1", lambdaH1},
    {"lambdaH2", lambdaH2},
    {"lambda12", lambda12}
};


Print[""];
Print["================ RGBeta T3 BASIC QUARTICS ================"];

Do[
    label = model[[1]];
    dims = model[[2]];
    alpha = model[[3]];

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
    Print[
        "T3-", label, " base build: ",
        If[AssociationQ[build], "SUCCESS", "FAILED"]
    ];

    If[AssociationQ[build],
        quarticBuild = CheckAbort[
            Quiet[
                AddBasicQuartics[dims[[1]], dims[[2]]];
                "Success"
            ],
            $Aborted
        ];

        Print["  quartic registration: ", quarticBuild];

        If[quarticBuild === "Success",
            Do[
                TryQuarticBeta[item[[1]], item[[2]]],
                {item, quartics}
            ];
        ];
    ],
    {model, models}
];

Print[""];
Print["RGBETA T3 BASIC QUARTICS CHECK: COMPLETE"];

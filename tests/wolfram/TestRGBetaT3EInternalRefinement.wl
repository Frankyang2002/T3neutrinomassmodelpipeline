(* TestRGBetaT3EInternalRefinement.wl
   T3-E refinement using RGBeta's actual package-scope function.
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

Print[""];
Print["================ RGBeta T3-E INTERNAL REFINEMENT ================"];

Print["Matching symbols:"];
Print[Names["RGBeta`*RefineGroupStructures*"]];
Print[Names["RGBeta`PackageScope`*RefineGroupStructures*"]];

build = Quiet @ T3RGBetaBuild[3, 3, 2, 0];

If[!AssociationQ[build],
    Print["base build: FAILED"];
    Exit[1]
];

Print["base build: SUCCESS"];

AddQuartic[
    lambdaT3,
    {H, H, S1, Bar @ S2},
    GroupInvariant -> (
        delS2[SU2L, a, #1, #2] *
        fStruct[SU2L, a, #3, #4] &
    ),
    SelfConjugate -> False
];

Print["lambdaT3 registration: SUCCESS"];

raw = Quiet @ BetaTerm[lambdaT3, 1];

Print[""];
Print["--- Raw beta ---"];
Print[InputForm[raw]];

refined = CheckAbort[
    Quiet[
        RGBeta`PackageScope`RefineGroupStructures[raw]
    ],
    $Aborted
];

Print[""];
Print[
    "internal refinement: ",
    If[refined === $Aborted, "ABORTED", "RETURNED"]
];

If[refined =!= $Aborted,
    refined = Simplify[refined];

    Print[""];
    Print["--- Refined beta ---"];
    Print[InputForm[refined]];

    residual = DeleteDuplicates @ Cases[
        refined,
        _del | _delS2 | _fStruct | _tGen,
        Infinity
    ];

    Print[""];
    Print["residual group tensors: ", residual];

    Print[
        "fully scalar: ",
        If[residual === {}, "YES", "NO"]
    ];
];

Print[""];
Print["RGBETA T3-E INTERNAL REFINEMENT CHECK: COMPLETE"];

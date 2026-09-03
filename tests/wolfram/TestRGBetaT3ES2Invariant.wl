(* TestRGBetaT3ES2Invariant.wl
   T3-E lambdaT3 using an invariant written purely in the S2 basis.

   The SU(2) triplet is represented by the symmetric product of two
   fundamentals.  Instead of fStruct[SU2L,...], which carries adjoint indices,
   construct the antisymmetric 3-triplet invariant entirely from delS2 and eps.

   This keeps S1, S2 and all triplet fermions/scalars in the same RGBeta
   representation label SU2L[S2].
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

ClearAll[T3ES2TripleInvariant];

(* Three-triplet invariant in the S2 basis.

   a,b,c are S2 indices.
   Each delS2 maps one triplet index to a symmetric pair of fundamentals.
   The three epsilon tensors close the six fundamental indices in a triangle.
*)
T3ES2TripleInvariant[a_, b_, c_] := Module[
    {i, j, k, l, m, n},
    delS2[SU2L, a, i, j] *
    delS2[SU2L, b, k, l] *
    delS2[SU2L, c, m, n] *
    eps[SU2L @ fund, i, k] *
    eps[SU2L @ fund, j, m] *
    eps[SU2L @ fund, l, n]
];


Print[""];
Print["================ RGBeta T3-E S2-BASIS INVARIANT ================"];

build = Quiet @ T3RGBetaBuild[3, 3, 2, 0];

Print[
    "base build: ",
    If[AssociationQ[build], "SUCCESS", "FAILED"]
];

If[AssociationQ[build],

    reg = CheckAbort[
        Quiet[
            AddQuartic[
                lambdaT3,
                {H, H, S1, Bar @ S2},
                GroupInvariant -> (
                    delS2[SU2L, a, #1, #2] *
                    T3ES2TripleInvariant[a, #3, #4] &
                ),
                SelfConjugate -> False
            ];
            "Success"
        ],
        $Aborted
    ];

    Print["lambdaT3 registration: ", reg];

    If[reg === "Success",

        beta = CheckAbort[
            Quiet[BetaTerm[lambdaT3, 1]],
            $Aborted
        ];

        Print[
            "beta(lambdaT3): ",
            If[beta === $Aborted, "ABORTED", "SUCCESS"]
        ];

        If[beta =!= $Aborted,
            beta = Simplify[beta];

            Print[""];
            Print["--- beta(lambdaT3) ---"];
            Print[InputForm[beta]];

            residual = DeleteDuplicates @ Cases[
                beta,
                _del | _delS2 | _fStruct | _tGen | _eps,
                Infinity
            ];

            Print[""];
            Print["residual group tensors: ", residual];
            Print[
                "fully scalar: ",
                If[residual === {}, "YES", "NO"]
            ];
        ];
    ];
];

Print[""];
Print["RGBETA T3-E S2-BASIS INVARIANT CHECK: COMPLETE"];

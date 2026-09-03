(* TestRGBetaT3ERefinement.wl
   Refine the T3-E lambdaT3 beta function.

   T3-E:
       S1 ~ 3_0
       S2 ~ 3_1
       F  ~ 2_{1/2}

   The raw beta function is generated successfully but contains residual
   adjoint/fundamental tensor structures.  This probe applies RGBeta's group
   structure refinement/simplification routines to determine whether those
   tensors reduce to ordinary scalar coefficients.
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

ClearAll[T3EMixInvariant, AddT3EMixing];

T3EMixInvariant :=
    (
        delS2[SU2L, a, #1, #2] *
        fStruct[SU2L, a, #3, #4] &
    );

AddT3EMixing[] := AddQuartic[
    lambdaT3,
    {H, H, S1, Bar @ S2},
    GroupInvariant -> T3EMixInvariant,
    SelfConjugate -> False
];


Print[""];
Print["================ RGBeta T3-E REFINEMENT ================"];

build = CheckAbort[
    Quiet[T3RGBetaBuild[3, 3, 2, 0]],
    $Aborted
];

Print[
    "base build: ",
    If[AssociationQ[build], "SUCCESS", "FAILED"]
];

If[AssociationQ[build],

    reg = CheckAbort[
        Quiet[
            AddT3EMixing[];
            "Success"
        ],
        $Aborted
    ];

    Print["lambdaT3 registration: ", reg];

    If[reg === "Success",

        raw = CheckAbort[
            Quiet[BetaTerm[lambdaT3, 1]],
            $Aborted
        ];

        Print[
            "raw beta: ",
            If[raw === $Aborted, "ABORTED", "SUCCESS"]
        ];

        If[raw =!= $Aborted,

            Print[""];
            Print["--- Raw expression ---"];
            Print[InputForm[raw]];

            refined1 = CheckAbort[
                Quiet[RefineGroupStructures[raw]],
                $Aborted
            ];

            Print[""];
            Print[
                "RefineGroupStructures: ",
                If[refined1 === $Aborted, "ABORTED", "SUCCESS"]
            ];

            If[refined1 =!= $Aborted,
                Print["--- Refined expression ---"];
                Print[InputForm[refined1]];
            ];

            refined2 = If[
                refined1 === $Aborted,
                $Aborted,
                CheckAbort[
                    Quiet[
                        Simplify[
                            RefineGroupStructures[refined1]
                        ]
                    ],
                    $Aborted
                ]
            ];

            Print[""];
            Print[
                "second refinement + Simplify: ",
                If[refined2 === $Aborted, "ABORTED", "SUCCESS"]
            ];

            If[refined2 =!= $Aborted,
                Print["--- Final expression ---"];
                Print[InputForm[refined2]];

                residual = Cases[
                    refined2,
                    _del | _delS2 | _fStruct | _tGen,
                    Infinity
                ];

                Print[""];
                Print[
                    "residual group tensors: ",
                    DeleteDuplicates[residual]
                ];
            ];
        ];
    ];
];

Print[""];
Print["RGBETA T3-E REFINEMENT CHECK: COMPLETE"];

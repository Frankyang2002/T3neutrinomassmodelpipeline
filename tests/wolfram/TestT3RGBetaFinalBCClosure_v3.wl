(* TestT3RGBetaFinalBCClosure_v3.wl
   Final B/C-only RG-closure check.

   Starting point:
     production T3RGBetaModel.wl with the 27-beta B/C scalar basis.

   Remaining omitted gauge+Z2 quartics from the original completeness scan:
     1. H Hdag S1dag S2dag + h.c. [contraction 1]

   The other previously omitted B/C structures were added in v3.

   We re-check the remaining protected contraction in the presence of the
   complete enlarged B/C scalar basis.
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

ClearAll[
    BuildFresh,
    Probe,
    HasInternalTensor,
    HasResidualGroupTensor
];

BuildFresh[dims_List, alpha_Integer] :=
    Quiet[T3RGBetaBuild[dims[[1]], dims[[2]], dims[[3]], alpha]];

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

Probe[
    label_String,
    dims_List,
    alpha_Integer,
    coupling_Symbol,
    addFunction_
] := Module[
    {build, reg, beta, zeroBeta, classification},

    Print[""];
    Print["    ", label];

    build = CheckAbort[BuildFresh[dims, alpha], $Aborted];

    If[!AssociationQ[build],
        Print["      build FAILED"];
        Return["INCONCLUSIVE"];
    ];

    reg = CheckAbort[
        Quiet[
            addFunction[];
            "Success"
        ],
        $Aborted
    ];

    If[reg =!= "Success",
        Print["      registration ABORTED"];
        Return["INCONCLUSIVE"];
    ];

    beta = CheckAbort[
        Quiet[BetaTerm[coupling, 1]],
        $Aborted
    ];

    If[beta === $Aborted,
        Print["      beta ABORTED"];
        Return["INCONCLUSIVE"];
    ];

    zeroBeta = Simplify[
        beta /. {
            coupling -> 0,
            Bar[coupling] -> 0
        }
    ];

    classification = Which[
        HasInternalTensor[zeroBeta],
            "INCONCLUSIVE",
        HasResidualGroupTensor[zeroBeta],
            "INCONCLUSIVE",
        TrueQ[zeroBeta === 0],
            "PROTECTED",
        True,
            "GENERATED"
    ];

    Print["      beta at candidate coupling = 0:"];
    Print["      ", InputForm[zeroBeta]];
    Print["      classification: ", classification];

    classification
];


Print[""];
Print["================ FINAL T3-B/C RG-CLOSURE TEST v3 ================"];

Do[
    dims = If[label === "B", {2,2,1}, {2,2,3}];
    alpha = -1;
    coupling = If[label === "B", kBFinal, kCFinal];

    Print[""];
    Print["T3-", label];

    Probe[
        "H Hdag S1dag S2dag + h.c. [remaining contraction]",
        dims,
        alpha,
        coupling,
        With[{cc = coupling},
            Function[{},
                AddQuartic[
                    cc,
                    {H, Bar @ H, Bar @ S1, Bar @ S2},
                    GroupInvariant -> (
                        del[SU2L @ fund, #1, #2] *
                        eps[SU2L @ fund, #3, #4] &
                    ),
                    SelfConjugate -> False
                ]
            ]
        ]
    ];

    ,
    {label, {"B", "C"}}
];

Print[""];
Print["FINAL T3-B/C RG-CLOSURE TEST v3: COMPLETE"];

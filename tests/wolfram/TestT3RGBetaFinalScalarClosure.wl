(* TestT3RGBetaFinalScalarClosure.wl
   Final scalar-potential RG-closure test.

   Starting point:
     the EXPANDED production T3RGBetaModel.wl, including all class-dependent
     quartics already found to be radiatively generated.

   This test probes every field structure that remains omitted from the
   gauge+Z2 quartic scan for the five d<=3 T3 benchmark classes.

   For structures with multiplicity > 1, linearly independent contractions
   are tested separately.

   Classification:
     GENERATED    beta(kappa)|_{kappa=0} != 0
     PROTECTED    beta(kappa)|_{kappa=0} == 0
     INCONCLUSIVE unresolved RGBeta internals/tensors remain
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
Print["================ FINAL T3 SCALAR RG-CLOSURE TEST ================"];


(* ================================================================== *)
(* T3-A : (dS1,dS2,dF) = (1,3,2), alpha=0                            *)
(* ================================================================== *)

dims = {1,3,2};
alpha = 0;

Print[""];
Print["T3-A"];

Probe[
    "S1dag^2 S2 S2dag + h.c.",
    dims, alpha, kA1,
    Function[{},
        AddQuartic[
            kA1,
            {Bar @ S1, Bar @ S1, S2, Bar @ S2},
            GroupInvariant -> (
                del[SU2L @ S2, #3, #4] &
            ),
            SelfConjugate -> False
        ]
    ]
];

Probe[
    "S1dag^4 + h.c.",
    dims, alpha, kA2,
    Function[{},
        AddQuartic[
            kA2,
            {Bar @ S1, Bar @ S1, Bar @ S1, Bar @ S1},
            GroupInvariant -> (1 &),
            SelfConjugate -> False
        ]
    ]
];

Probe[
    "S1 S1dag^3 + h.c.",
    dims, alpha, kA3,
    Function[{},
        AddQuartic[
            kA3,
            {S1, Bar @ S1, Bar @ S1, Bar @ S1},
            GroupInvariant -> (1 &),
            SelfConjugate -> False
        ]
    ]
];

Probe[
    "Hdag^2 S1 S2 + h.c.",
    dims, alpha, kA4,
    Function[{},
        AddQuartic[
            kA4,
            {Bar @ H, Bar @ H, S1, S2},
            GroupInvariant -> (
                delS2[SU2L, #4, #1, #2] &
            ),
            SelfConjugate -> False
        ]
    ]
];

Probe[
    "H Hdag S1dag^2 + h.c.",
    dims, alpha, kA5,
    Function[{},
        AddQuartic[
            kA5,
            {H, Bar @ H, Bar @ S1, Bar @ S1},
            GroupInvariant -> (
                del[SU2L @ fund, #1, #2] &
            ),
            SelfConjugate -> False
        ]
    ]
];


(* ================================================================== *)
(* T3-B/C : remaining omitted structures                              *)
(* ================================================================== *)

Do[
    dims = If[label === "B", {2,2,1}, {2,2,3}];
    alpha = -1;

    Print[""];
    Print["T3-", label];

    c1 = If[label === "B", kB1, kC1];
    c2 = If[label === "B", kB2, kC2];
    c3 = If[label === "B", kB3, kC3];
    c4 = If[label === "B", kB4, kC4];

    Probe[
        "S1dag S2 S2dag^2 + h.c.",
        dims, alpha, c1,
        With[{cc = c1},
            Function[{},
                AddQuartic[
                    cc,
                    {Bar @ S1, S2, Bar @ S2, Bar @ S2},
                    GroupInvariant -> (
                        delS2[SU2L, a, #1, #2] *
                        delS2[SU2L, a, #3, #4] &
                    ),
                    SelfConjugate -> False
                ]
            ]
        ]
    ];

    Probe[
        "S1 S1dag^2 S2dag + h.c.",
        dims, alpha, c2,
        With[{cc = c2},
            Function[{},
                AddQuartic[
                    cc,
                    {S1, Bar @ S1, Bar @ S1, Bar @ S2},
                    GroupInvariant -> (
                        delS2[SU2L, a, #1, #4] *
                        delS2[SU2L, a, #2, #3] &
                    ),
                    SelfConjugate -> False
                ]
            ]
        ]
    ];

    (* Two independent contractions for H Hdag S1dag S2dag. *)
    Probe[
        "H Hdag S1dag S2dag + h.c. [contraction 1]",
        dims, alpha, c3,
        With[{cc = c3},
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

    Probe[
        "H Hdag S1dag S2dag + h.c. [contraction 2]",
        dims, alpha, c4,
        With[{cc = c4},
            Function[{},
                AddQuartic[
                    cc,
                    {H, Bar @ H, Bar @ S1, Bar @ S2},
                    GroupInvariant -> (
                        eps[SU2L @ fund, #1, #3] *
                        eps[SU2L @ fund, #2, #4] &
                    ),
                    SelfConjugate -> False
                ]
            ]
        ]
    ];

    ,
    {label, {"B", "C"}}
];


(* ================================================================== *)
(* T3-D : (3,1,2), alpha=-2                                          *)
(* ================================================================== *)

dims = {3,1,2};
alpha = -2;

Print[""];
Print["T3-D"];

Probe[
    "S2dag^4 + h.c.",
    dims, alpha, kD1,
    Function[{},
        AddQuartic[
            kD1,
            {Bar @ S2, Bar @ S2, Bar @ S2, Bar @ S2},
            GroupInvariant -> (1 &),
            SelfConjugate -> False
        ]
    ]
];

Probe[
    "S2 S2dag^3 + h.c.",
    dims, alpha, kD2,
    Function[{},
        AddQuartic[
            kD2,
            {S2, Bar @ S2, Bar @ S2, Bar @ S2},
            GroupInvariant -> (1 &),
            SelfConjugate -> False
        ]
    ]
];

Probe[
    "S1 S1dag S2dag^2 + h.c.",
    dims, alpha, kD3,
    Function[{},
        AddQuartic[
            kD3,
            {S1, Bar @ S1, Bar @ S2, Bar @ S2},
            GroupInvariant -> (
                del[SU2L @ S2, #1, #2] &
            ),
            SelfConjugate -> False
        ]
    ]
];

Probe[
    "Hdag^2 S1dag S2dag + h.c.",
    dims, alpha, kD4,
    Function[{},
        AddQuartic[
            kD4,
            {Bar @ H, Bar @ H, Bar @ S1, Bar @ S2},
            GroupInvariant -> (
                delS2[SU2L, #3, #1, #2] &
            ),
            SelfConjugate -> False
        ]
    ]
];

Probe[
    "H Hdag S2dag^2 + h.c.",
    dims, alpha, kD5,
    Function[{},
        AddQuartic[
            kD5,
            {H, Bar @ H, Bar @ S2, Bar @ S2},
            GroupInvariant -> (
                del[SU2L @ fund, #1, #2] &
            ),
            SelfConjugate -> False
        ]
    ]
];


(* ================================================================== *)
(* T3-E : (3,3,2), alpha=0                                           *)
(* ================================================================== *)

dims = {3,3,2};
alpha = 0;

Print[""];
Print["T3-E"];

(* Two independent contractions for S1dag^2 S2 S2dag. *)
Probe[
    "S1dag^2 S2 S2dag + h.c. [singlet-like contraction]",
    dims, alpha, kE1,
    Function[{},
        AddQuartic[
            kE1,
            {Bar @ S1, Bar @ S1, S2, Bar @ S2},
            GroupInvariant -> (
                del[SU2L @ S2, #1, #2] *
                del[SU2L @ S2, #3, #4] &
            ),
            SelfConjugate -> False
        ]
    ]
];

Probe[
    "S1dag^2 S2 S2dag + h.c. [crossed symmetric contraction]",
    dims, alpha, kE2,
    Function[{},
        AddQuartic[
            kE2,
            {Bar @ S1, Bar @ S1, S2, Bar @ S2},
            GroupInvariant -> (
                del[SU2L @ S2, #1, #3] *
                del[SU2L @ S2, #2, #4] +
                del[SU2L @ S2, #1, #4] *
                del[SU2L @ S2, #2, #3] &
            ),
            SelfConjugate -> False
        ]
    ]
];

Probe[
    "S1dag^4 + h.c.",
    dims, alpha, kE3,
    Function[{},
        AddQuartic[
            kE3,
            {Bar @ S1, Bar @ S1, Bar @ S1, Bar @ S1},
            GroupInvariant -> (
                del[SU2L @ S2, #1, #2] *
                del[SU2L @ S2, #3, #4] &
            ),
            SelfConjugate -> False
        ]
    ]
];

Probe[
    "S1 S1dag^3 + h.c.",
    dims, alpha, kE4,
    Function[{},
        AddQuartic[
            kE4,
            {S1, Bar @ S1, Bar @ S1, Bar @ S1},
            GroupInvariant -> (
                del[SU2L @ S2, #1, #2] *
                del[SU2L @ S2, #3, #4] &
            ),
            SelfConjugate -> False
        ]
    ]
];

Probe[
    "Hdag^2 S1 S2 + h.c.",
    dims, alpha, kE5,
    Function[{},
        AddQuartic[
            kE5,
            {Bar @ H, Bar @ H, S1, S2},
            GroupInvariant -> (
                delS2[SU2L, a, #1, #2] *
                T3RGBetaS2TripleInvariant[a, #3, #4] &
            ),
            SelfConjugate -> False
        ]
    ]
];

Probe[
    "H Hdag S1dag^2 + h.c.",
    dims, alpha, kE6,
    Function[{},
        AddQuartic[
            kE6,
            {H, Bar @ H, Bar @ S1, Bar @ S1},
            GroupInvariant -> (
                del[SU2L @ fund, #1, #2] *
                del[SU2L @ S2, #3, #4] &
            ),
            SelfConjugate -> False
        ]
    ]
];


Print[""];
Print["FINAL T3 SCALAR RG-CLOSURE TEST: COMPLETE"];

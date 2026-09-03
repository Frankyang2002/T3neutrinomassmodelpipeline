(* TestT3RGBetaNewFieldStructureClosure_v2.wl
   Corrected new-field-structure RG-closure probe.

   IMPORTANT:
   Wolfram symbol names below contain only letters and digits.
   An underscore in Wolfram Language denotes a pattern and must NOT be used
   as part of a coupling name.

   A result is classified as:
       GENERATED     - beta at coupling=0 evaluated and is non-zero
       PROTECTED     - beta at coupling=0 evaluated and is exactly zero
       INCONCLUSIVE  - unresolved RGBeta internal functions remain
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
    HHTripletDoubletInvariant,
    TryClosure
];

HHTripletDoubletInvariant :=
    (
        delS2[SU2L, a, #1, #2] *
        delS2[SU2L, a, #3, #4] &
    );

BuildFresh[dims_List, alpha_Integer] :=
    Quiet[T3RGBetaBuild[dims[[1]], dims[[2]], dims[[3]], alpha]];


TryClosure[
    label_String,
    dims_List,
    alpha_Integer,
    coupling_Symbol,
    addFunction_
] := Module[
    {build, reg, beta, zeroBeta, internal, residual, classification},

    Print[""];
    Print["    ", label];

    build = CheckAbort[BuildFresh[dims, alpha], $Aborted];

    If[!AssociationQ[build],
        Print["      build FAILED"];
        Return[];
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
        Return[];
    ];

    beta = CheckAbort[
        Quiet[BetaTerm[coupling, 1]],
        $Aborted
    ];

    If[beta === $Aborted,
        Print["      beta ABORTED"];
        Return[];
    ];

    zeroBeta = Simplify[
        beta /. {
            coupling -> 0,
            Bar[coupling] -> 0
        }
    ];

    internal = Cases[
        zeroBeta,
        _RGBeta`PackageScope`QuarticTensors |
        _RGBeta`PackageScope`UpsilonQuarticTensors,
        Infinity
    ];

    residual = DeleteDuplicates @ Cases[
        zeroBeta,
        _del | _delS2 | _fStruct | _tGen | _eps,
        Infinity
    ];

    classification = Which[
        internal =!= {}, "INCONCLUSIVE",
        TrueQ[zeroBeta === 0], "PROTECTED",
        True, "GENERATED"
    ];

    Print["      beta at missing coupling = 0:"];
    Print["      ", InputForm[zeroBeta]];
    Print["      classification: ", classification];
    Print["      residual group tensors: ", residual];
];


models = <|
    "A" -> <|"dims" -> {1,3,2}, "alpha" -> 0|>,
    "B" -> <|"dims" -> {2,2,1}, "alpha" -> -1|>,
    "C" -> <|"dims" -> {2,2,3}, "alpha" -> -1|>,
    "D" -> <|"dims" -> {3,1,2}, "alpha" -> -2|>,
    "E" -> <|"dims" -> {3,3,2}, "alpha" -> 0|>
|>;


Print[""];
Print["================ T3 NEW FIELD-STRUCTURE CLOSURE v2 ================"];


(* ================================================================== *)
(* T3-A                                                               *)
(* ================================================================== *)

dims = models["A"]["dims"];
alpha = models["A"]["alpha"];

Print[""];
Print["T3-A"];

TryClosure[
    "S1^4 + h.c.",
    dims, alpha,
    kAS14,
    Function[{},
        AddQuartic[
            kAS14,
            {S1, S1, S1, S1},
            GroupInvariant -> (1 &),
            SelfConjugate -> False
        ]
    ]
];

TryClosure[
    "H Hdag S1dag^2 + h.c.",
    dims, alpha,
    kAHHS12,
    Function[{},
        AddQuartic[
            kAHHS12,
            {H, Bar @ H, Bar @ S1, Bar @ S1},
            GroupInvariant -> (
                del[SU2L @ fund, #1, #2] &
            ),
            SelfConjugate -> False
        ]
    ]
];

TryClosure[
    "Hdag^2 S1 S2 + h.c.",
    dims, alpha,
    kAH2S1S2,
    Function[{},
        AddQuartic[
            kAH2S1S2,
            {Bar @ H, Bar @ H, S1, S2},
            GroupInvariant -> (
                delS2[SU2L, #4, #1, #2] &
            ),
            SelfConjugate -> False
        ]
    ]
];


(* ================================================================== *)
(* T3-B                                                               *)
(* ================================================================== *)

dims = models["B"]["dims"];
alpha = models["B"]["alpha"];

Print[""];
Print["T3-B"];

TryClosure[
    "Hdag^2 S2^2 + h.c.",
    dims, alpha,
    kBH2S22,
    Function[{},
        AddQuartic[
            kBH2S22,
            {Bar @ H, Bar @ H, S2, S2},
            GroupInvariant -> HHTripletDoubletInvariant,
            SelfConjugate -> False
        ]
    ]
];

TryClosure[
    "Hdag^2 S1dag^2 + h.c.",
    dims, alpha,
    kBH2S112,
    Function[{},
        AddQuartic[
            kBH2S112,
            {Bar @ H, Bar @ H, Bar @ S1, Bar @ S1},
            GroupInvariant -> HHTripletDoubletInvariant,
            SelfConjugate -> False
        ]
    ]
];

TryClosure[
    "S1dag^2 S2dag^2 + h.c.",
    dims, alpha,
    kBS12S22,
    Function[{},
        AddQuartic[
            kBS12S22,
            {Bar @ S1, Bar @ S1, Bar @ S2, Bar @ S2},
            GroupInvariant -> (
                delS2[SU2L, a, #1, #2] *
                delS2[SU2L, a, #3, #4] &
            ),
            SelfConjugate -> False
        ]
    ]
];


(* ================================================================== *)
(* T3-C                                                               *)
(* ================================================================== *)

dims = models["C"]["dims"];
alpha = models["C"]["alpha"];

Print[""];
Print["T3-C"];

TryClosure[
    "Hdag^2 S2^2 + h.c.",
    dims, alpha,
    kCH2S22,
    Function[{},
        AddQuartic[
            kCH2S22,
            {Bar @ H, Bar @ H, S2, S2},
            GroupInvariant -> HHTripletDoubletInvariant,
            SelfConjugate -> False
        ]
    ]
];

TryClosure[
    "Hdag^2 S1dag^2 + h.c.",
    dims, alpha,
    kCH2S112,
    Function[{},
        AddQuartic[
            kCH2S112,
            {Bar @ H, Bar @ H, Bar @ S1, Bar @ S1},
            GroupInvariant -> HHTripletDoubletInvariant,
            SelfConjugate -> False
        ]
    ]
];

TryClosure[
    "S1dag^2 S2dag^2 + h.c.",
    dims, alpha,
    kCS12S22,
    Function[{},
        AddQuartic[
            kCS12S22,
            {Bar @ S1, Bar @ S1, Bar @ S2, Bar @ S2},
            GroupInvariant -> (
                delS2[SU2L, a, #1, #2] *
                delS2[SU2L, a, #3, #4] &
            ),
            SelfConjugate -> False
        ]
    ]
];


(* ================================================================== *)
(* T3-D                                                               *)
(* ================================================================== *)

dims = models["D"]["dims"];
alpha = models["D"]["alpha"];

Print[""];
Print["T3-D"];

TryClosure[
    "S2^4 + h.c.",
    dims, alpha,
    kDS24,
    Function[{},
        AddQuartic[
            kDS24,
            {S2, S2, S2, S2},
            GroupInvariant -> (1 &),
            SelfConjugate -> False
        ]
    ]
];

TryClosure[
    "H Hdag S2^2 + h.c.",
    dims, alpha,
    kDHHS22,
    Function[{},
        AddQuartic[
            kDHHS22,
            {H, Bar @ H, S2, S2},
            GroupInvariant -> (
                del[SU2L @ fund, #1, #2] &
            ),
            SelfConjugate -> False
        ]
    ]
];

TryClosure[
    "Hdag^2 S1dag S2dag + h.c.",
    dims, alpha,
    kDH2S1S2,
    Function[{},
        AddQuartic[
            kDH2S1S2,
            {Bar @ H, Bar @ H, Bar @ S1, Bar @ S2},
            GroupInvariant -> (
                delS2[SU2L, #3, #1, #2] &
            ),
            SelfConjugate -> False
        ]
    ]
];


(* ================================================================== *)
(* T3-E                                                               *)
(* ================================================================== *)

dims = models["E"]["dims"];
alpha = models["E"]["alpha"];

Print[""];
Print["T3-E"];

TryClosure[
    "S1^4 + h.c.",
    dims, alpha,
    kES14,
    Function[{},
        AddQuartic[
            kES14,
            {S1, S1, S1, S1},
            GroupInvariant -> (
                del[SU2L @ S2, #1, #2] *
                del[SU2L @ S2, #3, #4] &
            ),
            SelfConjugate -> False
        ]
    ]
];

TryClosure[
    "H Hdag S1dag^2 + h.c.",
    dims, alpha,
    kEHHS12,
    Function[{},
        AddQuartic[
            kEHHS12,
            {H, Bar @ H, Bar @ S1, Bar @ S1},
            GroupInvariant -> (
                del[SU2L @ fund, #1, #2] *
                del[SU2L @ S2, #3, #4] &
            ),
            SelfConjugate -> False
        ]
    ]
];

TryClosure[
    "Hdag^2 S1 S2 + h.c.",
    dims, alpha,
    kEH2S1S2,
    Function[{},
        AddQuartic[
            kEH2S1S2,
            {Bar @ H, Bar @ H, S1, S2},
            GroupInvariant -> (
                delS2[SU2L, a, #1, #2] *
                T3RGBetaS2TripleInvariant[a, #3, #4] &
            ),
            SelfConjugate -> False
        ]
    ]
];


Print[""];
Print["T3 NEW FIELD-STRUCTURE CLOSURE v2: COMPLETE"];

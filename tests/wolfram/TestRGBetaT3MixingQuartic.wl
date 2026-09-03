(* TestRGBetaT3MixingQuartic.wl
   Validate the topology-defining quartic

       lambdaT3 H H S1 S2^\dagger + h.c.

   for the five d <= 3 T3 classes.

   There are only three distinct SU(2) structures:
       A/D : singlet x triplet
       B/C : doublet x doublet
       E   : triplet x triplet

   The two Higgs fields are identical bosons, so their pair is projected onto
   the symmetric triplet channel.
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

ClearAll[T3MixInvariant, AddT3MixingQuartic, TryMixBeta];

(* A: S1 singlet, S2 triplet
   H x H -> triplet, contracted directly with S2^\dagger. *)
T3MixInvariant[1, 3] :=
    (delS2[SU2L, #4, #1, #2] &);

(* D: S1 triplet, S2 singlet. *)
T3MixInvariant[3, 1] :=
    (delS2[SU2L, #3, #1, #2] &);

(* B/C: S1 and S2 are doublets.
   Couple H H to a triplet and S1 S2^\dagger to the same triplet. *)
T3MixInvariant[2, 2] :=
    (
        delS2[SU2L, a, #1, #2] *
        delS2[SU2L, a, #3, #4] &
    );

(* E: all non-Higgs BSM scalars are triplets.
   3 x 3 contains a triplet through the antisymmetric SU(2) structure
   constant.  The H H pair supplies the other triplet. *)
T3MixInvariant[3, 3] :=
    (
        delS2[SU2L, a, #1, #2] *
        fStruct[SU2L, a, #3, #4] &
    );


AddT3MixingQuartic[dS1_Integer, dS2_Integer] := Module[{inv},
    inv = T3MixInvariant[dS1, dS2];

    AddQuartic[
        lambdaT3,
        {H, H, S1, Bar @ S2},
        GroupInvariant -> inv,
        SelfConjugate -> False
    ];
];


TryMixBeta[] := Module[{result},
    result = CheckAbort[
        Quiet[BetaTerm[lambdaT3, 1]],
        $Aborted
    ];

    Which[
        result === $Aborted,
            Print["  beta(lambdaT3): ABORTED"],
        Head[result] === BetaTerm,
            Print["  beta(lambdaT3): UNEVALUATED"],
        True,
            Print["  beta(lambdaT3): SUCCESS"];
            Print["    ", result]
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


Print[""];
Print["================ RGBeta T3 MIXING QUARTIC ================"];

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
    Print["T3-", label, " ", dims, " alpha=", alpha];

    If[!AssociationQ[build],
        Print["  base build: FAILED"];
        Continue[];
    ];

    Print["  base build: SUCCESS"];

    reg = CheckAbort[
        Quiet[
            AddT3MixingQuartic[dims[[1]], dims[[2]]];
            "Success"
        ],
        $Aborted
    ];

    Print["  lambdaT3 registration: ", reg];

    If[reg === "Success",
        TryMixBeta[];
    ],
    {model, models}
];

Print[""];
Print["RGBETA T3 MIXING QUARTIC CHECK: COMPLETE"];

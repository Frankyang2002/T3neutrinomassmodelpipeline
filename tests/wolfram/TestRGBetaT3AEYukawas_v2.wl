(* TestRGBetaT3AEYukawas_v2.wl
   Corrected A-E T3 Yukawa probe for RGBeta.

   Important translation from the Matchete model:
     - the two T3 Yukawas use opposite heavy-fermion chiralities;
     - a vectorlike heavy fermion is therefore represented in RGBeta by
       two left-Weyl fields F and Fc carrying opposite U(1) charges.

   Project hypercharges:
       Y(S1) = alpha/2
       Y(F)  = (alpha+1)/2
       Y(S2) = (alpha+2)/2

   RGBeta Yukawa choices:
       y1 : Bar[S1]  l  F
       y2 : S2       l  Fc

   Their U(1) sums are identically zero:
       -Y(S1) - 1/2 + Y(F) = 0
        Y(S2) - 1/2 - Y(F) = 0.
*)

ClearAll["Global`*"];
Needs["RGBeta`"];

ClearAll[SU2Rep, GaugeRepFor, T3Invariant, BuildClass];

SU2Rep[1] := None;
SU2Rep[2] := SU2L[fund];
SU2Rep[3] := SU2L[S2];

GaugeRepFor[d_Integer, y_] := If[
    d == 1,
    {U1Y[y]},
    {U1Y[y], SU2Rep[d]}
];

(* Slots are {scalar, lepton-doublet, heavy-fermion}. *)
T3Invariant[2, 1] :=
    (eps[SU2L @ fund, #1, #2] &);

T3Invariant[1, 2] :=
    (eps[SU2L @ fund, #2, #3] &);

T3Invariant[2, 3] :=
    (delS2[SU2L, #3, #1, #2] &);

T3Invariant[3, 2] :=
    (delS2[SU2L, #1, #2, #3] &);


BuildClass[label_String, dims_List, alpha_Integer] := Module[
    {
        dS1, dS2, dF,
        yS1, yS2, yF,
        inv1, inv2,
        reg1, reg2,
        beta1, beta2
    },

    {dS1, dS2, dF} = dims;

    yS1 = alpha/2;
    yF  = (alpha + 1)/2;
    yS2 = (alpha + 2)/2;

    ResetModel[];

    AddGaugeGroup[gY, U1Y, U1];
    AddGaugeGroup[g2, SU2L, SU[2]];

    Dim[gen] = 3;
    Dim[heavy] = 3;

    AddFermion[
        l,
        GaugeRep -> {U1Y[-1/2], SU2L[fund]},
        FlavorIndices -> {gen}
    ];

    (* Two left-Weyl components of the heavy vectorlike fermion.
       When Y(F)=0 they have identical gauge quantum numbers, but keeping
       both labels makes the generic A-E construction uniform. *)
    AddFermion[
        F,
        GaugeRep -> GaugeRepFor[dF, yF],
        FlavorIndices -> {heavy}
    ];

    AddFermion[
        Fc,
        GaugeRep -> GaugeRepFor[dF, -yF],
        FlavorIndices -> {heavy}
    ];

    AddScalar[
        S1,
        GaugeRep -> GaugeRepFor[dS1, yS1]
    ];

    AddScalar[
        S2,
        GaugeRep -> GaugeRepFor[dS2, yS2]
    ];

    inv1 = T3Invariant[dS1, dF];
    inv2 = T3Invariant[dS2, dF];

    reg1 = Check[
        AddYukawa[
            y1,
            {Bar @ S1, l, F},
            GroupInvariant -> inv1,
            CouplingIndices -> ({gen[#2], heavy[#3]} &),
            Chirality -> Right
        ];
        "Success",
        "Failed"
    ];

    reg2 = Check[
        AddYukawa[
            y2,
            {S2, l, Fc},
            GroupInvariant -> inv2,
            CouplingIndices -> ({gen[#2], heavy[#3]} &),
            Chirality -> Right
        ];
        "Success",
        "Failed"
    ];

    beta1 = If[
        reg1 === "Success",
        Check[BetaTerm[y1, 1] // Simplify, $Failed],
        $Failed
    ];

    beta2 = If[
        reg2 === "Success",
        Check[BetaTerm[y2, 1] // Simplify, $Failed],
        $Failed
    ];

    Print[""];
    Print["================ T3-", label, " ================"];
    Print["dims = ", dims, ", alpha = ", alpha];
    Print["Y(S1) = ", yS1, ", Y(S2) = ", yS2, ", Y(F) = ", yF];
    Print[
        "U(1) y1 check = ",
        Simplify[-yS1 - 1/2 + yF]
    ];
    Print[
        "U(1) y2 check = ",
        Simplify[yS2 - 1/2 - yF]
    ];
    Print["y1 registration: ", reg1];
    Print["y2 registration: ", reg2];
    Print[
        "beta(y1): ",
        If[beta1 === $Failed, "Failed", "Success"]
    ];
    Print[
        "beta(y2): ",
        If[beta2 === $Failed, "Failed", "Success"]
    ];
];


BuildClass["A", {1, 3, 2}, 0];
BuildClass["B", {2, 2, 1}, -1];
BuildClass["C", {2, 2, 3}, -1];
BuildClass["D", {3, 1, 2}, -2];
BuildClass["E", {3, 3, 2}, 0];

Print[""];
Print["RGBETA T3 A-E YUKAWA INVARIANT PROBE V2: COMPLETE"];

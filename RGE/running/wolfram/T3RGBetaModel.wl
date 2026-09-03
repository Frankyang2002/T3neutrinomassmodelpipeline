(* T3RGBetaModel.wl
   Reusable RGBeta model builder for the current T3 scope.

   Current supported SU(2) dimensions:
       1 = singlet
       2 = doublet
       3 = triplet

   Project convention:
       Q = T3 + Y
       Y(S1) = alpha/2
       Y(F)  = (alpha + 1)/2
       Y(S2) = (alpha + 2)/2

   The historical T3 paper labels hypercharge by 2Y.

   This file defines the renormalisable gauge/fermion/scalar field content and
   the SM + T3 Yukawa sector in RGBeta.  Scalar-potential registration is kept
   in separate helpers so it can be extended without changing the field logic.
*)

ClearAll[
    T3RGBetaSupportedDimensionQ,
    T3RGBetaRep,
    T3RGBetaGaugeRep,
    T3RGBetaYukawaInvariant,
    T3RGBetaBuild,
    T3RGBetaAddSM,
    T3RGBetaAddHeavyFermion,
    T3RGBetaAddT3Yukawas,
    T3RGBetaAddMasses,
    T3RGBetaPairInvariant,
    T3RGBetaAddBasicQuartics,
    T3RGBetaS2TripleInvariant,
    T3RGBetaMixInvariant,
    T3RGBetaAddMixingQuartic,
    T3RGBetaPortalAdjInvariant,
    T3RGBetaAddRGClosedQuartics,
    T3RGBetaOneLoopBetas
];


T3RGBetaSupportedDimensionQ[d_Integer] := MemberQ[{1, 2, 3}, d];


T3RGBetaRep[1] := None;
T3RGBetaRep[2] := SU2L[fund];
T3RGBetaRep[3] := SU2L[S2];


T3RGBetaGaugeRep[d_Integer, y_] := If[
    d === 1,
    {U1Y[y]},
    {U1Y[y], T3RGBetaRep[d]}
];


(* Slots are always {scalar, lepton doublet, heavy fermion}. *)
T3RGBetaYukawaInvariant[2, 1] :=
    (eps[SU2L @ fund, #1, #2] &);

T3RGBetaYukawaInvariant[1, 2] :=
    (eps[SU2L @ fund, #2, #3] &);

T3RGBetaYukawaInvariant[2, 3] :=
    (delS2[SU2L, #3, #1, #2] &);

T3RGBetaYukawaInvariant[3, 2] :=
    (delS2[SU2L, #1, #2, #3] &);


T3RGBetaAddSM[] := Module[{},
    AddGaugeGroup[gY, U1Y, U1];
    AddGaugeGroup[g2, SU2L, SU[2]];
    AddGaugeGroup[g3, SU3c, SU[3]];

    Dim[gen] = 3;

    AddFermion[
        q,
        GaugeRep -> {
            U1Y[1/6],
            SU2L[fund],
            SU3c[fund]
        },
        FlavorIndices -> {gen}
    ];

    AddFermion[
        u,
        GaugeRep -> {
            U1Y[-2/3],
            Bar @ SU3c[fund]
        },
        FlavorIndices -> {gen}
    ];

    AddFermion[
        d,
        GaugeRep -> {
            U1Y[1/3],
            Bar @ SU3c[fund]
        },
        FlavorIndices -> {gen}
    ];

    AddFermion[
        l,
        GaugeRep -> {
            U1Y[-1/2],
            SU2L[fund]
        },
        FlavorIndices -> {gen}
    ];

    AddFermion[
        e,
        GaugeRep -> {
            U1Y[1]
        },
        FlavorIndices -> {gen}
    ];

    AddScalar[
        H,
        GaugeRep -> {
            U1Y[1/2],
            SU2L[fund]
        }
    ];

    (* Standard-model Yukawa structures. *)
    AddYukawa[
        yu,
        {H, q, u},
        GroupInvariant -> (
            del[SU3c @ fund, #2, #3] *
            eps[SU2L @ fund, #1, #2] &
        ),
        CouplingIndices -> ({gen[#2], gen[#3]} &),
        Chirality -> Right
    ];

    AddYukawa[
        yd,
        {Bar @ H, q, d},
        GroupInvariant -> (
            del[SU3c @ fund, #2, #3] *
            del[SU2L @ fund, #1, #2] &
        ),
        CouplingIndices -> ({gen[#2], gen[#3]} &),
        Chirality -> Right
    ];

    AddYukawa[
        ye,
        {Bar @ H, l, e},
        GroupInvariant -> (
            del[SU2L @ fund, #1, #2] &
        ),
        CouplingIndices -> ({gen[#2], gen[#3]} &),
        Chirality -> Right
    ];
];


T3RGBetaAddHeavyFermion[dF_Integer, yF_] := Module[{},
    Dim[heavy] = 3;

    If[yF === 0,
        (* Real/neutral representation: one heavy Weyl field with a Majorana
           mass.  This is the T3-B/T3-C situation at alpha=-1. *)
        AddFermion[
            F,
            GaugeRep -> T3RGBetaGaugeRep[dF, 0],
            FlavorIndices -> {heavy}
        ],

        (* Vectorlike case: two left-Weyl fields in conjugate gauge
           representations. *)
        AddFermion[
            F,
            GaugeRep -> T3RGBetaGaugeRep[dF, yF],
            FlavorIndices -> {heavy}
        ];

        AddFermion[
            Fc,
            GaugeRep -> T3RGBetaGaugeRep[dF, -yF],
            FlavorIndices -> {heavy}
        ];
    ];
];


T3RGBetaAddT3Yukawas[
    dS1_Integer,
    dS2_Integer,
    dF_Integer,
    yF_
] := Module[
    {inv1, inv2, secondHeavy},

    inv1 = T3RGBetaYukawaInvariant[dS1, dF];
    inv2 = T3RGBetaYukawaInvariant[dS2, dF];

    secondHeavy = If[yF === 0, F, Fc];

    AddYukawa[
        y1,
        {Bar @ S1, l, F},
        GroupInvariant -> inv1,
        CouplingIndices -> ({gen[#2], heavy[#3]} &),
        Chirality -> Right
    ];

    AddYukawa[
        y2,
        {S2, l, secondHeavy},
        GroupInvariant -> inv2,
        CouplingIndices -> ({gen[#2], heavy[#3]} &),
        Chirality -> Right
    ];
];


T3RGBetaAddMasses[
    dS1_Integer,
    dS2_Integer,
    dF_Integer,
    yF_
] := Module[
    {scalarInvariant1, scalarInvariant2, fermionInvariant},

    scalarInvariant1 = If[
        dS1 === 1,
        (1 &),
        (del[T3RGBetaRep[dS1], #1, #2] &)
    ];

    scalarInvariant2 = If[
        dS2 === 1,
        (1 &),
        (del[T3RGBetaRep[dS2], #1, #2] &)
    ];

    AddScalarMass[
        mS1Sq,
        {Bar @ S1, S1},
        GroupInvariant -> scalarInvariant1
    ];

    AddScalarMass[
        mS2Sq,
        {Bar @ S2, S2},
        GroupInvariant -> scalarInvariant2
    ];

    If[yF === 0,
        fermionInvariant = If[
            dF === 1,
            (1 &),
            (del[T3RGBetaRep[dF], #1, #2] &)
        ];

        AddFermionMass[
            MF,
            {F, F},
            GroupInvariant -> fermionInvariant,
            MassIndices -> ({heavy[#1], heavy[#2]} &),
            Chirality -> Right
        ],

        fermionInvariant = If[
            dF === 1,
            (1 &),
            (del[T3RGBetaRep[dF], #1, #2] &)
        ];

        AddFermionMass[
            MF,
            {F, Fc},
            GroupInvariant -> fermionInvariant,
            MassIndices -> ({heavy[#1], heavy[#2]} &),
            Chirality -> Right
        ];
    ];
];


(* ---------------------------------------------------------------------- *)
(* Scalar potential                                                        *)
(* ---------------------------------------------------------------------- *)

T3RGBetaPairInvariant[1, i_, j_] := 1;
T3RGBetaPairInvariant[2, i_, j_] := del[SU2L @ fund, i, j];
T3RGBetaPairInvariant[3, i_, j_] := del[SU2L @ S2, i, j];


T3RGBetaAddBasicQuartics[
    dS1_Integer,
    dS2_Integer
] := Module[{},
    (* Coupling normalization follows the current Matchete Lagrangian:
         lambdaH /2 (H†H)^2
         lambdaS1/2 (S1†S1)^2
         lambdaS2/2 (S2†S2)^2
         lambdaH1 (H†H)(S1†S1)
         lambdaH2 (H†H)(S2†S2)
         lambda12 (S1†S1)(S2†S2)
    *)
    SetReal[
        lambdaH,
        lambdaS1,
        lambdaS2,
        lambdaH1,
        lambdaH2,
        lambda12
    ];

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
            T3RGBetaPairInvariant[dS1, #1, #2] *
            T3RGBetaPairInvariant[dS1, #3, #4] / 2 &
        )
    ];

    AddQuartic[
        lambdaS2,
        {Bar @ S2, S2, Bar @ S2, S2},
        GroupInvariant -> (
            T3RGBetaPairInvariant[dS2, #1, #2] *
            T3RGBetaPairInvariant[dS2, #3, #4] / 2 &
        )
    ];

    AddQuartic[
        lambdaH1,
        {Bar @ H, H, Bar @ S1, S1},
        GroupInvariant -> (
            del[SU2L @ fund, #1, #2] *
            T3RGBetaPairInvariant[dS1, #3, #4] &
        )
    ];

    AddQuartic[
        lambdaH2,
        {Bar @ H, H, Bar @ S2, S2},
        GroupInvariant -> (
            del[SU2L @ fund, #1, #2] *
            T3RGBetaPairInvariant[dS2, #3, #4] &
        )
    ];

    AddQuartic[
        lambda12,
        {Bar @ S1, S1, Bar @ S2, S2},
        GroupInvariant -> (
            T3RGBetaPairInvariant[dS1, #1, #2] *
            T3RGBetaPairInvariant[dS2, #3, #4] &
        )
    ];
];


(* Three-triplet antisymmetric invariant written entirely in the S2 basis.
   This is required for T3-E.  Using fStruct would introduce SU2L[adj]
   indices while our triplet fields are registered as SU2L[S2]. *)
T3RGBetaS2TripleInvariant[a_, b_, c_] := Module[
    {i, j, k, l, m, n},
    delS2[SU2L, a, i, j] *
    delS2[SU2L, b, k, l] *
    delS2[SU2L, c, m, n] *
    eps[SU2L @ fund, i, k] *
    eps[SU2L @ fund, j, m] *
    eps[SU2L @ fund, l, n]
];


(* Topology quartic:
       lambdaT3 H H S1 S2† + h.c.

   The two identical Higgs fields are in the symmetric triplet channel.
*)
T3RGBetaMixInvariant[1, 3] :=
    (delS2[SU2L, #4, #1, #2] &);

T3RGBetaMixInvariant[3, 1] :=
    (delS2[SU2L, #3, #1, #2] &);

T3RGBetaMixInvariant[2, 2] :=
    (
        delS2[SU2L, a, #1, #2] *
        delS2[SU2L, a, #3, #4] &
    );

T3RGBetaMixInvariant[3, 3] :=
    (
        delS2[SU2L, a, #1, #2] *
        T3RGBetaS2TripleInvariant[a, #3, #4] &
    );


T3RGBetaAddMixingQuartic[
    dS1_Integer,
    dS2_Integer
] := AddQuartic[
    lambdaT3,
    {H, H, S1, Bar @ S2},
    GroupInvariant -> T3RGBetaMixInvariant[dS1, dS2],
    SelfConjugate -> False
];



(* ---------------------------------------------------------------------- *)
(* Additional quartics required by one-loop RG closure                    *)
(* ---------------------------------------------------------------------- *)

T3RGBetaPortalAdjInvariant[2] :=
    (
        tGen[SU2L @ fund, A, #1, #2] *
        tGen[SU2L @ fund, A, #3, #4] &
    );

T3RGBetaPortalAdjInvariant[3] :=
    (
        tGen[SU2L @ fund, A, #1, #2] *
        tGen[SU2L @ S2, A, #3, #4] &
    );


T3RGBetaAddRGClosedQuartics[
    dS1_Integer,
    dS2_Integer,
    dF_Integer,
    alpha_Integer
] := Module[{},
    (* -------------------------------------------------------------- *)
    (* Extra Higgs-BSM portal contractions                            *)
    (* -------------------------------------------------------------- *)

    If[dS1 > 1,
        SetReal[lambdaH1Adj];

        AddQuartic[
            lambdaH1Adj,
            {Bar @ H, H, Bar @ S1, S1},
            GroupInvariant -> T3RGBetaPortalAdjInvariant[dS1]
        ];
    ];

    If[dS2 > 1,
        SetReal[lambdaH2Adj];

        AddQuartic[
            lambdaH2Adj,
            {Bar @ H, H, Bar @ S2, S2},
            GroupInvariant -> T3RGBetaPortalAdjInvariant[dS2]
        ];
    ];


    (* -------------------------------------------------------------- *)
    (* Second self quartic for complex triplets                       *)
    (* -------------------------------------------------------------- *)

    If[dS1 == 3,
        SetReal[lambdaS1Adj];

        AddQuartic[
            lambdaS1Adj,
            {Bar @ S1, S1, Bar @ S1, S1},
            GroupInvariant -> (
                tGen[SU2L @ S2, A, #1, #2] *
                tGen[SU2L @ S2, A, #3, #4] / 2 &
            )
        ];
    ];

    If[dS2 == 3,
        SetReal[lambdaS2Adj];

        AddQuartic[
            lambdaS2Adj,
            {Bar @ S2, S2, Bar @ S2, S2},
            GroupInvariant -> (
                tGen[SU2L @ S2, A, #1, #2] *
                tGen[SU2L @ S2, A, #3, #4] / 2 &
            )
        ];
    ];


    (* -------------------------------------------------------------- *)
    (* Additional S1-S2 contractions                                  *)
    (* -------------------------------------------------------------- *)

    If[dS1 > 1 && dS2 > 1,
        SetReal[lambda12Adj];

        AddQuartic[
            lambda12Adj,
            {Bar @ S1, S1, Bar @ S2, S2},
            GroupInvariant -> (
                tGen[T3RGBetaRep[dS1], A, #1, #2] *
                tGen[T3RGBetaRep[dS2], A, #3, #4] &
            )
        ];
    ];

    (* T3-E has a third independent mixed contraction. *)
    If[dS1 == 3 && dS2 == 3,
        SetReal[lambda12Cross];

        AddQuartic[
            lambda12Cross,
            {Bar @ S1, S1, Bar @ S2, S2},
            GroupInvariant -> (
                del[SU2L @ S2, #1, #3] *
                del[SU2L @ S2, #2, #4] +
                del[SU2L @ S2, #1, #4] *
                del[SU2L @ S2, #2, #3] &
            )
        ];
    ];


    (* -------------------------------------------------------------- *)
    (* New field structures required only for T3-B / T3-C            *)
    (* -------------------------------------------------------------- *)

    If[dS1 == 2 && dS2 == 2 && alpha == -1,
        AddQuartic[
            lambdaHHdagS2S2,
            {Bar @ H, Bar @ H, S2, S2},
            GroupInvariant -> (
                delS2[SU2L, a, #1, #2] *
                delS2[SU2L, a, #3, #4] &
            ),
            SelfConjugate -> False
        ];

        AddQuartic[
            lambdaHHdagS1barS1bar,
            {Bar @ H, Bar @ H, Bar @ S1, Bar @ S1},
            GroupInvariant -> (
                delS2[SU2L, a, #1, #2] *
                delS2[SU2L, a, #3, #4] &
            ),
            SelfConjugate -> False
        ];

        AddQuartic[
            lambdaS1bar2S2bar2,
            {Bar @ S1, Bar @ S1, Bar @ S2, Bar @ S2},
            GroupInvariant -> (
                delS2[SU2L, a, #1, #2] *
                delS2[SU2L, a, #3, #4] &
            ),
            SelfConjugate -> False
        ];

        (* Generated only after the first B/C closure enlargement. *)

        AddQuartic[
            lambdaS1barS2S2bar2,
            {Bar @ S1, S2, Bar @ S2, Bar @ S2},
            GroupInvariant -> (
                delS2[SU2L, a, #1, #2] *
                delS2[SU2L, a, #3, #4] &
            ),
            SelfConjugate -> False
        ];

        AddQuartic[
            lambdaS1S1bar2S2bar,
            {S1, Bar @ S1, Bar @ S1, Bar @ S2},
            GroupInvariant -> (
                delS2[SU2L, a, #1, #4] *
                delS2[SU2L, a, #2, #3] &
            ),
            SelfConjugate -> False
        ];

        AddQuartic[
            lambdaHHdagS1barS2barCross,
            {H, Bar @ H, Bar @ S1, Bar @ S2},
            GroupInvariant -> (
                eps[SU2L @ fund, #1, #3] *
                eps[SU2L @ fund, #2, #4] &
            ),
            SelfConjugate -> False
        ];
    ];
];



T3RGBetaBuild[
    dS1_Integer,
    dS2_Integer,
    dF_Integer,
    alpha_Integer
] := Module[
    {yS1, yS2, yF},

    If[
        !And @@ (T3RGBetaSupportedDimensionQ /@ {dS1, dS2, dF}),
        Return[$Failed]
    ];

    yS1 = alpha/2;
    yF  = (alpha + 1)/2;
    yS2 = (alpha + 2)/2;

    ResetModel[];

    T3RGBetaAddSM[];

    T3RGBetaAddHeavyFermion[dF, yF];

    AddScalar[
        S1,
        GaugeRep -> T3RGBetaGaugeRep[dS1, yS1]
    ];

    AddScalar[
        S2,
        GaugeRep -> T3RGBetaGaugeRep[dS2, yS2]
    ];

    T3RGBetaAddT3Yukawas[dS1, dS2, dF, yF];
    T3RGBetaAddMasses[dS1, dS2, dF, yF];

    T3RGBetaAddBasicQuartics[dS1, dS2];
    T3RGBetaAddMixingQuartic[dS1, dS2];
    T3RGBetaAddRGClosedQuartics[dS1, dS2, dF, alpha];

    T3RGBetaLastBuild = <|
        "dS1" -> dS1,
        "dS2" -> dS2,
        "dF" -> dF,
        "alpha" -> alpha
    |>;

    <|
        "dS1" -> dS1,
        "dS2" -> dS2,
        "dF" -> dF,
        "alpha" -> alpha,
        "YS1" -> yS1,
        "YS2" -> yS2,
        "YF" -> yF,
        "MajoranaFermion" -> (yF === 0)
    |>
];


T3RGBetaOneLoopBetas[] := Module[
    {result, meta, dS1, dS2, alpha},

    meta = T3RGBetaLastBuild;
    dS1 = meta["dS1"];
    dS2 = meta["dS2"];
    alpha = meta["alpha"];

    result = <|
        "gY" -> Quiet[BetaTerm[gY, 1]],
        "g2" -> Quiet[BetaTerm[g2, 1]],
        "g3" -> Quiet[BetaTerm[g3, 1]],
        "yu" -> Quiet[BetaTerm[yu, 1]],
        "yd" -> Quiet[BetaTerm[yd, 1]],
        "ye" -> Quiet[BetaTerm[ye, 1]],
        "y1" -> Quiet[BetaTerm[y1, 1]],
        "y2" -> Quiet[BetaTerm[y2, 1]],
        "MF" -> Quiet[BetaTerm[MF, 1]],
        "mS1Sq" -> Quiet[BetaTerm[mS1Sq, 1]],
        "mS2Sq" -> Quiet[BetaTerm[mS2Sq, 1]],
        "lambdaH" -> Quiet[BetaTerm[lambdaH, 1]],
        "lambdaS1" -> Quiet[BetaTerm[lambdaS1, 1]],
        "lambdaS2" -> Quiet[BetaTerm[lambdaS2, 1]],
        "lambdaH1" -> Quiet[BetaTerm[lambdaH1, 1]],
        "lambdaH2" -> Quiet[BetaTerm[lambdaH2, 1]],
        "lambda12" -> Quiet[BetaTerm[lambda12, 1]],
        "lambdaT3" -> Quiet[BetaTerm[lambdaT3, 1]]
    |>;

    If[dS1 > 1,
        AssociateTo[
            result,
            "lambdaH1Adj" -> Quiet[BetaTerm[lambdaH1Adj, 1]]
        ];
    ];

    If[dS2 > 1,
        AssociateTo[
            result,
            "lambdaH2Adj" -> Quiet[BetaTerm[lambdaH2Adj, 1]]
        ];
    ];

    If[dS1 == 3,
        AssociateTo[
            result,
            "lambdaS1Adj" -> Quiet[BetaTerm[lambdaS1Adj, 1]]
        ];
    ];

    If[dS2 == 3,
        AssociateTo[
            result,
            "lambdaS2Adj" -> Quiet[BetaTerm[lambdaS2Adj, 1]]
        ];
    ];

    If[dS1 > 1 && dS2 > 1,
        AssociateTo[
            result,
            "lambda12Adj" -> Quiet[BetaTerm[lambda12Adj, 1]]
        ];
    ];

    If[dS1 == 3 && dS2 == 3,
        AssociateTo[
            result,
            "lambda12Cross" -> Quiet[BetaTerm[lambda12Cross, 1]]
        ];
    ];

    If[dS1 == 2 && dS2 == 2 && alpha == -1,
        AssociateTo[
            result,
            <|
                "lambdaHHdagS2S2" ->
                    Quiet[BetaTerm[lambdaHHdagS2S2, 1]],
                "lambdaHHdagS1barS1bar" ->
                    Quiet[BetaTerm[lambdaHHdagS1barS1bar, 1]],
                "lambdaS1bar2S2bar2" ->
                    Quiet[BetaTerm[lambdaS1bar2S2bar2, 1]],
                "lambdaS1barS2S2bar2" ->
                    Quiet[BetaTerm[lambdaS1barS2S2bar2, 1]],
                "lambdaS1S1bar2S2bar" ->
                    Quiet[BetaTerm[lambdaS1S1bar2S2bar, 1]],
                "lambdaHHdagS1barS2barCross" ->
                    Quiet[BetaTerm[lambdaHHdagS1barS2barCross, 1]]
            |>
        ];
    ];

    result
];

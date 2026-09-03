(* ProbeRGBetaMaYukawaOrientation.wl

   Diagnose the SU(2) invariant used for the Ma/scotogenic Yukawa h eta L N.

   For two SU(2) doublets, the singlet contraction is antisymmetric:
       2 x 2 -> 1
   and is represented by eps[SU2L@fund, i, j].

   We compare:
     1. del orientation  -- the contraction previously used.
     2. eps orientation  -- the physical antisymmetric doublet singlet.

   We only print the h-Y_e mixed contributions to beta_lambda3 and
   beta_lambda4, because these terms distinguish the two orientations.

   Trusted Ma target:
       beta_lambda3 |_(h Ye) = -4 T_nue
       beta_lambda4 |_(h Ye) = +4 T_nue

   where T_nue = Tr[h h^\[Dagger] Ye Ye^\[Dagger]] up to matrix orientation.
*)

ClearAll["Global`*"];
Needs["RGBeta`"];

ClearAll[BuildAndProbe];

BuildAndProbe[label_, hInvariant_] := Module[
    {betaQ, mixed3, mixed4, zeroRules},

    ResetModel[];

    AddGaugeGroup[gY, U1Y, U1];
    AddGaugeGroup[g2, SU2L, SU[2]];
    AddGaugeGroup[g3, SU3c, SU[3]];

    Dim[gen] = 3;
    Dim[heavy] = 3;

    AddFermion[
        q,
        GaugeRep -> {U1Y[1/6], SU2L[fund], SU3c[fund]},
        FlavorIndices -> {gen}
    ];

    AddFermion[
        u,
        GaugeRep -> {U1Y[-2/3], Bar @ SU3c[fund]},
        FlavorIndices -> {gen}
    ];

    AddFermion[
        d,
        GaugeRep -> {U1Y[1/3], Bar @ SU3c[fund]},
        FlavorIndices -> {gen}
    ];

    AddFermion[
        l,
        GaugeRep -> {U1Y[-1/2], SU2L[fund]},
        FlavorIndices -> {gen}
    ];

    AddFermion[
        e,
        GaugeRep -> {U1Y[1]},
        FlavorIndices -> {gen}
    ];

    AddFermion[
        N,
        FlavorIndices -> {heavy}
    ];

    AddScalar[
        H,
        GaugeRep -> {U1Y[1/2], SU2L[fund]}
    ];

    AddScalar[
        eta,
        GaugeRep -> {U1Y[1/2], SU2L[fund]}
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

    AddYukawa[
        h,
        {eta, l, N},
        GroupInvariant -> hInvariant,
        CouplingIndices -> ({gen[#2], heavy[#3]} &),
        Chirality -> Right
    ];

    AddQuartic[
        lambda1,
        {Bar @ H, H, Bar @ H, H},
        GroupInvariant -> (
            del[SU2L @ fund, #1, #2]
            del[SU2L @ fund, #3, #4] / 2 &
        )
    ];

    AddQuartic[
        lambda2,
        {Bar @ eta, eta, Bar @ eta, eta},
        GroupInvariant -> (
            del[SU2L @ fund, #1, #2]
            del[SU2L @ fund, #3, #4] / 2 &
        )
    ];

    AddQuartic[
        lambda3,
        {Bar @ H, H, Bar @ eta, eta},
        GroupInvariant -> (
            del[SU2L @ fund, #1, #2]
            del[SU2L @ fund, #3, #4] &
        )
    ];

    AddQuartic[
        lambda4,
        {Bar @ H, eta, Bar @ eta, H},
        GroupInvariant -> (
            del[SU2L @ fund, #1, #2]
            del[SU2L @ fund, #3, #4] &
        )
    ];

    AddQuartic[
        lambda5,
        {Bar @ H, eta, Bar @ H, eta},
        GroupInvariant -> (
            del[SU2L @ fund, #1, #2]
            del[SU2L @ fund, #3, #4] / 2 &
        ),
        SelfConjugate -> False
    ];

    betaQ = BetaTerm[Quartic, 1] // Simplify;

    (* Remove gauge and quartic terms so only Yukawa-dependent pieces remain. *)
    zeroRules = {
        gY -> 0, g2 -> 0, g3 -> 0,
        lambda1 -> 0, lambda2 -> 0, lambda3 -> 0,
        lambda4 -> 0, lambda5 -> 0,
        Bar[lambda5] -> 0
    };

    mixed3 = Expand[betaQ[lambda3] /. zeroRules] // Simplify;
    mixed4 = Expand[betaQ[lambda4] /. zeroRules] // Simplify;

    Print[""];
    Print["================ ", label, " ================"];
    Print["beta_lambda3 Yukawa-only:"];
    Print[mixed3];
    Print[""];
    Print["beta_lambda4 Yukawa-only:"];
    Print[mixed4];
];

BuildAndProbe[
    "DEL ORIENTATION",
    (del[SU2L @ fund, #1, #2] &)
];

BuildAndProbe[
    "EPS ORIENTATION",
    (eps[SU2L @ fund, #1, #2] &)
];

Print[""];
Print["Trusted Ma mixed target:"];
Print["beta_lambda3 -> -4 Tr[h h^dagger Ye Ye^dagger]"];
Print["beta_lambda4 -> +4 Tr[h h^dagger Ye Ye^dagger]"];
Print[""];
Print["RGBETA MA YUKAWA ORIENTATION PROBE: COMPLETE"];

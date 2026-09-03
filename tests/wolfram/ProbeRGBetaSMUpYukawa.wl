(* ProbeRGBetaSMUpYukawa.wl
   Isolate the RGBeta convention for the SM up-type Yukawa.

   We vary:
     1. the anti-fundamental colour representation spelling;
     2. the SU(3) invariant spelling.

   The physical hypercharge assignment is fixed:
       q : (3,2)_(+1/6)
       u^c : (bar 3,1)_(-2/3)
       H : (1,2)_(+1/2)

   Hence H q u^c is U(1)-invariant and the SU(2) contraction is epsilon.
*)

ClearAll["Global`*"];
Needs["RGBeta`"];

ClearAll[BuildProbe];

BuildProbe[label_String, uRep_, colorInvariant_] := Module[{reg, beta},

    ResetModel[];

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
            uRep
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

    reg = CheckAbort[
        Quiet[
            AddYukawa[
                yu,
                {H, q, u},
                GroupInvariant -> (
                    eps[SU2L @ fund, #1, #2] *
                    colorInvariant[#3, #4] &
                ),
                CouplingIndices -> ({gen[#2], gen[#3]} &),
                Chirality -> Right
            ];
            "Success"
        ],
        $Aborted
    ];

    beta = If[
        reg === "Success",
        CheckAbort[
            Quiet[BetaTerm[yu, 1]],
            $Aborted
        ],
        "NotRun"
    ];

    Print[""];
    Print[label];
    Print["  registration: ", reg];
    Print[
        "  beta: ",
        Which[
            beta === $Aborted, "ABORTED",
            beta === "NotRun", "NOT RUN",
            True, "SUCCESS"
        ]
    ];

    If[beta =!= $Aborted && beta =!= "NotRun",
        Print["  ", beta]
    ];
];


(* Variant 1: the spelling used in our first production builder. *)
BuildProbe[
    "V1: SU3c[Bar@fund] + del[SU3c@fund]",
    SU3c[Bar @ fund],
    Function[{a, b}, del[SU3c @ fund, a, b]]
];

(* Variant 2: conjugate the full representation rather than the rep label. *)
BuildProbe[
    "V2: Bar[SU3c@fund] + del[SU3c@fund]",
    Bar[SU3c @ fund],
    Function[{a, b}, del[SU3c @ fund, a, b]]
];

(* Variant 3: conjugate the full rep and write the delta in that rep space. *)
BuildProbe[
    "V3: Bar[SU3c@fund] + del[Bar[SU3c@fund]]",
    Bar[SU3c @ fund],
    Function[{a, b}, del[Bar[SU3c @ fund], a, b]]
];

Print[""];
Print["RGBETA SM UP-YUKAWA PROBE: COMPLETE"];

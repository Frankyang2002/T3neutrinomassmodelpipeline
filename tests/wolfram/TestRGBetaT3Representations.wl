(* TestRGBetaT3Representations.wl
   Capability probe for RGBeta's SU(2) representation support.

   The generalized T3 project needs arbitrary SU(2) dimensions, including
   examples such as:
       (dS1,dS2,dF) = (3,3,4), (3,5,4), (4,6,5).

   RGBeta documents several built-in common representations:
       fund, adj, S2, A2.

   For SU(2):
       fund -> d=2
       adj  -> d=3
       S2   -> d=3
       A2   -> d=1

   This test checks those built-ins and then probes custom labels R4, R5, R6
   to determine whether RGBeta will accept arbitrary representation symbols
   after supplying their elementary group invariants by hand.
*)

ClearAll["Global`*"];
Needs["RGBeta`"];
ResetModel[];

AddGaugeGroup[g2, SU2L, SU[2]];

Print[""];
Print["================ RGBeta T3 REPRESENTATION PROBE ================"];

builtins = {
    "fund" -> SU2L[fund],
    "adj"  -> SU2L[adj],
    "S2"   -> SU2L[S2],
    "A2"   -> SU2L[A2]
};

Print[""];
Print["Built-in SU(2) representations:"];
Do[
    Print[
        name,
        ": rep=", rep,
        ", RepresentationCheck=", Quiet @ Check[RepresentationCheck[rep], $Failed],
        ", Dim=", Quiet @ Check[Dim[rep], $Failed],
        ", C2=", Quiet @ Check[Casimir2[rep], $Failed],
        ", T=", Quiet @ Check[TraceNormalization[rep], $Failed]
    ],
    {entry, builtins},
    {name, {First[entry]}},
    {rep, {Last[entry]}}
];

(* ------------------------------------------------------------------------- *)
(* Custom SU(2) representation labels                                        *)
(* ------------------------------------------------------------------------- *)

(* SU(2) irrep formulas for dimension d:
       C2 = (d^2 - 1)/4
       T  = d(d^2 - 1)/12
*)
Clear[R4, R5, R6];

customData = {
    {R4, 4},
    {R5, 5},
    {R6, 6}
};

Do[
    rep = SU2L[label];

    Dim[rep] = d;
    Casimir2[rep] = (d^2 - 1)/4;
    TraceNormalization[rep] = d (d^2 - 1)/12;

    Print[""];
    Print[
        "Custom d=", d,
        ": rep=", rep,
        ", RepresentationCheck=", Quiet @ Check[RepresentationCheck[rep], $Failed],
        ", Dim=", Dim[rep],
        ", C2=", Casimir2[rep],
        ", T=", TraceNormalization[rep]
    ];

    (* Try registering a scalar and fermion separately. *)
    scalarResult = Quiet @ Check[
        AddScalar[
            Symbol["ProbeScalar" <> ToString[d]],
            GaugeRep -> {rep}
        ];
        "Success",
        "Failed"
    ];

    fermionResult = Quiet @ Check[
        AddFermion[
            Symbol["ProbeFermion" <> ToString[d]],
            GaugeRep -> {rep}
        ];
        "Success",
        "Failed"
    ];

    Print[
        "  AddScalar=", scalarResult,
        ", AddFermion=", fermionResult
    ],
    {item, customData},
    {label, {item[[1]]}},
    {d, {item[[2]]}}
];

Print[""];
Print["Interpretation:"];
Print["  PASS requires RGBeta to accept d=4,5,6 representations, not merely"];
Print["  allow Dim/Casimir assignments. For production T3 use we ultimately"];
Print["  also need invariant tensors coupling these higher reps."];
Print[""];
Print["RGBETA T3 REPRESENTATION PROBE: COMPLETE"];

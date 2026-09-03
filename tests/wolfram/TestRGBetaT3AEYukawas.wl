(* TestRGBetaT3AEYukawas.wl
   Validate RGBeta SU(2) invariant choices for the five T3 classes A-E.

   Current project scope:
       only SU(2) dimensions 1, 2, 3.

   RGBeta representation map:
       d=1 : no non-trivial SU(2) representation
       d=2 : SU2L[fund]
       d=3 : SU2L[S2]

   Relevant SU(2) singlet tensors:
       eps[SU2L@fund, i, j]
           antisymmetric 2 x 2 -> 1

       delS2[SU2L, a, i, j]
           symmetric triplet-doublet-doublet invariant
           3 x 2 x 2 -> 1

   This probe only asks whether RGBeta accepts the T3 Yukawa structures.
   It does NOT yet build the full scalar potential.

   Hypercharge convention used internally by this project:
       Q = T3 + Y
       Y(S1) = alpha/2
       Y(F)  = (alpha+1)/2
       Y(S2) = (alpha+2)/2

   The historical paper labels hypercharge by 2Y, hence its integer alpha.
*)

ClearAll["Global`*"];
Needs["RGBeta`"];

ClearAll[SU2Rep, GaugeRepFor, BuildClass];

SU2Rep[1] := None;
SU2Rep[2] := SU2L[fund];
SU2Rep[3] := SU2L[S2];

GaugeRepFor[d_Integer, y_] := If[
    d == 1,
    {U1Y[y]},
    {U1Y[y], SU2Rep[d]}
];

(* Return the SU(2) invariant for L x F x S -> 1.
   Slot convention below is {scalar, lepton, fermion}. *)
ClearAll[T3YukawaInvariant];

(* scalar doublet, fermion singlet: S x L -> 1 *)
T3YukawaInvariant[2, 1] :=
    (eps[SU2L @ fund, #1, #2] &);

(* scalar singlet, fermion doublet: L x F -> 1 *)
T3YukawaInvariant[1, 2] :=
    (eps[SU2L @ fund, #2, #3] &);

(* scalar doublet, fermion triplet: S x L x F -> 1 *)
T3YukawaInvariant[2, 3] :=
    (delS2[SU2L, #3, #1, #2] &);

(* scalar triplet, fermion doublet: S x L x F -> 1 *)
T3YukawaInvariant[3, 2] :=
    (delS2[SU2L, #1, #2, #3] &);


BuildClass[label_String, dims_List, alpha_Integer] := Module[
    {
        dS1, dS2, dF,
        yS1, yS2, yF,
        inv1, inv2,
        result1, result2,
        betaY
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

    AddFermion[
        F,
        GaugeRep -> GaugeRepFor[dF, -yF],
        FlavorIndices -> {heavy}
    ];

    (* We use the conjugate-Weyl field convention for F here, hence -Y(F).
       The purpose of this test is the SU(2) invariant structure. *)

    AddScalar[
        S1,
        GaugeRep -> GaugeRepFor[dS1, yS1]
    ];

    AddScalar[
        S2,
        GaugeRep -> GaugeRepFor[dS2, yS2]
    ];

    inv1 = T3YukawaInvariant[dS1, dF];
    inv2 = T3YukawaInvariant[dS2, dF];

    result1 = Quiet @ Check[
        AddYukawa[
            y1,
            {S1, l, F},
            GroupInvariant -> inv1,
            CouplingIndices -> ({gen[#2], heavy[#3]} &),
            Chirality -> Right,
            CheckInvariance -> True
        ];
        If[KeyExistsQ[$yukawas, y1], "Success", "Rejected"],
        "Failed"
    ];

    result2 = Quiet @ Check[
        AddYukawa[
            y2,
            {Bar @ S2, l, F},
            GroupInvariant -> inv2,
            CouplingIndices -> ({gen[#2], heavy[#3]} &),
            Chirality -> Right,
            CheckInvariance -> True
        ];
        If[KeyExistsQ[$yukawas, y2], "Success", "Rejected"],
        "Failed"
    ];

    Print[""];
    Print["================ T3-", label, " ================"];
    Print["dims = ", dims, ", alpha = ", alpha];
    Print["Y(S1) = ", yS1, ", Y(S2) = ", yS2, ", Y(F) = ", yF];
    Print["y1 registration: ", result1];
    Print["y2 registration: ", result2];

    If[result1 === "Success" && result2 === "Success",
        betaY = Quiet @ Check[BetaTerm[Yukawa, 1] // Simplify, $Failed];
        Print["one-loop Yukawa beta generation: ",
            If[betaY === $Failed, "Failed", "Success"]
        ],
        Print["one-loop Yukawa beta generation: NotRun"]
    ];
];


(* Historical benchmark alpha choices already used by the pipeline. *)
BuildClass["A", {1, 3, 2}, 0];
BuildClass["B", {2, 2, 1}, -1];
BuildClass["C", {2, 2, 3}, -1];
BuildClass["D", {3, 1, 2}, -2];
BuildClass["E", {3, 3, 2}, 0];

Print[""];
Print["RGBETA T3 A-E YUKAWA INVARIANT PROBE: COMPLETE"];

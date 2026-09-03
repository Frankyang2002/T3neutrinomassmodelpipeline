(* TestRGBetaSMCouplings.wl
   RGBeta migration diagnostic for the renormalisable SM matter sector.

   Purpose:
     1. Verify that the official RGBeta SM model definition works headlessly.
     2. Print the one-loop Yukawa, Higgs-quartic and Higgs-mass beta terms.
     3. Establish the precise RGBeta normalization before the T3/Ma model is added.

   This deliberately follows RGBeta's own Sample_models.nb conventions.
*)

ClearAll["Global`*"];
Needs["RGBeta`"];
ResetModel[];

(* ------------------------------------------------------------------------- *)
(* Gauge groups                                                              *)
(* ------------------------------------------------------------------------- *)

AddGaugeGroup[gY, U1Y, U1];
AddGaugeGroup[g2, SU2L, SU[2]];
AddGaugeGroup[g3, SU3c, SU[3]];

(* ------------------------------------------------------------------------- *)
(* Standard Model left-handed Weyl fermions                                  *)
(* ------------------------------------------------------------------------- *)

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

Dim[gen] = 3;

(* ------------------------------------------------------------------------- *)
(* Higgs                                                                     *)
(* ------------------------------------------------------------------------- *)

AddScalar[
    H,
    GaugeRep -> {
        U1Y[1/2],
        SU2L[fund]
    }
];

(* ------------------------------------------------------------------------- *)
(* Yukawa couplings                                                          *)
(* ------------------------------------------------------------------------- *)

AddYukawa[
    yu,
    {H, q, u},
    GroupInvariant -> (
        del[SU3c @ fund, #2, #3] del[SU2L @ fund, #1, #2] &
    ),
    CouplingIndices -> ({gen[#2], gen[#3]} &),
    Chirality -> Right
];

AddYukawa[
    yd,
    {Bar @ H, q, d},
    GroupInvariant -> (
        del[SU3c @ fund, #2, #3] del[SU2L @ fund, #1, #2] &
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

(* ------------------------------------------------------------------------- *)
(* Higgs quartic                                                             *)
(* ------------------------------------------------------------------------- *)

(* This is copied from RGBeta's official SM sample.  We intentionally keep
   RGBeta's normalization here and infer its relation to the project's
   lambda_H convention from the printed beta function. *)
AddQuartic[
    lam,
    {Bar @ H, H, Bar @ H, H},
    GroupInvariant -> (
        del[SU2L @ fund, #1, #2] del[SU2L @ fund, #3, #4] / 2 &
    )
];

(* ------------------------------------------------------------------------- *)
(* Higgs mass parameter                                                      *)
(* ------------------------------------------------------------------------- *)

AddScalarMass[
    mH2,
    {Bar @ H, H},
    GroupInvariant -> (
        del[SU2L @ fund, #1, #2] &
    )
];

(* ------------------------------------------------------------------------- *)
(* Diagnostics                                                               *)
(* ------------------------------------------------------------------------- *)

Print[""];
Print["================ RGBeta SM COUPLING DIAGNOSTIC ================"];

Print[""];
Print["Gauge beta, converted to 16*pi^2 dg/dln(mu):"];
rawGauge = BetaTerm[Gauge, 1] // Simplify;
Print[<|
    gY -> Simplify[rawGauge[gY] / (2 gY)],
    g2 -> Simplify[rawGauge[g2] / (2 g2)],
    g3 -> Simplify[rawGauge[g3] / (2 g3)]
|>];

Print[""];
Print["One-loop Yukawa beta terms:"];
betaY = BetaTerm[Yukawa, 1] // Simplify;
Print[betaY];

Print[""];
Print["One-loop quartic beta terms:"];
betaQ = BetaTerm[Quartic, 1] // Simplify;
Print[betaQ];

Print[""];
Print["One-loop scalar-mass beta terms:"];
betaM2 = BetaTerm[ScalarMass, 1] // Simplify;
Print[betaM2];

Print[""];
Print["Defined coupling types:"];
Print[$couplings];

Print[""];
Print["RGBETA SM COUPLING DIAGNOSTIC: COMPLETE"];

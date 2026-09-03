(* TestRGBetaMaScalar.wl
   RGBeta migration diagnostic for the scalar/gauge sector of the
   Ma (scotogenic / T3-B alpha=-1) model.

   Potential convention:
     V =
       mH2 H^\[Dagger]H + mEta2 eta^\[Dagger]eta
       + (lambda1/2) (H^\[Dagger]H)^2
       + (lambda2/2) (eta^\[Dagger]eta)^2
       + lambda3 (H^\[Dagger]H)(eta^\[Dagger]eta)
       + lambda4 (H^\[Dagger]eta)(eta^\[Dagger]H)
       + (lambda5/2) [(H^\[Dagger]eta)^2 + h.c.]

   This is the same quartic normalization used by RGE/running/MaUVRGE.py.

   We first test the gauge + scalar sector only.  No SM Yukawas and no
   scotogenic Yukawa h are included in this file.  Their contributions are
   added in the next Ma UV test after the scalar convention is verified.
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

(* SM fermions are required for the correct SM gauge running. *)
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

Dim[gen] = 3;

(* ------------------------------------------------------------------------- *)
(* Two scalar doublets                                                       *)
(* ------------------------------------------------------------------------- *)

AddScalar[
    H,
    GaugeRep -> {U1Y[1/2], SU2L[fund]}
];

AddScalar[
    eta,
    GaugeRep -> {U1Y[1/2], SU2L[fund]}
];

(* ------------------------------------------------------------------------- *)
(* Ma / inert-doublet quartics                                               *)
(* ------------------------------------------------------------------------- *)

AddQuartic[
    lambda1,
    {Bar @ H, H, Bar @ H, H},
    GroupInvariant -> (
        del[SU2L @ fund, #1, #2] del[SU2L @ fund, #3, #4] / 2 &
    )
];

AddQuartic[
    lambda2,
    {Bar @ eta, eta, Bar @ eta, eta},
    GroupInvariant -> (
        del[SU2L @ fund, #1, #2] del[SU2L @ fund, #3, #4] / 2 &
    )
];

AddQuartic[
    lambda3,
    {Bar @ H, H, Bar @ eta, eta},
    GroupInvariant -> (
        del[SU2L @ fund, #1, #2] del[SU2L @ fund, #3, #4] &
    )
];

AddQuartic[
    lambda4,
    {Bar @ H, eta, Bar @ eta, H},
    GroupInvariant -> (
        del[SU2L @ fund, #1, #2] del[SU2L @ fund, #3, #4] &
    )
];

AddQuartic[
    lambda5,
    {Bar @ H, eta, Bar @ H, eta},
    GroupInvariant -> (
        del[SU2L @ fund, #1, #2] del[SU2L @ fund, #3, #4] / 2 &
    ),
    SelfConjugate -> False
];

(* ------------------------------------------------------------------------- *)
(* Scalar masses                                                             *)
(* ------------------------------------------------------------------------- *)

AddScalarMass[
    mH2,
    {Bar @ H, H},
    GroupInvariant -> (
        del[SU2L @ fund, #1, #2] &
    )
];

AddScalarMass[
    mEta2,
    {Bar @ eta, eta},
    GroupInvariant -> (
        del[SU2L @ fund, #1, #2] &
    )
];

(* ------------------------------------------------------------------------- *)
(* Diagnostics                                                               *)
(* ------------------------------------------------------------------------- *)

Print[""];
Print["================ RGBeta MA SCALAR DIAGNOSTIC ================"];

rawGauge = BetaTerm[Gauge, 1] // Simplify;

Print[""];
Print["Gauge beta, converted to 16*pi^2 dg/dln(mu):"];
Print[<|
    gY -> Simplify[rawGauge[gY] / (2 gY)],
    g2 -> Simplify[rawGauge[g2] / (2 g2)],
    g3 -> Simplify[rawGauge[g3] / (2 g3)]
|>];

Print[""];
Print["One-loop quartic beta terms (gauge + scalar only):"];
betaQ = BetaTerm[Quartic, 1] // Simplify;
Print[betaQ];

Print[""];
Print["One-loop scalar-mass beta terms (gauge + scalar only):"];
betaM2 = BetaTerm[ScalarMass, 1] // Simplify;
Print[betaM2];

Print[""];
Print["Defined coupling types:"];
Print[$couplings];

Print[""];
Print["RGBETA MA SCALAR DIAGNOSTIC: COMPLETE"];

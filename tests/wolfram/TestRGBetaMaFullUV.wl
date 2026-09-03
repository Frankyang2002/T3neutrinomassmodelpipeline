(* TestRGBetaMaFullUV.wl
   Full RGBeta diagnostic for the Ma / scotogenic / T3-B alpha=-1 UV theory.

   This is the final migration gate before the old generic Python UV-RGE
   generators are removed from the production pipeline.

   Conventions:
     Q = T3 + Y

     V =
       mH2 H^\[Dagger]H + mEta2 eta^\[Dagger]eta
       + (lambda1/2) (H^\[Dagger]H)^2
       + (lambda2/2) (eta^\[Dagger]eta)^2
       + lambda3 (H^\[Dagger]H)(eta^\[Dagger]eta)
       + lambda4 (H^\[Dagger]eta)(eta^\[Dagger]H)
       + (lambda5/2) [(H^\[Dagger]eta)^2 + h.c.]

   N is a gauge-singlet left-handed Weyl field with a symmetric Majorana
   mass matrix M.  h is the scotogenic Yukawa coupling eta L N.

   The file prints all one-loop renormalisable beta functions.  We compare
   the resulting structures with RGE/running/MaUVRGE.py.
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
(* Flavor spaces                                                             *)
(* ------------------------------------------------------------------------- *)

Dim[gen] = 3;
Dim[heavy] = 3;

(* ------------------------------------------------------------------------- *)
(* Standard Model left-handed Weyl fermions                                  *)
(* ------------------------------------------------------------------------- *)

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

(* Ma singlet fermions. *)
AddFermion[
    N,
    FlavorIndices -> {heavy}
];

(* ------------------------------------------------------------------------- *)
(* Scalars                                                                   *)
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
(* Standard Model Yukawa couplings                                           *)
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
(* Scotogenic Yukawa                                                         *)
(* ------------------------------------------------------------------------- *)

(* eta and l are both SU(2) doublets.  Their SU(2)-singlet contraction
   is antisymmetric, 2 x 2 -> 1, so the physical Ma Yukawa uses eps. *)
AddYukawa[
    h,
    {eta, l, N},
    GroupInvariant -> (
        eps[SU2L @ fund, #1, #2] &
    ),
    CouplingIndices -> ({gen[#2], heavy[#3]} &),
    Chirality -> Right
];

(* ------------------------------------------------------------------------- *)
(* Majorana mass                                                             *)
(* ------------------------------------------------------------------------- *)

AddFermionMass[
    M,
    {N, N},
    GroupInvariant -> (1 &),
    MassIndices -> ({heavy[#1], heavy[#2]} &),
    Chirality -> Right
];

(* ------------------------------------------------------------------------- *)
(* Scalar quartics                                                           *)
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
(* One-loop beta functions                                                   *)
(* ------------------------------------------------------------------------- *)

Print[""];
Print["================ RGBeta FULL MA UV DIAGNOSTIC ================"];

rawGauge = BetaTerm[Gauge, 1] // Simplify;
betaGauge = <|
    gY -> Simplify[rawGauge[gY] / (2 gY)],
    g2 -> Simplify[rawGauge[g2] / (2 g2)],
    g3 -> Simplify[rawGauge[g3] / (2 g3)]
|>;

Print[""];
Print["Gauge beta: 16*pi^2 dg/dln(mu)"];
Print[betaGauge];

Print[""];
Print["Yukawa beta terms:"];
betaY = BetaTerm[Yukawa, 1] // Simplify;
Print[betaY];

Print[""];
Print["Fermion-mass beta terms:"];
betaFM = BetaTerm[FermionMass, 1] // Simplify;
Print[betaFM];

Print[""];
Print["Quartic beta terms:"];
betaQ = BetaTerm[Quartic, 1] // Simplify;
Print[betaQ];

Print[""];
Print["Scalar-mass beta terms:"];
betaSM = BetaTerm[ScalarMass, 1] // Simplify;
Print[betaSM];

Print[""];
Print["Defined coupling types:"];
Print[$couplings];

Print[""];
Print["Expected gauge checkpoint: <|gY -> 7 gY^3, g2 -> -3 g2^3, g3 -> -7 g3^3|>"];
Print[""];

gaugePass =
    Simplify[betaGauge[gY] - 7 gY^3] === 0 &&
    Simplify[betaGauge[g2] + 3 g2^3] === 0 &&
    Simplify[betaGauge[g3] + 7 g3^3] === 0;

If[gaugePass,
    Print["RGBETA FULL MA GAUGE CHECK: PASS"],
    Print["RGBETA FULL MA GAUGE CHECK: FAIL"]
];

Print["RGBETA FULL MA UV DIAGNOSTIC: COMPLETE"];

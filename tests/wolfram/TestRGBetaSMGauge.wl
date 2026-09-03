(* TestRGBetaSMGauge.wl
   RGBeta migration regression test.

   Verifies RGBeta against the project's Standard Model gauge convention

       16 pi^2 d g_a / d ln(mu) = b_a g_a^3

   with ordinary hypercharge Q = T3 + Y.

   IMPORTANT:
   RGBeta represents gauge couplings through the gauge-coupling-squared
   structure. At one loop,

       BetaTerm[g, 1] = 2 g * (16 pi^2 d g / d ln(mu)).

   Therefore we divide the RGBeta gauge term by 2 g before comparing with
   the usual beta function for g itself.
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

(* u, d and e are the left-handed conjugates u^c, d^c and e^c. *)
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
(* One-loop gauge beta functions                                             *)
(* ------------------------------------------------------------------------- *)

rawGauge = BetaTerm[Gauge, 1] // Simplify;

(* Convert RGBeta's gauge-coupling-squared beta to beta_g itself. *)
betaGauge = <|
    gY -> Simplify[rawGauge[gY] / (2 gY)],
    g2 -> Simplify[rawGauge[g2] / (2 g2)],
    g3 -> Simplify[rawGauge[g3] / (2 g3)]
|>;

expected = <|
    gY -> (41/6) gY^3,
    g2 -> (-19/6) g2^3,
    g3 -> -7 g3^3
|>;

checks = AssociationMap[
    Simplify[betaGauge[#] - expected[#]] === 0 &,
    Keys[expected]
];

Print["RGBeta raw one-loop gauge terms:"];
Print[rawGauge];

Print[""];
Print["Converted to 16*pi^2 dg/dln(mu):"];
Print[betaGauge];

Print[""];
Print["Expected:"];
Print[expected];

Print[""];
Print["Checks:"];
Print[checks];

If[
    And @@ Values[checks],
    Print["RGBETA SM GAUGE CHECK: PASS"];
    Exit[0],
    Print["RGBETA SM GAUGE CHECK: FAIL"];
    Exit[1]
];

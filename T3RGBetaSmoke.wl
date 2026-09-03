(* T3RGBetaSmoke.wl
   First migration test: verify that RGBeta reproduces the one-loop
   Standard Model gauge beta functions in the project's hypercharge convention.

   Project convention:
       Q = T3 + Y
       16 pi^2 d g / d ln(mu) = beta^(1)_g

   RGBeta's BetaTerm[coupling, 1] returns the one-loop numerator, i.e. the
   quantity to compare directly with the convention above.
*)

ClearAll["Global`*"];

(* Load RGBeta from the Mathematica Applications directory. *)
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

(* u, d and e denote the left-handed conjugates u^c, d^c and e^c. *)
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
        U1Y[1],
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

betaGauge = BetaTerm[Gauge, 1] // Simplify;

Print["RGBeta one-loop SM gauge beta numerators:"];
Print[betaGauge];

expected = <|
    gY -> (41/6) gY^3,
    g2 -> (-19/6) g2^3,
    g3 -> -7 g3^3
|>;

checks = AssociationMap[
    Simplify[betaGauge[#] - expected[#]] === 0 &,
    Keys[expected]
];

Print[""];
Print["Expected convention: 16*pi^2 dg/dln(mu) = beta^(1)_g"];
Print["Expected: ", expected];
Print["Checks:   ", checks];

If[And @@ Values[checks],
    Print["RGBETA SM GAUGE CHECK: PASS"];
    Exit[0],
    Print["RGBETA SM GAUGE CHECK: FAIL"];
    Exit[1]
];

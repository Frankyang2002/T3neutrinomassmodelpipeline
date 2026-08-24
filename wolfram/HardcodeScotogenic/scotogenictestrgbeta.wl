<< RGBeta`
ResetModel[];

(* Define our stuff*)

AddGaugeGroup[g1, U1Y, U1]
AddGaugeGroup[g2, SU2L, SU[2]]
AddGaugeGroup[g3, SU3c, SU[3]]

AddFermion[q, GaugeRep-> {U1Y[1/6], SU2L[fund], SU3c[fund]}, FlavorIndices-> {gen}]
AddFermion[u, GaugeRep-> {U1Y[-2/3], Bar@ SU3c[fund]}, FlavorIndices-> {gen}]
AddFermion[d, GaugeRep-> {U1Y[1/3], Bar@ SU3c[fund]}, FlavorIndices-> {gen}]
AddFermion[l, GaugeRep-> {U1Y[-1/2], SU2L[fund]}, FlavorIndices-> {gen}]
AddFermion[e, GaugeRep-> {U1Y[1]}, FlavorIndices-> {gen}]
(*Neutrino Neutral*)
AddFermion[N, GaugeRep -> {}, FlavorIndices -> {nN}]
Dim[gen] = 3; 
Dim[nN] = 3;


AddScalar[H, GaugeRep-> {U1Y[1/2], SU2L[fund]}]

(*Scalar Neutral*)
AddScalar[eta, GaugeRep -> {U1Y[1/2], SU2L[fund]}]

AddYukawa[
  yu, {H, q, u}, 
  GroupInvariant -> (del[SU3c@ fund, #2, #3] eps[SU2L@ fund, #1, #2] &), 
  CouplingIndices -> ({gen[#2], gen[#3]} &), 
  Chirality -> Right
]

AddYukawa[
  yd, {Bar@ H, q, d},
  GroupInvariant -> (del[SU3c@ fund, #2, #3] del[SU2L@ fund, #1, #2] &),
  CouplingIndices -> ({gen[#2], gen[#3]} &),
  Chirality -> Right
]

AddYukawa[
  ye, {Bar@ H, l, e},
  GroupInvariant -> (del[SU2L@ fund, #1, #2] &),
  CouplingIndices -> ({gen[#2], gen[#3]} &),
  Chirality -> Right
]

(*from the ylnN term*)
AddYukawa[yn, {eta, l, N},
  GroupInvariant -> (
    eps[SU2L @ fund, #2, #1] &
  ),
  CouplingIndices -> ({gen[#2], nN[#3]} &),
  Chirality -> Right
]

AddFermionMass[MN, {N, N},
  GroupInvariant -> (1 &),
  MassIndices -> ({nN[#1], nN[#2]} &),
  Chirality -> Right
]


(* Note to self fermions do not have mass yet due to no EWSB yet *)

AddScalarMass[mH2, {Bar@ H, H},
  GroupInvariant -> (del[SU2L@ fund, #1, #2] &)]

AddScalarMass[meta2, { Bar @ eta, eta},
  GroupInvariant -> (
    del[SU2L @ fund, #1, #2] &
  )
]


(* Potential terms *)

AddQuartic[lambda1, {H, Bar @ H, H, Bar @ H},
  GroupInvariant -> (
    del[SU2L @ fund, #1, #2] del[SU2L @ fund, #3, #4] &
  )
]

AddQuartic[lambda2, {eta, Bar @ eta, eta, Bar @ eta},
  GroupInvariant -> (
    del[SU2L @ fund, #1, #2] del[SU2L @ fund, #3, #4] &
  )
]

AddQuartic[lambda3, {H, Bar @ H, eta, Bar @ eta},
  GroupInvariant -> (
    del[SU2L @ fund, #1, #2] del[SU2L @ fund, #3, #4] &
  )
]

AddQuartic[lambda4, {H, Bar @ eta, eta, Bar @ H},
  GroupInvariant -> (
    del[SU2L @ fund, #1, #2] del[SU2L @ fund, #3, #4] &
  )
]

AddQuartic[lambda5, {Bar @ H, eta, Bar @ H, eta},
  GroupInvariant -> (
    del[SU2L @ fund, #1, #2] del[SU2L @ fund, #3, #4] &
  )
]



betag1 = Finalize @ BetaFunction[g1, 1];
betag2 = Finalize @ BetaFunction[g2, 1];
betag3 = Finalize @ BetaFunction[g3, 1];

betayu = Finalize @ BetaFunction[yu, 1];
betayd = Finalize @ BetaFunction[yd, 1];
betaye = Finalize @ BetaFunction[ye, 1];
betayn = Finalize @ BetaFunction[yn, 1];


betalambda1 = Finalize @ BetaFunction[lambda1, 1];
betalambda2 = Finalize @ BetaFunction[lambda2, 1];
betalambda3 = Finalize @ BetaFunction[lambda3, 1];
betalambda4 = Finalize @ BetaFunction[lambda4, 1];
betalambda5 = Finalize @ BetaFunction[lambda5, 1];

betamH2 = Finalize @ BetaFunction[mH2, 1];
betameta2 = Finalize @ BetaFunction[meta2, 1];
betaMN = Finalize @ BetaFunction[MN, 1];


Print[betag1];

Print["--- Gauge beta functions ---"];
Print["beta(g1) = "];
Print[betag1];

Print["beta(g2) = "];
Print[betag2];

Print["beta(g3) = "];
Print[betag3];

Print["--- SM Yukawa beta functions ---"];
Print["beta(yu) = "];
Print[betayu];

Print["beta(yd) = "];
Print[betayd];

Print[""];
Print["beta(ye) = "];
Print[betaye];

Print["--- Scotogenic Yukawa beta function ---"];
Print["beta(yn) = "];
Print[betayn];

Print["--- Quartic beta functions ---"];
Print["beta(lambda1) = "];
Print[betalambda1];

Print["beta(lambda2) = "];
Print[betalambda2];


Print["beta(lambda3) = "];
Print[betalambda3];

Print["beta(lambda4) = "];
Print[betalambda4];


Print["beta(lambda5) = "];
Print[betalambda5];

Print["--- Scalar mass beta functions ---"];
Print["beta(mH2) = "];
Print[betamH2];

Print["beta(meta2) = "];
Print[betameta2];

Print["--- Majorana mass beta function ---"];
Print["beta(MN) = "];
Print[betaMN];

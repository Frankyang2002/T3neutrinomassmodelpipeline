(* ::Package:: *)

Needs["Matchete`"];

ResetAll[];

Print["Load SM"];

LSM = LoadModel["SM"];

(*--Defining Couplings--*)
Print["Define Couplings"];
DefineCoupling[
  yN,
  Indices -> {Flavor},
  SelfConjugate -> False
];

DefineCoupling[
  lambda3,
  SelfConjugate -> True
];

DefineCoupling[
  lambda4,
  SelfConjugate -> True
];

DefineCoupling[
  lambda5,
  SelfConjugate -> False
];

DefineCoupling[
  lambdaEta,
  SelfConjugate -> True
];


(*--Defining Fields--*)
Print["Defining the heavy fields..."];
DefineField[
  NR,
  Fermion,
  SelfConjugate -> True,
  Mass -> {Heavy, MNR}
];

DefineField[
  Eta,
  Scalar,
  SelfConjugate -> False,
  Indices -> {SU2L[fund]},
  Charges -> {U1Y[-1/2]},
  Mass -> {Heavy, MEta}
];

Print["Heavy fields defined."];


(*--Constructing Scotogenic Lagrangian--*)
Print["Constructing the scotogenic Lagrangian"];

ClearAll[LBSM, LYukawa, V3, V4, V5, VEta];

(* Beyond Standard Model Lagrangian *)
(* p is lepton flavour, ijkm is SU(2) indices *)
LBSM = Module[
  {i, j, k, m, p},
  (* PR is right projection operator defined in Matchete *)
  LYukawa =
    PlusHc @ (
      yN[p]
      Bar @ l[i, p] ** PR ** NR[]
      Eta[i]
    );

  V3 =
    lambda3[]
    Bar @ H[i] H[i]
    Bar @ Eta[j] Eta[j];


  V4 =
    lambda4[]
    H[i] Eta[j]
    Bar @ H[k] Bar @ Eta[m]
    Bar @ CG[eps[SU2L], {i, j}]
    CG[eps[SU2L], {k, m}];


  V5 =
    PlusHc @ (
      lambda5[]/2
      H[i] Eta[j]
      H[k] Eta[m]
      Bar @ CG[eps[SU2L], {i, j}]
      Bar @ CG[eps[SU2L], {k, m}]
    );

  VEta =
    lambdaEta[]/2
    Bar @ Eta[i] Eta[i]
    Bar @ Eta[j] Eta[j];

  (* Free Lagrangian gives mass and kinetic term *)
  (
    FreeLag[NR, Eta]
    + LYukawa
    - V3
    - V4
    - V5
    - VEta
  ) // RelabelIndices
];

LScotogenic = (LSM + LBSM) // RelabelIndices;


(*--Checks on Lagrangian Validity--*)
Print["Checking the BSM Lagrangian..."];

If[
  CheckLagrangian[LBSM] === $Aborted,
  Print["BSM Lagrangian check failed."];
  Abort[],
  Print["BSM Lagrangian check passed."]
];

Print["Checking the complete UV Lagrangian..."];

If[
  CheckLagrangian[LScotogenic] === $Aborted,
  Print["Complete UV Lagrangian check failed."];
  Abort[],
  Print["Complete UV Lagrangian check passed."]
];

LBSM // NiceForm

(*--Integrate out and matching--*)
Print[
  "Integrating out NR and Eta through dimension five..."
];

(*
  LoopOrder -> 1 performs matching through one loop.
  EFTOrder -> 5 retains terms through EFT dimension five.
*)
LEFTraw = Match[
  LScotogenic,
  EFTOrder -> 5,
  LoopOrder -> 1
];

Print["Matching completed."];

LEFTgreen = GreensSimplify[LEFTraw]; (* Gets rid of redundancies *)

Print["Green-basis simplification completed."];

LEFToutput = EOMSimplify[LEFTgreen]; (* More redundancies removed *)

Print["EOM simplification completed."];

LEFTloop = EvaluateLoopFunctions[LEFToutput];

LEFTrep = ReplaceEffectiveCouplings[LEFTrep];


(* Save the EFT Lagrangian *)
outputDirectory =
  FileNameJoin[
    {
      Directory[],
      "scotogenic_output"
    }
  ];

If[
  !DirectoryQ[outputDirectory],
  CreateDirectory[outputDirectory]
];

Put[
  LEFTrep,
  FileNameJoin[
    {
      outputDirectory,
      "scotogenic_eft_lagrangian.wl"
    }
  ]
];

Print["Saved matched EFT Lagrangian to: ",
  FileNameJoin[
    {
      outputDirectory,
      "scotogenic_eft_lagrangian.wl"
    }
  ]
];


 

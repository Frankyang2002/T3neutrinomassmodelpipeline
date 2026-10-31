(* Ordinary T3 Yukawa and topology vertices, plus symmetry-allowed scalar
   operator families. The established two topology Yukawas remain unchanged.

   Workstream A extension (2026-10-31): with
     T3_EXPERIMENTAL_GENERAL_FERMIONS=1
   also enumerate independent Yukawa couplings to the conjugate of a
   hypercharge-zero complex BSM scalar. These additional interactions are
   absent from the original restricted topology basis. The opt-in avoids
   silently altering RGE, matching and numerical workstreams.

   T3_EXPERIMENTAL_GENERAL_QUADRATICS=1 requests BOTH extra quartics and
   extra quadratic operators; the general-scalar beta functions must be
   exported separately with RunGeneralT3QuadraticRGBeta.wl.
   A non-diagonal heavy mass is NOT a valid Matchete matching basis. *)

Get[FileNameJoin[{DirectoryName @ ExpandFileName[$InputFileName],
  "GeneralScalarOperators.wl"}]];

T3ExperimentalGeneralQuadraticsQ[] :=
  Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] === "1";

T3ExperimentalGeneralFermionsQ[] :=
  Environment["T3_EXPERIMENTAL_GENERAL_FERMIONS"] === "1";

(* Only one neutral scalar can occur because Y(S2)-Y(S1)=1. *)
T3GeneralExtraYukawaCount[model_Association] :=
  Boole[TrueQ[model["Scalar1", "Y"] === 0]] +
  Boole[TrueQ[model["Scalar2", "Y"] === 0]];

BuildT3YukawaCandidates[model_Association, which_Integer] := Module[
  {s, dS, dF, y, scalarBar, useCConj, baseName, p, r, i, j, k, f,
   scalar, labels, cg},
  s = model[If[which === 1, "Scalar1", "Scalar2"]];
  dS = s["SU2"];
  dF = model["Fermion", "SU2"];
  useCConj = which === 1;
  scalarBar = which === 2;
  y = If[which === 1, y1, y2];
  baseName = "Yukawa" <> ToString[which];
  cg = If[which === 1, T3Y1CG, T3Y2CG];
  f = If[dF === 1,
    If[useCConj, CConj[NewFermion[r]], NewFermion[r]],
    If[useCConj, CConj[NewFermion[j, r]], NewFermion[j, r]]
  ];
  scalar = If[dS === 1,
    ScalarFieldValue[which, scalarBar],
    ScalarFieldValue[which, scalarBar, k]
  ];
  labels = Join[{i}, If[dF === 1, {}, {j}], If[dS === 1, {}, {k}]];
  {
    <|
      "Name" -> baseName <> "_GenericCG",
      "Class" -> baseName,
      "AlternativeGroup" -> baseName,
      "Expression" -> PlusHc[
        y[p, r] NCM[Bar[l[i, p]], f] scalar CG[cg, labels]
      ]
    |>
  }
];

(* The new Yukawas are not HC duplicates: each uses the same fermion
   chirality as its topology partner, but the OPPOSITE complex scalar.
   SU(2) supplies a unique singlet for each 2 x dF x dS product.
   These candidates still require a live Matchete CheckLagrangian test. *)
BuildT3ExtraYukawaCandidates[model_Association] := Module[
  {dF, dS, p, r, i, j, k, f, scalar, labels, cg, coupling, which,
   conjugated, name},
  If[!TrueQ[T3ExperimentalGeneralFermionsQ[]], Return[{}]];
  If[T3GeneralExtraYukawaCount[model] === 0, Return[{}]];

  dF = model["Fermion", "SU2"];
  which = If[TrueQ[model["Scalar1", "Y"] === 0], 1, 2];
  dS = model[If[which === 1, "Scalar1", "Scalar2"], "SU2"];
  conjugated = which === 1;
  cg = If[which === 1, T3Y1ConjCG, T3Y2ConjCG];
  coupling = If[which === 1, y1Conj, y2Conj];
  name = If[which === 1, "Yukawa1ConjugateScalar", "Yukawa2ConjugateScalar"];

  (* Orient the CG to the actual field objects: Bar[L], F^c/F, S^dag/S. *)
  If[DefineSU2InvariantCG[
      cg, {2, dF, dS},
      If[which === 1, {True, True, True}, {True, False, False}],
      {}, {1}
    ] === $Failed,
    Print["ERROR: extra Yukawa SU(2) tensor could not be defined for ", name];
    Return[$Failed]
  ];
  If[Check[
      Apply[DefineCoupling, {coupling, Indices -> {Flavor, NFlavor},
        SelfConjugate -> False}]; True,
      $Failed
    ] === $Failed,
    Print["ERROR: extra Yukawa coupling could not be registered for ", name];
    Return[$Failed]
  ];

  f = If[dF === 1,
    If[which === 1, CConj[NewFermion[r]], NewFermion[r]],
    If[which === 1, CConj[NewFermion[j, r]], NewFermion[j, r]]
  ];
  scalar = If[dS === 1,
    ScalarFieldValue[which, conjugated],
    ScalarFieldValue[which, conjugated, k]
  ];
  labels = Join[{i}, If[dF === 1, {}, {j}], If[dS === 1, {}, {k}]];

  {<|
    "Name" -> name <> "_GenericCG",
    "Class" -> "ExtraYukawa",
    "AlternativeGroup" -> name,
    "Expression" -> PlusHc[
      coupling[p, r] NCM[Bar[l[i, p]], f] scalar CG[cg, labels]
    ]
  |>}
];

BuildT3MixingCandidates[model_Association] := Module[
  {d1, d2, i, j, a, b, s1, s2, labels},
  d1 = model["Scalar1", "SU2"];
  d2 = model["Scalar2", "SU2"];
  s1 = If[d1 === 1, NewScalar1[], NewScalar1[a]];
  s2 = If[d2 === 1, Bar[NewScalar2[]], Bar[NewScalar2[b]]];
  labels = Join[{i, j}, If[d1 === 1, {}, {a}], If[d2 === 1, {}, {b}]];
  {
    <|
      "Name" -> "T3ScalarMix_GenericCG",
      "Class" -> "T3ScalarMix",
      "AlternativeGroup" -> "T3ScalarMix",
      "Expression" -> PlusHc[
        lambdaT3[] H[i] H[j] s1 s2 CG[T3MixCG, labels]
      ]
    |>
  }
];

BuildT3InteractionCandidates[model_Association] := Module[
  {withMass, withQuartics, extraQ = {}, extraM = {}, extraY = {},
   structures, expectedMass},
  withMass = TrueQ[T3ExperimentalGeneralQuadraticsQ[]];
  withQuartics = withMass || TrueQ[T3ExperimentalGeneralQuarticsQ[]];
  extraY = BuildT3ExtraYukawaCandidates[model];
  If[!ListQ[extraY] ||
     (TrueQ[T3ExperimentalGeneralFermionsQ[]] &&
        Length[extraY] =!= T3GeneralExtraYukawaCount[model]),
    Print["ERROR: independent extra Yukawa construction failed."];
    Return[$Failed]
  ];
  If[withQuartics,
    extraQ = T3ExtendedBuildExtraQuarticCandidates[model];
    If[!ListQ[extraQ],
      Print["ERROR: general-scalar quartic construction failed."];
      Return[$Failed]
    ]
  ];
  If[withMass,
    structures = T3ExtendedExtraQuadraticStructures[model];
    expectedMass = Length[structures];
    extraM = T3ExtendedBuildCandidates[model, 2, structures];
    If[!ListQ[extraM] || Length[extraM] =!= expectedMass ||
       !AllTrue[extraM, Lookup[#, "Class", ""] === "ExtraQuadratic" &],
      Print["ERROR: expected quadratic operator was not constructed."];
      Return[$Failed]
    ];
    Print["GENERAL SCALAR: extra quartics=", Length[extraQ],
      ", extra bilinears=", Length[extraM],
      ". Matchete mass basis requires separate validation."]
  ];
  Join[
    BuildT3YukawaCandidates[model, 1],
    BuildT3YukawaCandidates[model, 2],
    extraY,
    BuildGeneralScalarPotentialCandidates[model],
    BuildT3MixingCandidates[model],
    extraQ,
    extraM
  ]
];

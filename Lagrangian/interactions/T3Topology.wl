(* Ordinary T3 Yukawa and topology vertices, plus the symmetry-allowed
   scalar operator families.  The original restricted model remains the
   default.  T3_EXPERIMENTAL_GENERAL_QUADRATICS=1 requests BOTH extra
   quartics and extra quadratic operators; the general-scalar beta functions
   must be exported separately with RunGeneralT3QuadraticRGBeta.wl.
   A non-diagonal heavy mass term is not automatically a valid Matchete
   matching basis.  BuildT3Lagrangian reports this obstruction explicitly. *)

Get[FileNameJoin[{DirectoryName @ ExpandFileName[$InputFileName],
  "GeneralScalarOperators.wl"}]];

T3ExperimentalGeneralQuadraticsQ[] :=
  Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] === "1";

BuildT3YukawaCandidates[model_Association, which_Integer] := Module[
  {s, dS, dF, y, scalarBar, useCConj, baseName, p, r, i, j, k, f, scalar, labels, cg},
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
  {withMass, withQuartics, extraQ = {}, extraM = {}, structures,
   expectedMass},
  withMass = TrueQ[T3ExperimentalGeneralQuadraticsQ[]];
  withQuartics = withMass || TrueQ[T3ExperimentalGeneralQuarticsQ[]];
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
    BuildGeneralScalarPotentialCandidates[model],
    BuildT3MixingCandidates[model],
    extraQ,
    extraM
  ]
];

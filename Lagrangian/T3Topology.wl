(* T3 topology interaction construction.
   Builds the two lepton-fermion-scalar Yukawa vertices, the HH S1 S2^dagger
   mixing quartic, and combines them with the general scalar potential. *)

(* Build the topology Yukawa with the generic Matchete-generated CG. Remember its yukawa 
1. We load whichever scalar is being used for yukawa here
2. Load its dimensions, whether its conjugated or how it is contracted etc and its couplings
3. Convert our fermion and scalar into conjugated forms dependending on which scalar we use
4. Label indices for CG indices contraction, and give it a name*)
BuildT3YukawaCandidates[model_Association, which_Integer] := Module[
  {s, dS, dF, y, scalarBar, useCConj, baseName, p, r, i, j, k, f, scalar, labels, cg},
  (* Eg: 2x1x2 -> {i,k}, and 2x4x3 -> {i,j,k}  due to singlets getting out*)
  s = model[If[which === 1, "Scalar1", "Scalar2"]];
  dS = s["SU2"];
  dF = model["Fermion", "SU2"];
  useCConj = which === 1;
  scalarBar = which === 2;
  y = If[which === 1, y1, y2];
  baseName = "Yukawa" <> ToString[which];
  cg = If[which === 1, T3Y1CG, T3Y2CG];

  (* We decide on the charge conjugation based on which scalar we use, 
  This is designated on the fact that Y(S1)=alpha/2 while Y(S2)=alpha+2/2
  This gives us the correct contraction for our interactions 
  We have Lbar Fc S1 and Lbar F S2dag
  *)

  (* Create new Fermion that becomes charge conjugated
  based on scalar choice for yukawa*)
  f = If[dF === 1,
    If[useCConj, CConj[NewFermion[r]], NewFermion[r]],
    If[useCConj, CConj[NewFermion[j, r]], NewFermion[j, r]]
  ];

  (* Create new scalar that becomes conjugated based on scalar choice for yukawa *)
  scalar = If[dS === 1,
    ScalarFieldValue[which, scalarBar],
    ScalarFieldValue[which, scalarBar, k]
  ];

  labels = Join[
    {i},
    If[dF === 1, {}, {j}],
    If[dS === 1, {}, {k}]
  ];

  (* Plus HC adds the hermitian conjugate in, NCM is product of dirac spinors *)
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

(* Assemble all interaction candidates without validating them. *)
BuildT3InteractionCandidates[model_Association] := Join[
  BuildT3YukawaCandidates[model, 1],
  BuildT3YukawaCandidates[model, 2],
  BuildGeneralScalarPotentialCandidates[model],
  BuildT3MixingCandidates[model]
];

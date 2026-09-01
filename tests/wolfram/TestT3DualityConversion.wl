
load = UsingFrontEnd[Needs["Matchete`"]; True];
If[load =!= True, Print["Matchete load failed"]; Exit[1]];

ResetAll[];
LSM = LoadModel["SM"];

DefineRepresentation[T3Probe3, SU2L, {2}, IndexAlphabet -> {"a","b","c","d"}];
DefineRepresentation[T3Probe4, SU2L, {3}, IndexAlphabet -> {"u","v","w","x"}];

(* Canonical generators *)
gf = Normal @ Generators[SU[2], {1}];
g4 = Normal @ Generators[SU[2], {3}];
g3 = Normal @ Generators[SU[2], {2}];

(* Raw invariant: this was found to transform with G^T on all three axes. *)
Traw = Normal @ First @ InvariantTensors[SU[2], {{1},{3},{2}}];

(* SU(2) self-duality intertwiners.  For an irrep R, J obeys
       G^T J + J G = 0
   and converts a dual/covariant index into the representation convention
   used by Matchete's registered index type. *)
J2 = Normal @ First @ InvariantTensors[SU[2], {{1},{1}}];
J4 = Normal @ First @ InvariantTensors[SU[2], {{3},{3}}];
J3 = Normal @ First @ InvariantTensors[SU[2], {{2},{2}}];

Print["============================================================"];
Print["SELF-DUALITY MATRICES"];
Print["============================================================"];
Print["J2 = ", InputForm[J2]];
Print["J4 = ", InputForm[J4]];
Print["J3 = ", InputForm[J3]];
Print["det(J2), det(J4), det(J3) = ",
  InputForm[{Det[J2], Det[J4], Det[J3]}]];

relationResidual[gens_, J_] := Table[
  FullSimplify[Transpose[gens[[A]]] . J + J . gens[[A]]],
  {A, Length[gens]}
];

zeroMatrixQ[m_] := And @@ Flatten[Map[TrueQ[FullSimplify[# == 0]] &, m, {2}]];

Print["J relations: ",
  {
    And @@ (zeroMatrixQ /@ relationResidual[gf,J2]),
    And @@ (zeroMatrixQ /@ relationResidual[g4,J4]),
    And @@ (zeroMatrixQ /@ relationResidual[g3,J3])
  }
];

applyOnAxis[t_, a_, axis_Integer] := Module[
  {dims = Dimensions[t], out, inds, old},
  out = ConstantArray[0, dims];
  Do[
    inds = {i,j,k};
    Do[
      old = ReplacePart[inds, axis -> q];
      out[[i,j,k]] += a[[inds[[axis]],q]] Extract[t,old],
      {q,dims[[axis]]}
    ],
    {i,dims[[1]]},{j,dims[[2]]},{k,dims[[3]]}
  ];
  FullSimplify[out]
];

(* Convert every raw covariant/dual index with J^{-1}. *)
Tmat = applyOnAxis[Traw, Inverse[J2], 1];
Tmat = applyOnAxis[Tmat, Inverse[J4], 2];
Tmat = applyOnAxis[Tmat, Inverse[J3], 3];
Tmat = FullSimplify[Tmat];

Print[""];
Print["============================================================"];
Print["CONVERTED TENSOR"];
Print["============================================================"];
Print["Tmat = ", InputForm[Tmat]];

zeroTensorQ[x_] := And @@ Flatten[Map[TrueQ[FullSimplify[# == 0]] &, x, {3}]];

residuals = Table[
  FullSimplify[
    applyOnAxis[Tmat, gf[[A]], 1] +
    applyOnAxis[Tmat, g4[[A]], 2] +
    applyOnAxis[Tmat, g3[[A]], 3]
  ],
  {A,3}
];

Print["Natural G-action invariant after conversion? ",
  And @@ (zeroTensorQ /@ residuals)];

Do[
  nz = Count[Flatten[residuals[[A]]], x_ /; !TrueQ[FullSimplify[x == 0]]];
  Print["generator ", A, ": nonzero residual components = ", nz],
  {A,3}
];

Print[""];
Print["============================================================"];
Print["DEFINECG TEST"];
Print["============================================================"];

(* Only attempt DefineCG if the direct generator test succeeded. *)
If[And @@ (zeroTensorQ /@ residuals),
  status = CheckAbort[
    Check[
      UsingFrontEnd[
        DefineCG[
          T3ProbeCG243,
          {SU2L[fund], T3Probe4, T3Probe3},
          Tmat
        ];
        True
      ],
      False
    ],
    False
  ];
  Print["DefineCG accepted converted tensor? ", status],
  Print["Skipping DefineCG: converted tensor did not pass direct invariance."]
];


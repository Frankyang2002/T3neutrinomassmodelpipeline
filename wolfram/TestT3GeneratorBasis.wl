
load = UsingFrontEnd[Needs["Matchete`"]; True];
If[load =!= True, Print["Matchete load failed"]; Exit[1]];

ResetAll[];
LSM = LoadModel["SM"];

DefineRepresentation[T3Probe3, SU2L, {2}, IndexAlphabet -> {"a","b","c","d"}];
DefineRepresentation[T3Probe4, SU2L, {3}, IndexAlphabet -> {"u","v","w","x"}];

Print["============================================================"];
Print["CANONICAL SU(2) GENERATORS"];
Print["============================================================"];

gf = Normal @ Generators[SU[2], {1}];
g4 = Normal @ Generators[SU[2], {3}];
g3 = Normal @ Generators[SU[2], {2}];

Print["fund dimensions = ", Dimensions[gf]];
Print["quartet dimensions = ", Dimensions[g4]];
Print["triplet dimensions = ", Dimensions[g3]];

Print[""];
Print["============================================================"];
Print["RAW 2 x 4 x 3 INVARIANT"];
Print["============================================================"];

tensors = CheckAbort[
  Check[InvariantTensors[SU[2], {{1},{3},{2}}], $Failed],
  $Failed
];

If[tensors === $Failed || !ListQ[tensors] || Length[tensors] == 0,
  Print["InvariantTensors failed."];
  Exit[2];
];

T = Normal @ First[tensors];

Print["Tensor dimensions = ", Dimensions[T]];
Print["Tensor = ", InputForm[T]];

applyOnAxis[t_, g_, axis_Integer] := Module[
  {dims = Dimensions[t], out, inds, old},
  out = ConstantArray[0, dims];

  Do[
    inds = {i,j,k};

    Do[
      old = ReplacePart[inds, axis -> a];
      out[[i,j,k]] += g[[inds[[axis]], a]] * Extract[t, old],
      {a, dims[[axis]]}
    ],
    {i, dims[[1]]},
    {j, dims[[2]]},
    {k, dims[[3]]}
  ];

  FullSimplify[out]
];

zeroTensorQ[x_] := And @@ Flatten[Map[TrueQ[FullSimplify[# == 0]] &, x, {3}]];

variants[g_] := <|
  "G" -> g,
  "GT" -> Transpose[g],
  "-G" -> -g,
  "-GT" -> -Transpose[g],
  "Gc" -> Conjugate[g],
  "-Gc" -> -Conjugate[g],
  "Gdag" -> ConjugateTranspose[g],
  "-Gdag" -> -ConjugateTranspose[g]
|>;

Print[""];
Print["============================================================"];
Print["INFINITESIMAL INVARIANCE TEST"];
Print["============================================================"];

vf = variants /@ gf;
v4 = variants /@ g4;
v3 = variants /@ g3;

names = Keys[First[vf]];
passes = {};

Do[
  residuals = Table[
    FullSimplify[
      applyOnAxis[T, vf[[A]][n1], 1] +
      applyOnAxis[T, v4[[A]][n2], 2] +
      applyOnAxis[T, v3[[A]][n3], 3]
    ],
    {A, 3}
  ];

  If[And @@ (zeroTensorQ /@ residuals),
    AppendTo[passes, {n1,n2,n3}];
    Print["PASS ", InputForm[{n1,n2,n3}]]
  ],
  {n1,names},
  {n2,names},
  {n3,names}
];

Print[""];
Print["============================================================"];
Print["RESULT"];
Print["============================================================"];

If[passes === {},
  Print["NO PASSING GENERATOR CONVENTION."];
  Print["Interpretation: the raw InvariantTensors tensor is in a basis that does not directly match the canonical generator action tested here."],
  Print["Passing convention count = ", Length[passes]];
  Scan[Print[InputForm[#]] &, passes]
];

(* Also directly test the mathematically natural convention first, and print
   the maximum number of nonzero residual components for each generator. *)
Print[""];
Print["============================================================"];
Print["NATURAL-CONVENTION RESIDUAL"];
Print["============================================================"];

natural = Table[
  FullSimplify[
    applyOnAxis[T, gf[[A]], 1] +
    applyOnAxis[T, g4[[A]], 2] +
    applyOnAxis[T, g3[[A]], 3]
  ],
  {A,3}
];

Do[
  nz = Count[Flatten[natural[[A]]], x_ /; !TrueQ[FullSimplify[x == 0]]];
  Print["generator ", A, ": nonzero residual components = ", nz],
  {A,3}
];

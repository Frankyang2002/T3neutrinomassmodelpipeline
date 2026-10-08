(* Export the SU(2) generators used by Matchete's Lie-algebra interface.
   This does NOT establish that the matrices use the same index basis as
   registered T3Y1CG/T3Y2CG/T3MixCG tensors. That requires a covariance test.

   Usage:
   wolframscript -file Lagrangian/interactions/ProbeT3Generators.wl output.json
*)
ClearAll["Global`*"];

If[Length[$ScriptCommandLine] < 2,
  Print["Usage: wolframscript -file ProbeT3Generators.wl output.json"];
  Exit[2]
];

outputPath = ExpandFileName[$ScriptCommandLine[[2]]];
If[!DirectoryQ[DirectoryName[outputPath]],
  CreateDirectory[DirectoryName[outputPath], CreateIntermediateDirectories -> True]
];

Print["Loading Matchete..."];
loaded = Quiet @ Check[UsingFrontEnd[Needs["Matchete`"]; True], $Failed];
If[loaded === $Failed || !TrueQ[loaded],
  Print["ERROR: Matchete load failed"];
  Exit[3]
];

ClearAll[exportGenerators];
exportGenerators[d_Integer?Positive] := Module[
  {generators, matrices, hermiticity, lieResiduals, casimir, j, checks},
  generators = Quiet @ Check[Generators[SU[2], {d - 1}], $Failed];
  If[generators === $Failed || Dimensions[generators] =!= {3, d, d},
    Return[<|"status" -> "ERROR", "dimension" -> d,
      "message" -> "Generators returned incorrect dimensions or failed"|>]
  ];
  matrices = Table[Normal[generators[[a]]], {a, 1, 3}];
  hermiticity = Table[
    TrueQ[Simplify[matrices[[a]] == ConjugateTranspose[matrices[[a]]]]],
    {a, 1, 3}
  ];
  lieResiduals = {
    Simplify[matrices[[1]].matrices[[2]] - matrices[[2]].matrices[[1]] - I matrices[[3]]],
    Simplify[matrices[[2]].matrices[[3]] - matrices[[3]].matrices[[2]] - I matrices[[1]]],
    Simplify[matrices[[3]].matrices[[1]] - matrices[[1]].matrices[[3]] - I matrices[[2]]]
  };
  j = (d - 1)/2;
  casimir = Simplify[Total[Map[#.# &, matrices]] - j (j + 1) IdentityMatrix[d]];
  checks = <|
    "Hermitian" -> And @@ hermiticity,
    "su2_commutators" -> And @@ (TrueQ[# == ConstantArray[0, {d, d}]] & /@ lieResiduals),
    "Casimir" -> TrueQ[casimir == ConstantArray[0, {d, d}]]
  |>;
  <|
    "status" -> If[And @@ Values[checks], "Success", "CheckFailed"],
    "dimension" -> d,
    "dynkin_label" -> {d - 1},
    "generator_order" -> {"T1", "T2", "T3"},
    "generators_input_form" -> Table[
      Table[ToString[InputForm[matrices[[a, r, c]]]], {r, d}, {c, d}],
      {a, 3}
    ],
    "checks" -> checks,
    "basis_relation_to_exported_CG" -> "unverified"
  |>
];

results = Table[exportGenerators[d], {d, {2, 3}}];
payload = <|
  "status" -> If[AllTrue[results, Lookup[#, "status", ""] === "Success" &], "Success", "CheckFailed"],
  "method" -> "Generators[SU[2], {d-1}] (Matchete Lie-algebra interface)",
  "basis_relation_to_T3_CG" -> "unverified; covariance check required",
  "representations" -> results
|>;

Export[outputPath, payload, "RawJSON"];
Print["Generator export status: ", payload["status"]];
Print["Generator JSON: ", outputPath];
If[payload["status"] =!= "Success", Exit[6]];

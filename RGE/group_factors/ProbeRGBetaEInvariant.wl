(* Focused probe for the T3-E (dS1=dS2=3) RGBeta invariant.

Run from repository root:
  wolframscript -file RGE/group_factors/ProbeRGBetaEInvariant.wl

This prints the exact RGBeta-refined forms of:
  1. the triplet triple invariant,
  2. the T3-E HH S1 Bar[S2] mixing invariant,
  3. the triplet portal adjoint invariant,
  4. the S1-S2 adjoint invariant,
  5. the exact T3-E Cross invariant.

The output is intentionally restricted to these tensors only.
*)

ClearAll["Global`*"];

Needs["RGBeta`"];

modelPath = FileNameJoin[{
    DirectoryName[$InputFileName],
    "..", "running", "wolfram", "T3RGBetaModel.wl"
}];

If[!FileExistsQ[modelPath],
    modelPath = FileNameJoin[{
        Directory[],
        "RGE", "running", "wolfram", "T3RGBetaModel.wl"
    }]
];

If[!FileExistsQ[modelPath],
    Print["ERROR: could not locate T3RGBetaModel.wl"];
    Exit[1];
];

Get[modelPath];

Print["=== T3-E RGBeta invariant probe ==="];
Print["Model source: ", modelPath];

triple =
    T3RGBetaS2TripleInvariant[a, b, c];

mix33 =
    delS2[SU2L, a, h1, h2] *
    T3RGBetaS2TripleInvariant[a, s1, s2];

portalAdj =
    tGen[SU2L @ fund, A, hbar, h] *
    tGen[SU2L @ S2, A, sbar, s];

scalarAdj =
    tGen[SU2L @ S2, A, s1bar, s1] *
    tGen[SU2L @ S2, A, s2bar, s2];

crossExact =
    del[SU2L @ S2, s1bar, s2bar] *
    del[SU2L @ S2, s1, s2] +
    del[SU2L @ S2, s1bar, s2] *
    del[SU2L @ S2, s1, s2bar];

Print["\n=== Triple invariant raw ==="];
Print[InputForm[triple]];

Print["\n=== Triple invariant refined ==="];
Print[InputForm[RGBeta`PackageScope`RefineGroupStructures[triple]]];

Print["\n=== T3-E mixing invariant raw ==="];
Print[InputForm[mix33]];

Print["\n=== T3-E mixing invariant refined ==="];
Print[InputForm[RGBeta`PackageScope`RefineGroupStructures[mix33]]];

Print["\n=== H-triplet portal adjoint refined ==="];
Print[InputForm[RGBeta`PackageScope`RefineGroupStructures[portalAdj]]];

Print["\n=== S1-S2 adjoint refined ==="];
Print[InputForm[RGBeta`PackageScope`RefineGroupStructures[scalarAdj]]];

Print["\n=== T3-E Cross raw ==="];
Print[InputForm[crossExact]];

Print["\n=== T3-E Cross refined ==="];
Print[InputForm[RGBeta`PackageScope`RefineGroupStructures[crossExact]]];

Print["\nDONE"];
Exit[0];

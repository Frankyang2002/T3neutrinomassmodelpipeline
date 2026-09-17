(* Probe RGBeta's symmetric-two-index S2 normalization and contractions.

Run from repository root:
  wolframscript -file RGE/group_factors/ProbeRGBetaS2Contractions.wl

This uses RGBeta's own RefineGroupStructures routine.  The objective is to
determine the normalization induced by delS2 and tGen[S2], and to test the
generator action on delS2 before constructing the full lambdaT3 Ward identity.
*)

ClearAll["Global`*"];
Needs["RGBeta`"];

modelFile = FileNameJoin[{
    Directory[],
    "RGE", "running", "wolfram", "T3RGBetaModel.wl"
}];
If[FileExistsQ[modelFile], Get[modelFile]];

refine = RGBeta`PackageScope`RefineGroupStructures;

Print["=== RGBeta S2 contraction / intertwiner probe ==="];

tests = {
    "delS2 pair, shared S2 index" ->
        delS2[SU2L, a, i, j] delS2[SU2L, a, k, l],

    "delS2 pair, explicit S2 delta" ->
        delS2[SU2L, a, i, j] del[SU2L @ S2, a, b] delS2[SU2L, b, k, l],

    "S2 generator sandwiched by delS2" ->
        delS2[SU2L, a, i, j] *
        tGen[SU2L @ S2, A, a, b] *
        delS2[SU2L, b, k, l],

    "fund generator acting on first symmetric leg" ->
        tGen[SU2L @ fund, A, i, p] *
        delS2[SU2L, a, p, j],

    "fund generator acting on second symmetric leg" ->
        tGen[SU2L @ fund, A, j, p] *
        delS2[SU2L, a, i, p],

    "S2 generator acting on triplet leg" ->
        tGen[SU2L @ S2, A, a, b] *
        delS2[SU2L, b, i, j]
};

Do[
    Print["\n--- ", label, " ---"];
    Print["raw:     ", InputForm[expr]];
    result = Quiet @ Check[refine[expr], $Failed];
    Print["refined: ", InputForm[result]],
    {label -> expr, tests}
];

Print["\n=== Candidate intertwiner residuals ==="];

rPlus =
    tGen[SU2L @ S2, A, a, b] delS2[SU2L, b, i, j]
    - tGen[SU2L @ fund, A, i, p] delS2[SU2L, a, p, j]
    - tGen[SU2L @ fund, A, j, p] delS2[SU2L, a, i, p];

rMinus =
    tGen[SU2L @ S2, A, a, b] delS2[SU2L, b, i, j]
    + tGen[SU2L @ fund, A, i, p] delS2[SU2L, a, p, j]
    + tGen[SU2L @ fund, A, j, p] delS2[SU2L, a, i, p];

Print["R_plus raw: ", InputForm[rPlus]];
Print["R_plus refined: ", InputForm @ Quiet @ Check[refine[rPlus], $Failed]];
Print["R_minus raw: ", InputForm[rMinus]];
Print["R_minus refined: ", InputForm @ Quiet @ Check[refine[rMinus], $Failed]];

Print["\n=== T3-B mixing invariant refinement ==="];
mix22 =
    delS2[SU2L, a, h1, h2] *
    delS2[SU2L, a, s1, s2];

Print["I22 raw: ", InputForm[mix22]];
Print["I22 refined: ", InputForm @ Quiet @ Check[refine[mix22], $Failed]];

Print["\nDONE"];
Exit[0];

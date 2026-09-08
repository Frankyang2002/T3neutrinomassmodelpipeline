(* TestT3Conventions.wl
   Consolidated representation/convention regression tests.
   Replaces the former SU2, generic-representation, CG/duality and
   representation-API probe files.
*)

ClearAll["Global`*"];

scriptDirectory = DirectoryName @ ExpandFileName[$InputFileName];
projectRoot = ExpandFileName @ FileNameJoin[{scriptDirectory, "..", ".."}];

Fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);
Assert[label_, condition_] := If[TrueQ[condition], Print["PASS: ", label], Fail[label]];
AssertEqual[label_, actual_, expected_] := Assert[label <> " (got " <> ToString[actual, InputForm] <> ")", actual === expected];

(* ---------------------------------------------------------------------- *)
(* Project SU(2) representation logic                                      *)
(* ---------------------------------------------------------------------- *)

gaugeFile = FileNameJoin[{projectRoot, "wolfram", "core", "GaugeInvariance.wl"}];
catalogFile = FileNameJoin[{projectRoot, "wolfram", "t3", "T3ModelCatalog.wl"}];

If[!FileExistsQ[gaugeFile], Fail["GaugeInvariance.wl not found at " <> gaugeFile]];
If[!FileExistsQ[catalogFile], Fail["T3ModelCatalog.wl not found at " <> catalogFile]];

Get[gaugeFile];
Get[catalogFile];

AssertEqual["2 x 2", SU2Combine[2, 2], {1, 3}];
AssertEqual["2 x 3", SU2Combine[2, 3], {2, 4}];
AssertEqual["3 x 3", SU2Combine[3, 3], {1, 3, 5}];
Assert["2 x 2 contains singlet", SU2SingletQ[{2, 2}]];
Assert["2 x 3 has no singlet", !SU2SingletQ[{2, 3}]];
Assert["2 x 2 x 3 contains singlet", SU2SingletQ[{2, 2, 3}]];

(* Benchmark classes and representative higher-dimensional generalisations. *)
Do[
  Assert["allowed dimensions " <> ToString[dims, InputForm], T3DimensionsAllowedQ @@ dims],
  {dims, {{1,3,2}, {2,2,1}, {2,2,3}, {3,1,2}, {3,3,2}, {3,3,4}, {3,5,4}, {5,3,4}, {5,5,4}, {4,6,5}}}
];
Assert["forbid non-closing (1,1,2)", !T3DimensionsAllowedQ[1,1,2]];
Assert["forbid non-adjacent scalar (2,6,5)", !T3DimensionsAllowedQ[2,6,5]];

(* ---------------------------------------------------------------------- *)
(* Matchete representation / SU(2) invariant convention                    *)
(* ---------------------------------------------------------------------- *)

matcheteLoaded = UsingFrontEnd[Needs["Matchete`"]; True];
If[!TrueQ[matcheteLoaded], Fail["Matchete failed to load"]];

ResetAll[];
LoadModel["SM"];

DefineRepresentation[T3Probe3, SU2L, {2}, IndexAlphabet -> {"a","b","c","d"}];
DefineRepresentation[T3Probe4, SU2L, {3}, IndexAlphabet -> {"u","v","w","x"}];

reps = GetRepresentations[];
Assert["custom triplet representation registered", !FreeQ[reps, T3Probe3, Infinity]];
Assert["custom quartet representation registered", !FreeQ[reps, T3Probe4, Infinity]];

rawInvariant = Quiet @ Check[InvariantTensors[SU[2], {{1}, {3}, {2}}], $Failed];
Assert["2 x 4 x 3 invariant exists", ListQ[rawInvariant] && Length[rawInvariant] >= 1];
AssertEqual["2 x 4 x 3 invariant dimensions", Dimensions[Normal @ First[rawInvariant]], {2,4,3}];

(* SU(2) irreps are self-dual.  The invariant bilinear must exist for each
   representation used in the convention conversion. *)
Do[
  bilinear = Quiet @ Check[InvariantTensors[SU[2], {rep, rep}], $Failed];
  Assert["self-duality invariant " <> ToString[rep, InputForm], ListQ[bilinear] && Length[bilinear] >= 1],
  {rep, {{{1}}, {{2}}, {{3}}}}
];

Print["ALL T3 CONVENTION TESTS PASSED"];
Exit[0];

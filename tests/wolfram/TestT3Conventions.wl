(* TestT3Conventions.wl
   Consolidated representation/convention regression tests.
   Replaces the former SU2, generic-representation, CG/duality and
   representation-API probe files.
*)

ClearAll["Global`*"];

scriptDirectory = DirectoryName @ ExpandFileName[$InputFileName];
projectRoot = ExpandFileName @ FileNameJoin[{scriptDirectory, "..", ".."}];

TestFail[msg_] := (Print["FAIL: ", msg]; Exit[1]);
Assert[label_, condition_] := If[TrueQ[condition], Print["PASS: ", label], TestFail[label]];
AssertEqual[label_, actual_, expected_] := Assert[label <> " (got " <> ToString[actual, InputForm] <> ")", actual === expected];

(* ---------------------------------------------------------------------- *)
(* Project SU(2) representation logic                                      *)
(* ---------------------------------------------------------------------- *)

su2File = FileNameJoin[{projectRoot, "Lagrangian", "model", "SU2Invariants.wl"}];
catalogFile = FileNameJoin[{projectRoot, "Lagrangian", "model", "T3ModelCatalog.wl"}];

If[!FileExistsQ[su2File], TestFail["SU2Invariants.wl not found at " <> su2File]];
If[!FileExistsQ[catalogFile], TestFail["T3ModelCatalog.wl not found at " <> catalogFile]];

Get[catalogFile];

(* Current production SU(2) representation convention: dimension d maps to
   highest-weight Dynkin label {d-1}. *)
AssertEqual["d=2 Dynkin label", {2 - 1}, {1}];
AssertEqual["d=3 Dynkin label", {3 - 1}, {2}];
AssertEqual["d=4 Dynkin label", {4 - 1}, {3}];

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
If[!TrueQ[matcheteLoaded], TestFail["Matchete failed to load"]];

Get[su2File];
AssertEqual["SU2DynkinLabel[2]", SU2DynkinLabel[2], {1}];
AssertEqual["SU2DynkinLabel[3]", SU2DynkinLabel[3], {2}];
AssertEqual["SU2DynkinLabel[4]", SU2DynkinLabel[4], {3}];

representationProbe = UsingFrontEnd[
  ResetAll[];
  LoadModel["SM"];

  tripletRegistered = EnsureBSMRepresentation[3];
  quartetRegistered = EnsureBSMRepresentation[4];

  <|
    "TripletStatus" -> tripletRegistered,
    "QuartetStatus" -> quartetRegistered,
    "Representations" -> GetRepresentations[]
  |>
];

If[!AssociationQ[representationProbe],
  TestFail["could not read Matchete representation registry"]
];

Assert[
  "production triplet representation registered",
  TrueQ[representationProbe["TripletStatus"]]
];
Assert[
  "production quartet representation registered",
  TrueQ[representationProbe["QuartetStatus"]]
];
Assert[
  "triplet registry key present",
  KeyExistsQ[representationProbe["Representations"], T3BSMd3]
];
Assert[
  "quartet registry key present",
  KeyExistsQ[representationProbe["Representations"], T3BSMd4]
];

(* The full production invariant-tensor path is exercised by
   TestT3RGEExport.wl through BuildT3Lagrangian.  Keep this conventions test
   focused on the representation and model conventions themselves; direct
   standalone InvariantTensors probes are front-end fragile under wolframscript. *)

Print["ALL T3 CONVENTION TESTS PASSED"];
Exit[0];

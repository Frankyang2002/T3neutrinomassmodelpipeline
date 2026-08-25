(* TestGenericT3Reps.wl
   Fast representation-level tests.  No matching is performed here. *)
ClearAll[assert];
assert[name_, cond_] := If[TrueQ[cond], Print["PASS: ", name], Print["FAIL: ", name]; Exit[1]];

scriptDirectory = DirectoryName@ExpandFileName[$InputFileName];
Get[FileNameJoin@{scriptDirectory, "T3ModelCatalog.wl"}];

assert["B = (2,2,1)", T3DimensionsAllowedQ[2,2,1]];
assert["A = (1,3,2)", T3DimensionsAllowedQ[1,3,2]];
assert["D = (3,1,2)", T3DimensionsAllowedQ[3,1,2]];
assert["E = (3,3,2)", T3DimensionsAllowedQ[3,3,2]];
assert["C = (2,2,3)", T3DimensionsAllowedQ[2,2,3]];

assert["forbid (1,1,2): HH cannot close", !T3DimensionsAllowedQ[1,1,2]];
assert["quartet F: (3,3,4)", T3DimensionsAllowedQ[3,3,4]];
assert["quartet F: (3,5,4)", T3DimensionsAllowedQ[3,5,4]];
assert["quartet F: (5,3,4)", T3DimensionsAllowedQ[5,3,4]];
assert["quartet F: (5,5,4)", T3DimensionsAllowedQ[5,5,4]];
assert["quintet F: (4,6,5)", T3DimensionsAllowedQ[4,6,5]];
assert["forbid non-adjacent scalar", !T3DimensionsAllowedQ[2,6,5]];

Print["All generic T3 representation tests passed."];
Exit[0];

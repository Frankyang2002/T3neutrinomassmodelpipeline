(*
  TestSU2Invariance.wl
  Representation-layer tests for GaugeInvariance.wl.
  Run from the directory containing GaugeInvariance.wl with:
      wolframscript -file TestSU2Invariance.wl
*)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "core", "GaugeInvariance.wl"}]];

ClearAll[AssertEqual];
AssertEqual[label_String, actual_, expected_] := Module[{},
  If[actual === expected,
    Print["PASS: ", label, " -> ", actual],
    Print["FAIL: ", label, " -> got ", actual, ", expected ", expected];
    Exit[1]
  ]
];

AssertEqual["2 x 2", SU2Combine[2, 2], {1, 3}];
AssertEqual["2 x 3", SU2Combine[2, 3], {2, 4}];
AssertEqual["3 x 3", SU2Combine[3, 3], {1, 3, 5}];

AssertEqual["2 x 2 contains singlet", SU2SingletQ[{2, 2}], True];
AssertEqual["2 x 3 does not contain singlet", SU2SingletQ[{2, 3}], False];
AssertEqual["2 x 2 x 3 contains singlet", SU2SingletQ[{2, 2, 3}], True];

baseTemplate = <|
  "Name" -> "Yukawa",
  "Fields" -> {
    <|"Name" -> "L", "Conjugated" -> True|>,
    <|"Name" -> "Scalar", "Conjugated" -> True|>,
    <|"Name" -> "Fermion", "Conjugated" -> False|>
  }
|>;

fieldDataForScalar[d_Integer] := <|
  "L" -> <|"SU2" -> 2, "Y" -> -1/2|>,
  "Scalar" -> <|"SU2" -> d, "Y" -> 1/2|>,
  "Fermion" -> <|"SU2" -> 1, "Y" -> 0|>
|>;

AssertEqual[
  "scalar singlet -> fermion doublet",
  InferFieldSU2Representations[baseTemplate, fieldDataForScalar[1], "Fermion"],
  {2}
];
AssertEqual[
  "scalar doublet -> fermion singlet or triplet",
  InferFieldSU2Representations[baseTemplate, fieldDataForScalar[2], "Fermion"],
  {1, 3}
];
AssertEqual[
  "scalar triplet -> fermion doublet or quartet",
  InferFieldSU2Representations[baseTemplate, fieldDataForScalar[3], "Fermion"],
  {2, 4}
];

Print["All SU(2) representation-layer tests passed."];


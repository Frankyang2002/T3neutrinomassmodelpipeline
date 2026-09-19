(* Final regression gate for the validated T3 UV -> one-loop matching -> C5 stage.
   Called by regression.py after the smoke scan has produced output directories. *)

ClearAll["Global`*"];
scriptDirectory = DirectoryName @ ExpandFileName[$InputFileName];
args = Rest[$ScriptCommandLine];
outputRoot = If[Length[args] >= 1, ExpandFileName @ args[[1]], FileNameJoin[{scriptDirectory, "output"}]];


projectRoot = ExpandFileName @ FileNameJoin[{scriptDirectory, "..", ".."}];
conventionFile = FileNameJoin[{projectRoot, "tests", "wolfram", "T3CouplingConventions.wl"}];
If[!FileExistsQ[conventionFile], Print["Missing coupling convention file: ", conventionFile]; Exit[2]];
Get[conventionFile];
matcheteLoadResult = UsingFrontEnd[Needs["Matchete`"]; True];
If[!TrueQ[matcheteLoadResult], Print["Matchete loading failed in regression script."]; Exit[2]];

benchmarks = {
  <|"Class" -> "B", "Alpha" -> -1, "Dir" -> "T3_B_alpha_m1"|>,
  <|"Class" -> "C", "Alpha" -> -1, "Dir" -> "T3_C_alpha_m1"|>,
  <|"Class" -> "A", "Alpha" -> 0,  "Dir" -> "T3_A_alpha_p0"|>,
  <|"Class" -> "D", "Alpha" -> -2, "Dir" -> "T3_D_alpha_m2"|>,
  <|"Class" -> "E", "Alpha" -> 0,  "Dir" -> "T3_E_alpha_p0"|>
};

results = {};

AddResult[name_, pass_, detail_: ""] := AppendTo[results, <|
  "Test" -> name,
  "Pass" -> TrueQ[pass],
  "Detail" -> ToString[detail, InputForm]
|>];

SafeZeroQ[x_] := TrueQ[Quiet @ Check[PossibleZeroQ[Together[x]], False]];

ReadSummary[dir_] := Module[{path = FileNameJoin[{outputRoot, dir, "comparison_summary.json"}]},
  If[!FileExistsQ[path], Return[$Failed]];
  Quiet @ Check[Import[path, "RawJSON"], $Failed]
];

ReadC5[dir_] := Module[
  {path = FileNameJoin[
    {outputRoot, dir, "data", "c5_coefficient.txt"}
  ], text},

  If[!FileExistsQ[path], Return[$Failed]];

  text = Import[path, "Text"];
  Quiet @ Check[ToExpression[text, InputForm], $Failed]
];

Do[
  b = benchmark;
  summary = ReadSummary[b["Dir"]];
  prefix = "T3-" <> b["Class"] <> " alpha=" <> ToString[b["Alpha"]];

  AddResult[prefix <> " summary exists", AssociationQ[summary]];
  If[AssociationQ[summary],
    AddResult[prefix <> " build", Lookup[summary, "BuildStatus", ""] === "Success"];
    AddResult[prefix <> " matching", Lookup[summary, "MatchingStatus", ""] === "Success"];
    AddResult[prefix <> " T3 ingredients", TrueQ[Lookup[summary, "T3IngredientsPresent", False]]];
    AddResult[prefix <> " Weinberg present", TrueQ[Lookup[summary, "WeinbergOperatorPresent", False]]];
    AddResult[prefix <> " C5 extraction", Lookup[summary, "WeinbergExtractionStatus", ""] === "Success"];
    AddResult[prefix <> " total Weinberg terms", Lookup[summary, "WeinbergTermCount", -1] === 8];
    AddResult[prefix <> " holomorphic terms", Lookup[summary, "WeinbergHolomorphicTermCount", -1] === 4];
    AddResult[prefix <> " HC terms", Lookup[summary, "WeinbergConjugateTermCount", -1] === 4];
  ];

  c5 = ReadC5[b["Dir"]];
  AddResult[prefix <> " C5 file parses", c5 =!= $Failed];
  If[c5 =!= $Failed,
    AddResult[prefix <> " C5 nonzero", !SafeZeroQ[c5]];

    (* Required topology-closing couplings must all be present. *)
    AddResult[prefix <> " contains lambdaT3", !FreeQ[c5, Coupling[lambdaT3, {}, 0]]];
    AddResult[prefix <> " contains y1", !FreeQ[c5, Bar[Coupling[y1, ___]]]];
    AddResult[prefix <> " contains y2", !FreeQ[c5, Bar[Coupling[y2, ___]]]];

    (* Turning off any required interaction must kill C5. *)
    AddResult[prefix <> " vanishes for lambdaT3=0",
      SafeZeroQ[c5 /. Coupling[lambdaT3, {}, 0] -> 0]];
    AddResult[prefix <> " vanishes for y1=0",
      SafeZeroQ[c5 /. Bar[Coupling[y1, ___]] -> 0]];
    AddResult[prefix <> " vanishes for y2=0",
      SafeZeroQ[c5 /. Bar[Coupling[y2, ___]] -> 0]];

    (* Mass-dimension test: scaling every heavy mass M -> s M should scale C5 -> C5/s. *)
    scaled = c5 /. {
      Coupling[MF, {}, 0] -> scale Coupling[MF, {}, 0],
      Coupling[MS1, {}, 0] -> scale Coupling[MS1, {}, 0],
      Coupling[MS2, {}, 0] -> scale Coupling[MS2, {}, 0]
    };
    dimensionCheck = Quiet @ Check[FullSimplify[scaled - c5/scale, Assumptions -> scale > 0], $Failed];
    AddResult[prefix <> " mass dimension -1", dimensionCheck =!= $Failed && SafeZeroQ[dimensionCheck]];
  ];
,
{benchmark, benchmarks}];

(* ------------------------------------------------------------------- *)
(* Explicit low-dimensional convention bridge. *)
AddResult[
  "T3-B legacy/general lambdaT3 conversion",
  SafeZeroQ[T3LambdaT3GeneralFromLegacyFactor[{2, 2, 1}] + Sqrt[3]/2],
  T3LambdaT3GeneralFromLegacyFactor[{2, 2, 1}]
];

(* T3-B scotogenic calibration                                         *)
(* ------------------------------------------------------------------- *)

c5B = ReadC5["T3_B_alpha_m1"];
If[c5B =!= $Failed,
  (* Remove flavour/group-independent names by mapping the exact Matchete objects
     to ordinary algebraic symbols.  Both Yukawas occur once in the holomorphic C5. *)
  reducedB = c5B /. {
    Bar[Coupling[y1, ___]] -> yy1,
    Bar[Coupling[y2, ___]] -> yy2,
    Coupling[lambdaT3, {}, 0] -> T3LambdaT3GeneralFromLegacyFactor[{2, 2, 1}] lam,
    Coupling[MF, {}, 0] -> mf,
    Coupling[MS1, {}, 0] -> ms1,
    Coupling[MS2, {}, 0] -> ms2
  };

  (* The extracted coefficient should carry precisely hbar*lam*yy1*yy2. *)
  loopB = Quiet @ Check[Cancel[Together[reducedB/(hbar lam yy1 yy2)]], $Failed];
  AddResult["T3-B common coupling factor", loopB =!= $Failed && FreeQ[loopB, lam | yy1 | yy2]];

  If[loopB =!= $Failed,
    (* Degenerate scalar limit ms2 -> ms1.  Use squared-mass variable x to avoid
       ambiguity from the two positive mass parameters approaching one another. *)
    loopX = loopB /. {ms1 -> Sqrt[x1], ms2 -> Sqrt[x2]};
    degenerate = Quiet @ Check[
      FullSimplify[
        Limit[loopX, x2 -> x1],
        Assumptions -> {x1 > 0, mf > 0, x1 != mf^2}
      ],
      $Failed
    ];

    expected = mf/(x1 - mf^2) * (1 - mf^2/(x1 - mf^2) Log[x1/mf^2]);
    difference = If[degenerate === $Failed, $Failed,
      Quiet @ Check[FullSimplify[degenerate - expected,
        Assumptions -> {x1 > 0, mf > 0, x1 != mf^2}], $Failed]
    ];
    AddResult["T3-B scotogenic degenerate loop function",
      difference =!= $Failed && SafeZeroQ[difference],
      If[difference === $Failed, "limit/simplification failed", difference]
    ];
  ];
];

passed = Count[Lookup[results, "Pass"], True];
failed = Length[results] - passed;
report = <|
  "Status" -> If[failed === 0, "PASS", "FAIL"],
  "Passed" -> passed,
  "Failed" -> failed,
  "Tests" -> results
|>;

Export[FileNameJoin[{outputRoot, "t3_regression_report.json"}], report, "RawJSON"];

Print["T3 REGRESSION: ", report["Status"], " (", passed, "/", Length[results], " passed)"];
If[failed > 0,
  Print["Failed checks:"];
  Scan[(If[!TrueQ[#Pass], Print["  - ", #Test, If[#Detail =!= "\"\"", ": " <> #Detail, ""]]]) &, results]
];

Exit[If[failed === 0, 0, 1]];


(* TestT3RGBeta.wl
   Consolidated strict UV RGBeta regression for the five d<=3 T3 benchmarks.
   Replaces the former builder, Yukawa, quartic, representation, closure,
   refinement and diagnostic RGBeta tests.
*)

ClearAll["Global`*"];
Needs["RGBeta`"];

scriptDirectory = DirectoryName @ ExpandFileName[$InputFileName];
projectRoot = ExpandFileName @ FileNameJoin[{scriptDirectory, "..", ".."}];
modelFile = FileNameJoin[{projectRoot, "RGE", "running", "wolfram", "T3RGBetaModel.wl"}];

Fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);
If[!FileExistsQ[modelFile], Fail["T3RGBetaModel.wl not found at " <> modelFile]];
Get[modelFile];

models = {
  {"A", {1,3,2}, 0},
  {"B", {2,2,1}, -1},
  {"C", {2,2,3}, -1},
  {"D", {3,1,2}, -2},
  {"E", {3,3,2}, 0}
};

baseNames = {
  "gY", "g2", "g3", "yu", "yd", "ye", "y1", "y2",
  "MF", "mS1Sq", "mS2Sq",
  "lambdaH", "lambdaS1", "lambdaS2", "lambdaH1", "lambdaH2",
  "lambda12", "lambdaT3"
};

expectedExtras = <|
  "A" -> {"lambdaH2Adj", "lambdaS2Adj"},
  "B" -> {
    "lambdaH1Adj", "lambdaH2Adj", "lambda12Adj",
    "lambdaHHdagS2S2", "lambdaHHdagS1barS1bar", "lambdaS1bar2S2bar2",
    "lambdaS1barS2S2bar2", "lambdaS1S1bar2S2bar",
    "lambdaHHdagS1barS2barCross"
  },
  "C" -> {
    "lambdaH1Adj", "lambdaH2Adj", "lambda12Adj",
    "lambdaHHdagS2S2", "lambdaHHdagS1barS1bar", "lambdaS1bar2S2bar2",
    "lambdaS1barS2S2bar2", "lambdaS1S1bar2S2bar",
    "lambdaHHdagS1barS2barCross"
  },
  "D" -> {"lambdaH1Adj", "lambdaS1Adj"},
  "E" -> {"lambdaH1Adj", "lambdaH2Adj", "lambdaS1Adj", "lambdaS2Adj", "lambda12Adj", "lambda12Cross"}
|>;

HasInternalTensor[expr_] := Length @ Cases[
  expr,
  _RGBeta`PackageScope`QuarticTensors | _RGBeta`PackageScope`UpsilonQuarticTensors,
  Infinity
] > 0;

HasResidualGroupTensor[expr_] := Length @ Cases[
  expr, _del | _delS2 | _fStruct | _tGen | _eps, Infinity
] > 0;

allPassed = True;

Do[
  label = model[[1]];
  dims = model[[2]];
  alpha = model[[3]];
  Print["\nT3-", label, "  dims=", dims, "  alpha=", alpha];

  build = CheckAbort[Quiet @ T3RGBetaBuild[dims[[1]], dims[[2]], dims[[3]], alpha], $Aborted];
  If[!AssociationQ[build],
    Print["  FAIL: model build"];
    allPassed = False;
    Continue[];
  ];

  betas = CheckAbort[Quiet @ T3RGBetaOneLoopBetas[], $Aborted];
  If[!AssociationQ[betas],
    Print["  FAIL: beta association"];
    allPassed = False;
    Continue[];
  ];

  expected = Join[baseNames, expectedExtras[label]];
  missing = Complement[expected, Keys[betas]];
  unexpected = Complement[Keys[betas], expected];

  If[missing =!= {}, Print["  FAIL missing betas: ", missing]; allPassed = False];
  If[unexpected =!= {}, Print["  FAIL unexpected betas: ", unexpected]; allPassed = False];

  KeyValueMap[
    Function[{name, expr},
      If[HasInternalTensor[expr] || HasResidualGroupTensor[expr],
        Print["  FAIL unresolved beta: ", name, " -> ", InputForm[expr]];
        allPassed = False
      ]
    ],
    betas
  ];

  If[missing === {} && unexpected === {} &&
     And @@ (Function[e, !HasInternalTensor[e] && !HasResidualGroupTensor[e]] /@ Values[betas]),
    Print["  PASS: ", Length[betas], " beta functions"]
  ];
,
{model, models}];

If[!TrueQ[allPassed], Fail["one or more T3 RGBeta regressions failed"]];
Print["\nALL T3 RGBETA TESTS PASSED"];
Exit[0];

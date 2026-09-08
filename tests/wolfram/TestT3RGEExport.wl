(* TestT3RGEExport.wl
   Consolidated T3 -> RGE tensor/export regression.
   Covers a simple scotogenic-like benchmark (B) and the richer triplet case (E).
*)

ClearAll["Global`*"];

scriptDirectory = DirectoryName @ ExpandFileName[$InputFileName];
projectRoot = ExpandFileName @ FileNameJoin[{scriptDirectory, "..", ".."}];

Fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);

catalogPath = FileNameJoin[{projectRoot, "wolfram", "t3", "T3ModelCatalog.wl"}];
builderPath = FileNameJoin[{projectRoot, "wolfram", "t3", "LagrangianBuilder.wl"}];
exporterPath = FileNameJoin[{projectRoot, "wolfram", "rge", "T3RGETensorExport.wl"}];

Do[If[!FileExistsQ[path], Fail["missing production file: " <> path]], {path, {catalogPath, builderPath, exporterPath}}];

matcheteLoaded = UsingFrontEnd[Needs["Matchete`"]; True];
If[!TrueQ[matcheteLoaded], Fail["Matchete failed to load"]];

Get[catalogPath];
Get[builderPath];
Get[exporterPath];

benchmarks = {{"B", -1}, {"E", 0}};
allPassed = True;

tmpRoot = CreateDirectory[];

Do[
  label = benchmark[[1]];
  alpha = benchmark[[2]];
  Print["\nT3-", label, " alpha=", alpha];

  model = T3ModelFromClass[label, alpha];
  If[model === $Failed, Print["  FAIL: model construction"]; allPassed = False; Continue[]];

  build = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
  If[!AssociationQ[build] || Lookup[build, "Status", ""] =!= "Success",
    Print["  FAIL: UV build"];
    allPassed = False;
    Continue[];
  ];

  data = CheckAbort[UsingFrontEnd[T3RGEExportAssociation[model]], $Aborted];
  If[!AssociationQ[data], Print["  FAIL: export association"]; allPassed = False; Continue[]];

  requiredKeys = {"quartic_components", "raw_yukawa_components"};
  missingKeys = Select[requiredKeys, !KeyExistsQ[data, #] &];
  If[missingKeys =!= {}, Print["  FAIL missing keys: ", missingKeys]; allPassed = False; Continue[]];

  quarticCount = Length @ Lookup[data, "quartic_components", {}];
  yukawaCount = Length @ Lookup[data, "raw_yukawa_components", {}];
  If[quarticCount <= 0, Print["  FAIL: no quartic components"]; allPassed = False];
  If[yukawaCount <= 0, Print["  FAIL: no Yukawa components"]; allPassed = False];

  outputPath = FileNameJoin[{tmpRoot, "t3_" <> label <> "_rge_exchange.json"}];
  exportResult = Quiet @ Check[ExportT3RGETensors[model, outputPath], $Failed];
  If[exportResult === $Failed || !FileExistsQ[outputPath],
    Print["  FAIL: JSON export"];
    allPassed = False,
    imported = Quiet @ Check[Import[outputPath, "RawJSON"], $Failed];
    If[!AssociationQ[imported], Print["  FAIL: exported JSON cannot be re-imported"]; allPassed = False]
  ];

  If[quarticCount > 0 && yukawaCount > 0, Print["  PASS: ", quarticCount, " quartic, ", yukawaCount, " Yukawa components"]];
,
{benchmark, benchmarks}];

Quiet @ Check[DeleteDirectory[tmpRoot, DeleteContents -> True], Null];

If[!TrueQ[allPassed], Fail["one or more T3 RGE export regressions failed"]];
Print["\nALL T3 RGE EXPORT TESTS PASSED"];
Exit[0];

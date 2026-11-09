(* A(-4) regression for the *installed* full-general scalar-potential builder.
   No in-memory overrides, no source modifications, no old output replacement.
   This checks retention through F-only tree Matchete matching, not RGBeta
   physical coupling normalization or physical Weinberg matching. *)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 1,
  Print["Usage: TestAminus4HermitianPortalProduction.wl OUT_DIR"]; Exit[2]
];
outDir = ExpandFileName[First[args]];
If[FileExistsQ[outDir] || DirectoryQ[outDir],
  Print["Refusing to overwrite: ", outDir]; Exit[3]
];
flags = {"T3_EXPERIMENTAL_GENERAL_FERMIONS",
  "T3_EXPERIMENTAL_GENERAL_QUARTICS",
  "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
  "T3_EXPERIMENTAL_FLAVOR_MASSES"};
If[!And @@ (Environment[#] === "1" & /@ flags),
  Print["All four general flags must equal 1"]; Exit[4]
];
CreateDirectory[outDir, CreateIntermediateDirectories -> True];
sourceDir = DirectoryName @ ExpandFileName[$InputFileName];
Needs["Matchete`"];
Get[FileNameJoin[{sourceDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{sourceDir, "model", "LagrangianBuilder.wl"}]];
Get[FileNameJoin[{sourceDir, "RunMatching.wl"}]];

ClearAll[portalCount, snapshot, writeReport, require];
portalCount[txt_String, inv_Integer] := StringCount[
  txt, "Coupling[lambdaH2Inv" <> ToString[inv] <> ","];
snapshot[expr_, label_String] := Module[{txt, path},
  txt = ToString[InputForm[expr]];
  path = FileNameJoin[{outDir, label <> ".inputform.txt"}];
  Export[path, txt <> "\n", "Text"];
  <|"file" -> FileNameTake[path],
    "lambdaH2Inv1" -> portalCount[txt, 1],
    "lambdaH2Inv2" -> portalCount[txt, 2],
    "characters" -> StringLength[txt]|>
];
writeReport[r_Association] := Export[
  FileNameJoin[{outDir, "production_portal_regression.json"}], r, "RawJSON"];
report = <|
  "model" -> "A_-4",
  "scope" -> "Installed general scalar-potential builder through F-only tree Match",
  "general_flags" -> AssociationMap[Environment, flags],
  "physical_matchete_rgbeta_basis_certified" -> False,
  "physical_gamma_cc_weyl_normalization_certified" -> False,
  "physical_C5_certified" -> False
|>;
require[condition_, status_String, exitCode_Integer] := If[!TrueQ[condition],
  report["status"] = status;
  writeReport[report];
  Print["REGRESSION FAILED: ", status];
  Exit[exitCode]
];
model = T3ModelFromClass["A", -4, 3];
require[AssociationQ[model], "InvalidAMinus4Model", 10];
require[SetT3HeavyFields[{"F"}] =!= $Failed,
  "FOnlyThresholdSetupFailed", 11];
built = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
require[AssociationQ[built], "GeneralUVBuildFailed", 12];
require[Lookup[built, "Status", ""] === "Success" &&
  TrueQ[Lookup[built, "FullValidation", False]] &&
  TrueQ[Lookup[built, "AllScalarExtensionsRequested", False]],
  "GeneralUVValidationFailed", 13];
report["build_status"] = built["Status"];
report["allowed_interactions"] = Lookup[built, "AllowedInteractions", {}];
require[And @@ (MemberQ[report["allowed_interactions"], #] & /@
  {"HiggsPortal2Inv1", "HiggsPortal2Inv2"}),
  "GeneralHiggsS2BasisIncomplete", 14];
snapshots = <|"LUV" -> snapshot[built["LUV"], "LUV"]|>;
report["snapshots"] = snapshots;
require[And @@ (snapshots["LUV"][#] > 0 & /@
  {"lambdaH2Inv1", "lambdaH2Inv2"}), "UVPortalMissing", 15];

matched = CheckAbort[UsingFrontEnd[RunT3Matching[built["LUV"], 5, 0]], $Aborted];
require[AssociationQ[matched] && Lookup[matched, "Status", ""] === "Success",
  "FOnlyMatcheteMatchFailed", 16];
Do[
  snapshots[key] = snapshot[matched[key], key],
  {key, {"RawEFT", "GreenEFT", "EOMEFT", "MatchedEFT"}}
];
report["snapshots"] = snapshots;
require[And @@ (snapshots["MatchedEFT"][#] > 0 & /@
  {"lambdaH2Inv1", "lambdaH2Inv2"}),
  "BothGeneralPortalsNotRetainedInMatchedEFT", 17];
report["status"] = "BothGeneralAminus4PortalsRetained_InventoryOnly";
writeReport[report];
Print["Production general A(-4) portal retention VERIFIED."];
Print["Report: ", FileNameJoin[{outDir,
  "production_portal_regression.json"}]];
Exit[0];

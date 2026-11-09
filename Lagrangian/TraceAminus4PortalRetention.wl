(* Trace both independent Higgs-S2 portals through the ACTUAL A(-4)
   general-Lagrangian builder and the same Matchete tree matching steps as
   RunT3GeneralFThreshold.wl, without changing any production source.

   Run from the repository root in a clean Wolfram kernel:
     wolframscript -file Lagrangian/TraceAminus4PortalRetention.wl \
       output/eft1_general_lagrangian_trace/A_m4_fresh_01

   Requires the FOUR general-source environment toggles to be "1".
   Never overwrites an existing output directory. Does not certify
   physical Matchete-to-RGBeta or GammaCC/Weyl conventions.
*)
ClearAll["Global`*"];

args = Rest[$ScriptCommandLine];
If[Length[args] =!= 1,
  Print["Usage: wolframscript -file Lagrangian/TraceAminus4PortalRetention.wl OUT_DIR"];
  Exit[2]
];
outDir = ExpandFileName[First[args]];
If[FileExistsQ[outDir] || DirectoryQ[outDir],
  Print["Refusing to overwrite existing trace directory: ", outDir];
  Exit[3]
];
requiredFlags = {
  "T3_EXPERIMENTAL_GENERAL_FERMIONS",
  "T3_EXPERIMENTAL_GENERAL_QUARTICS",
  "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
  "T3_EXPERIMENTAL_FLAVOR_MASSES"
};
If[!And @@ (Environment[#] === "1" & /@ requiredFlags),
  Print["All four general-Lagrangian flags must equal 1: ", requiredFlags];
  Exit[4]
];

sourceDir = DirectoryName @ ExpandFileName[$InputFileName];
repoRoot = DirectoryName[sourceDir];
CreateDirectory[outDir, CreateIntermediateDirectories -> True];

(* The same initialization order as the F-only matching driver. *)
Needs["Matchete`"];
Get[FileNameJoin[{sourceDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{sourceDir, "model", "LagrangianBuilder.wl"}]];
Get[FileNameJoin[{sourceDir, "RunMatching.wl"}]];

(* Export the FULL InputForm at every named stage. The first appearance of an
   absent *symbol* is not necessarily a removed invariant: Matchete may hold
   renormalizable operators in effective-coupling form. The explicit compact
   tree and full native invariant tensor must be compared before fixing code. *)
ClearAll[portalSnapshot];
portalSnapshot[label_String, expression_] := Module[{text, path},
  text = ToString[InputForm[expression]];
  path = FileNameJoin[{outDir, label <> ".inputform.txt"}];
  Export[path, text <> "\n", "Text"];
  <|
    "file" -> FileNameTake[path],
    "characters" -> StringLength[text],
    "lambdaH2Inv1_mentions" -> StringCount[text, "lambdaH2Inv1"],
    "lambdaH2Inv2_mentions" -> StringCount[text, "lambdaH2Inv2"],
    "lambdaH2Inv1_as_Coupling" -> StringContainsQ[text, "Coupling[lambdaH2Inv1"],
    "lambdaH2Inv2_as_Coupling" -> StringContainsQ[text, "Coupling[lambdaH2Inv2"]
  |>
];
ClearAll[saveReport];
saveReport[report_Association] := (
  Export[FileNameJoin[{outDir, "portal_retention_trace.json"}], report, "RawJSON"];
  Print["Trace report: ", FileNameJoin[{outDir, "portal_retention_trace.json"}]]
);

(* Compare the old F-only matching snapshot to the CURRENT local source.
   A difference means an old saved F-tree must not be treated as current. *)
sourcePairs = {
  {"Lagrangian/RunT3GeneralFThreshold.wl", "Lagrangian\\RunT3GeneralFThreshold.wl"},
  {"Lagrangian/RunMatching.wl", "Lagrangian\\RunMatching.wl"},
  {"Lagrangian/model/LagrangianBuilder.wl", "Lagrangian\\model\\LagrangianBuilder.wl"},
  {"Lagrangian/interactions/ScalarPotential.wl", "Lagrangian\\interactions\\ScalarPotential.wl"},
  {"Lagrangian/interactions/T3Topology.wl", "Lagrangian\\interactions\\T3Topology.wl"}
};
oldStatusPath = FileNameJoin[{
  repoRoot, "output", "general_observable", "fermion_first_all16_v1",
  "A_m4", "status.json"
}];
oldStatus = If[FileExistsQ[oldStatusPath], Import[oldStatusPath, "RawJSON"], <||>];
savedHashes = Lookup[oldStatus, "source_sha256", <||>];
sourceComparison = Association @ Map[Function[pair,
  Module[{rel, savedKey, currentPath, expected, actual},
    rel = pair[[1]]; savedKey = pair[[2]];
    currentPath = FileNameJoin[{repoRoot, Sequence @@ StringSplit[rel, "/"]}];
    expected = Lookup[savedHashes, savedKey, Missing["NotRecorded"]];
    actual = If[FileExistsQ[currentPath],
      IntegerString[FileHash[currentPath, "SHA256"], 16, 64], Missing["MissingFile"]];
    rel -> <|"saved" -> If[StringQ[expected], expected, Null],
      "current" -> If[StringQ[actual], actual, Null],
      "equal" -> TrueQ[expected === actual]|>
  ]
], sourcePairs];

model = T3ModelFromClass["A", -4, 3];
If[!AssociationQ[model] || SetT3HeavyFields[{"F"}] === $Failed,
  saveReport[<|"status" -> "ModelOrHeavyFieldSetupFailed", "source_hashes" -> sourceComparison|>];
  Exit[5]
];
built = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
If[!AssociationQ[built],
  saveReport[<|"status" -> "LagrangianBuildFailed", "source_hashes" -> sourceComparison|>];
  Exit[6]
];
lagSnap = portalSnapshot["LUV", built["LUV"]];
allowed = Lookup[built, "AllowedInteractions", {}];
sourceReady = (Lookup[built, "Status", ""] === "Success" &&
  TrueQ[Lookup[built, "FullValidation", False]] &&
  TrueQ[Lookup[built, "AllScalarExtensionsRequested", False]] &&
  MemberQ[allowed, "HiggsPortal2Inv1"] &&
  MemberQ[allowed, "HiggsPortal2Inv2"] &&
  lagSnap["lambdaH2Inv1_as_Coupling"] &&
  lagSnap["lambdaH2Inv2_as_Coupling"]);

report = <|
  "model" -> "A_-4",
  "scope" -> "Read-only full-general LUV to Matchete F-only tree retention trace",
  "general_flags" -> AssociationMap[Environment, requiredFlags],
  "source_hashes_vs_saved_run" -> sourceComparison,
  "build_status" -> Lookup[built, "Status", "Missing"],
  "full_validation" -> TrueQ[Lookup[built, "FullValidation", False]],
  "all_scalar_extensions_requested" -> TrueQ[Lookup[built, "AllScalarExtensionsRequested", False]],
  "additional_quartics_expected" -> Lookup[built, "AdditionalQuarticsExpected", Null],
  "additional_quartics_included" -> Lookup[built, "AdditionalQuarticsIncluded", Null],
  "allowed_interactions" -> allowed,
  "both_portals_present_in_assembled_LUV" -> sourceReady,
  "stages" -> <|"LUV" -> lagSnap|>,
  "caveat" -> "Raw/Green/EOM stages may hide native couplings in effective-coupling symbols. " <>
    "A direct-string absence alone does not establish loss of an independent tensor. " <>
    "No physical GammaCC/Weyl or Matchete/RGBeta conversion is certified."
|>;
If[!sourceReady,
  report["status"] = "Blocked_BothPortalsNotFoundInAssembledGeneralLUV";
  saveReport[report];
  Print["STOP: the general LUV is already missing a requested portal"];
  Exit[7]
];
Print["Both independently validated Higgs-S2 portals present in the assembled LUV."];

(* Exactly the tree-stage function used by RunSequentialT3Matching at level 1. *)
matching = CheckAbort[UsingFrontEnd[RunT3Matching[built["LUV"], 5, 0]], $Aborted];
If[!AssociationQ[matching] || Lookup[matching, "Status", ""] =!= "Success",
  report["status"] = "MatcheteTreeMatchingFailed";
  report["matching_status"] = If[AssociationQ[matching],
    Lookup[matching, "Status", "Missing"], ToString[InputForm[matching]]];
  saveReport[report];
  Exit[8]
];

snapshots = report["stages"];
Do[
  snapshots[key] = portalSnapshot[key, matching[key]],
  {key, {"RawEFT", "GreenEFT", "EOMEFT", "LoopEFT", "MatchedEFT"}}
];
report["stages"] = snapshots;
compact = CompactExplicitEFT[matching];
transition = CanonicalTransitionEFT[matching];
If[compact === $Failed || transition === $Failed,
  report["status"] = "CompactOrTransitionEFTFailed";
  saveReport[report];
  Exit[9]
];
snapshots = report["stages"];
snapshots["CompactExplicitEFT"] = portalSnapshot["CompactExplicitEFT", compact];
snapshots["CanonicalTransitionEFT"] = portalSnapshot["CanonicalTransitionEFT", transition];
report["stages"] = snapshots;

oldMatchedPath = FileNameJoin[{repoRoot, "output", "general_observable",
  "fermion_first_all16_v1", "A_m4", "data", "eft1_after_F_tree_inputform.txt"}];
If[FileExistsQ[oldMatchedPath],
  oldText = Import[oldMatchedPath, "Text"];
  report["old_saved_match"] = <|
    "file" -> oldMatchedPath,
    "lambdaH2Inv1_mentions" -> StringCount[oldText, "lambdaH2Inv1"],
    "lambdaH2Inv2_mentions" -> StringCount[oldText, "lambdaH2Inv2"]|>
];
final1 = report["stages"]["CompactExplicitEFT"]["lambdaH2Inv1_mentions"] > 0;
final2 = report["stages"]["CompactExplicitEFT"]["lambdaH2Inv2_mentions"] > 0;
report["status"] = Which[
  final1 && final2, "BothPortalsRetainedInFreshCompactEFT",
  !final1 && final2, "Inv1AbsentFromFreshCompactEFT_InspectIntermediateStages",
  final1 && !final2, "Inv2AbsentFromFreshCompactEFT_InspectIntermediateStages",
  True, "BothPortalSymbolsAbsentFromFreshCompactEFT_InspectIntermediateStages"
];
report["needs_native_tensor_basis_check"] = True;
report["physical_c5_certified"] = False;
saveReport[report];
Print["Status: ", report["status"]];
Print["Fresh compact EFT: lambdaH2Inv1 mentions=",
  report["stages"]["CompactExplicitEFT"]["lambdaH2Inv1_mentions"],
  ", lambdaH2Inv2 mentions=",
  report["stages"]["CompactExplicitEFT"]["lambdaH2Inv2_mentions"]];
Exit[0];

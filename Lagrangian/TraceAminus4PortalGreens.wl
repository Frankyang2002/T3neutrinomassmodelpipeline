(*
  Read-only A(-4) general-portal source isolation at Matchete GreensSimplify.
  Does not modify production matching code or assume a renormalisable
  Matchete/RGBeta basis dictionary. Each output directory must be new.

  PowerShell (repository root):
    $env:T3_EXPERIMENTAL_GENERAL_FERMIONS = '1'
    $env:T3_EXPERIMENTAL_GENERAL_QUARTICS = '1'
    $env:T3_EXPERIMENTAL_GENERAL_QUADRATICS = '1'
    $env:T3_EXPERIMENTAL_FLAVOR_MASSES = '1'
    wolframscript -file Lagrangian\TraceAminus4PortalGreens.wl `
      output\eft1_general_lagrangian_trace\A_m4_greens_01
*)
ClearAll["Global`*"];

args = Rest[$ScriptCommandLine];
If[Length[args] =!= 1,
  Print["Usage: wolframscript -file Lagrangian/TraceAminus4PortalGreens.wl OUT_DIR"];
  Exit[2]
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
  Print["All four general-Lagrangian switches must equal 1."]; Exit[4]
];
CreateDirectory[outDir, CreateIntermediateDirectories -> True];
sourceDir = DirectoryName @ ExpandFileName[$InputFileName];

Needs["Matchete`"];
Get[FileNameJoin[{sourceDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{sourceDir, "model", "LagrangianBuilder.wl"}]];
Get[FileNameJoin[{sourceDir, "RunMatching.wl"}]];

ClearAll[save, writeReport, termList, termsWithCoupling, attemptGreens];
save[label_String, expr_] := Module[{s, path},
  s = ToString[InputForm[expr]];
  path = FileNameJoin[{outDir, label <> ".inputform.txt"}];
  Export[path, s <> "\n", "Text"];
  <|"file" -> FileNameTake[path], "characters" -> StringLength[s],
    "inv1_mentions" -> StringCount[s, "lambdaH2Inv1"],
    "inv2_mentions" -> StringCount[s, "lambdaH2Inv2"],
    "inv1_as_explicit_coupling" -> StringContainsQ[s, "Coupling[lambdaH2Inv1"],
    "inv2_as_explicit_coupling" -> StringContainsQ[s, "Coupling[lambdaH2Inv2"],
    "exactly_zero" -> TrueQ[expr === 0]|>
];
writeReport[report_Association] := (
  Export[FileNameJoin[{outDir, "greens_portal_isolation.json"}],
    report, "RawJSON"];
  Print["Report: ", FileNameJoin[{outDir, "greens_portal_isolation.json"}]]
);
termList[expr_] := Module[{e = Expand[expr]},
  If[Head[e] === Plus, List @@ e, {e}]
];
(* Select whole expanded summands; retain CGs and every flavour index. *)
termsWithCoupling[expr_, couplingName_String] := Select[
  termList[expr],
  StringContainsQ[ToString[InputForm[#]], "Coupling[" <> couplingName <> ","] &
];
attemptGreens[expr_] := CheckAbort[UsingFrontEnd[GreensSimplify[expr]], $Aborted];

model = T3ModelFromClass["A", -4, 3];
If[!AssociationQ[model] || SetT3HeavyFields[{"F"}] === $Failed,
  Print["Cannot construct model or set all-fermion-first threshold"]; Exit[5]
];
built = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
If[!AssociationQ[built] || Lookup[built, "Status", ""] =!= "Success" ||
  !TrueQ[Lookup[built, "FullValidation", False]] ||
  !TrueQ[Lookup[built, "AllScalarExtensionsRequested", False]],
  Print["General UV Lagrangian did not pass its existing source checks"]; Exit[6]
];
allowed = Lookup[built, "AllowedInteractions", {}];
If[!And @@ (MemberQ[allowed, #] & /@
    {"HiggsPortal2Inv1", "HiggsPortal2Inv2"}),
  Print["One of the two independently validated UV portals is not assembled"];
  Exit[7]
];

(* Independently inspect the two CG arrays. This does NOT check their
   normalization against RGBeta or establish their operator-level map. *)
cgInfo = <|"status" -> "Unavailable"|>;
cg1 = Quiet @ Check[Matchete`CGManipulations`PackagePrivate`$CGtensors[
  T3HiggsPortal2CGInv1], $Failed];
cg2 = Quiet @ Check[Matchete`CGManipulations`PackagePrivate`$CGtensors[
  T3HiggsPortal2CGInv2], $Failed];
If[cg1 =!= $Failed && cg2 =!= $Failed &&
  ArrayQ[Normal[cg1]] && ArrayQ[Normal[cg2]] &&
  Dimensions[Normal[cg1]] === Dimensions[Normal[cg2]],
  cgInfo = <|
    "status" -> "TensorArraysRead",
    "dimensions" -> Dimensions[Normal[cg1]],
    "linear_rank" -> MatrixRank[{Flatten[Normal[cg1]], Flatten[Normal[cg2]]}]
  |>
];

matching = CheckAbort[UsingFrontEnd[RunT3Matching[built["LUV"], 5, 0]],
  $Aborted];
If[!AssociationQ[matching] || Lookup[matching, "Status", ""] =!= "Success",
  writeReport[<|"model" -> "A_-4", "status" -> "MatchingFailed",
    "cg_tensors" -> cgInfo|>]; Exit[8]
];
raw = matching["RawEFT"];
green = matching["GreenEFT"];
rawTerms = termList[raw];
inv1Terms = termsWithCoupling[raw, "lambdaH2Inv1"];
inv2Terms = termsWithCoupling[raw, "lambdaH2Inv2"];
raw1 = Expand[Total[inv1Terms]];
raw2 = Expand[Total[inv2Terms]];
rawWithout1 = Expand[raw - raw1];

results = <|
  "assembled_LUV" -> save["LUV", built["LUV"]],
  "raw_full" -> save["RawEFT", raw],
  "greens_full" -> save["GreenEFT", green],
  "raw_inv1_only" -> save["RawEFT_lambdaH2Inv1_only", raw1],
  "raw_inv2_only" -> save["RawEFT_lambdaH2Inv2_only", raw2],
  "raw_without_inv1" -> save["RawEFT_without_lambdaH2Inv1", rawWithout1]
|>;
report = <|
  "model" -> "A_-4",
  "status" -> "PendingSeparateGreensProbes",
  "scope" -> "Independent linear-source probes of Matchete GreensSimplify",
  "general_flags" -> AssociationMap[Environment, flags],
  "source_validation_passed" -> True,
  "both_portals_in_original_general_LUV" -> True,
  "cg_tensors" -> cgInfo,
  "raw_summand_counts" -> <|
    "all" -> Length[rawTerms],
    "inv1" -> Length[inv1Terms],
    "inv2" -> Length[inv2Terms]|>,
  "snapshots" -> results,
  "warning" -> "Coupling-name disappearance need not mean operator loss. " <>
    "Probes preserve CGs and indices but are not physical basis certificates.",
  "physical_C5_certified" -> False
|>;
If[Length[inv1Terms] == 0 || Length[inv2Terms] == 0,
  report["status"] = "Blocked_RawEFTMissingExplicitPortalTerms";
  writeReport[report]; Exit[9]
];

(* Use the same native simplifier as the verified full matching chain.
   Isolating source coefficients distinguishes an identity that already
   kills one invariant from full-input cross-term simplification. *)
Print["Testing GreensSimplify on isolated lambdaH2Inv1 source..."];
green1 = attemptGreens[raw1];
Print["Testing GreensSimplify on isolated lambdaH2Inv2 source..."];
green2 = attemptGreens[raw2];
Print["Testing GreensSimplify on RawEFT without lambdaH2Inv1..."];
greenWithout1 = attemptGreens[rawWithout1];

If[MemberQ[{$Failed, $Aborted}, green1] ||
  MemberQ[{$Failed, $Aborted}, green2] ||
  MemberQ[{$Failed, $Aborted}, greenWithout1],
  report["status"] = "OneOrMoreIsolatedGreensCallsFailed";
  report["isolated_Greens_status"] = <|
    "inv1_failed" -> MemberQ[{$Failed, $Aborted}, green1],
    "inv2_failed" -> MemberQ[{$Failed, $Aborted}, green2],
    "without_inv1_failed" -> MemberQ[{$Failed, $Aborted}, greenWithout1]|>;
  writeReport[report]; Exit[10]
];

results["greens_inv1_only"] = save["GreenEFT_lambdaH2Inv1_only", green1];
results["greens_inv2_only"] = save["GreenEFT_lambdaH2Inv2_only", green2];
results["greens_without_inv1"] = save["GreenEFT_without_lambdaH2Inv1", greenWithout1];
(* Compare in the same syntactic basis, without declaring mathematical
   equivalence of expressions with differently named effective couplings. *)
delta = Expand[green - greenWithout1];
results["greens_full_minus_greens_without_inv1"] =
  save["GreenEFT_difference_with_without_inv1", delta];
report["snapshots"] = results;
report["isolated_inv1_greens_exactly_zero"] = TrueQ[green1 === 0];
report["full_minus_without_inv1_exactly_zero"] = TrueQ[delta === 0];
report["status"] = Which[
  TrueQ[green1 === 0], "Inv1IsolatedGreensZero_CheckSU2CGIdentity",
  TrueQ[delta === 0], "FullGreensIndependentOfInv1_CheckRebasisOrCancellation",
  True, "Inv1SurvivesIsolatedReduction_CheckFullRecombination"
];
writeReport[report];
Print["Status: ", report["status"]];
Print["CG tensor rank: ", Lookup[cgInfo, "linear_rank", "Unavailable"]];
Print["Raw inv1 terms: ", Length[inv1Terms],
  "; raw inv2 terms: ", Length[inv2Terms]];
Exit[0];

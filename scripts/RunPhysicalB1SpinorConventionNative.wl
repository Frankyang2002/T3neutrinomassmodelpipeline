(* Restricted B1 spinor convention: interrogate the INSTALLED Matchete action.
   This is not an invented numerical GammaCC and does not call Match/EFT/RGE.
   Usage: wolframscript -file scripts/RunPhysicalB1SpinorConventionNative.wl OUTPUT.json
   A previous output is never overwritten.
*)
ClearAll["Global`*"];
$bcArgs = Rest[$ScriptCommandLine];
If[Length[$bcArgs] =!= 1, Print["Usage: OUTPUT.json"]; Exit[2]];
$bcOut = ExpandFileName[First[$bcArgs]];
If[FileExistsQ[$bcOut], Print["BLOCKED: refusing to overwrite: ", $bcOut]; Exit[2]];
$bcRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$bcFlags = {"T3_EXPERIMENTAL_GENERAL_FERMIONS", "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
  "T3_EXPERIMENTAL_GENERAL_QUARTICS", "T3_EXPERIMENTAL_FLAVOR_MASSES"};
$bcPrior = AssociationMap[Function[x, With[{v = Environment[x]}, If[StringQ[v], v, "UNSET"]]], $bcFlags];
SetEnvironment[Thread[$bcFlags -> ConstantArray["0", Length[$bcFlags]]]];
If[!AllTrue[$bcFlags, Environment[#] === "0" &], Print["BLOCKED: flags are not zero"]; Exit[3]];
If[CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]]; True, $Failed], $Aborted] =!= True,
 Print["BLOCKED: installed Matchete unavailable"]; Exit[4]];
Get[FileNameJoin[{$bcRoot, "Lagrangian", "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{$bcRoot, "Lagrangian", "model", "LagrangianBuilder.wl"}]];
ClearAll[BCFail, BCText, BCTerms, BCNumber, BCPick, BCRun, BCSourceHashes];
BCFail[msg_String] := Throw[<|"status" -> "blocked", "reason" -> msg|>, "BCStop"];
BCText[x_] := ToString[x, InputForm];
BCTerms[e_, strings_List] := Select[
  If[Head[Expand[e]] === Plus, List @@ Expand[e], {Expand[e]}],
  Function[t, With[{s = BCText[t]}, AllTrue[strings, StringContainsQ[s, #] &]]]];
BCNumber[t_] := Times @@ Select[If[Head[t] === Times, List @@ t, {t}], NumberQ];
BCPick[e_, tokens_List, n_Integer, expected_List] := Module[{items, coeff},
 items = BCTerms[e, tokens]; coeff = BCNumber /@ items;
 If[Length[items] =!= n || coeff =!= expected,
  BCFail["Unexpected native terms for " <> BCText[tokens] <> ": coefficients=" <> BCText[coeff]]];
 <|"coefficients" -> (BCText /@ coeff), "native_terms_inputform" -> (BCText /@ items)|>];
BCSourceHashes[] := Module[{paths, f},
 paths = {"Lagrangian/model/T3ModelCatalog.wl", "Lagrangian/model/LagrangianBuilder.wl",
   "Lagrangian/model/T3Fields.wl", "Lagrangian/model/SU2Invariants.wl",
   "Lagrangian/interactions/T3Topology.wl", "Lagrangian/interactions/ScalarPotential.wl"};
 Association@Table[
   f = FileNameJoin[Prepend[FileNameSplit[path], $bcRoot]];
   If[!FileExistsQ[f], BCFail["Required source missing: " <> path]];
   path -> IntegerString[FileHash[f, "SHA256"], 16, 64], {path, paths}]];
BCRun[] := Module[{cases = <||>, specs, id, family, alpha, built, lf, lint,
  y1orig, y2orig, y1hc, y2hc, mass, checks, sample, explicit4, sym, gamma,
  sources, gith, frc, info},
 specs = {{"A_p0", "A", 0}, {"B_m1", "B", -1}, {"C_m1", "C", -1}};
 Do[
  {id, family, alpha} = entry;
  built = CheckAbort[Check[UsingFrontEnd[BuildT3Lagrangian[T3ModelFromClass[family, alpha]]],
    $Failed], $Aborted];
  If[!AssociationQ[built] || Lookup[built, "Status", ""] =!= "Success" ||
      !TrueQ[Lookup[built, "FullValidation", False]],
    BCFail["Native validated action unavailable: " <> id]];
  lf = built["LFree"]; lint = built["LInt"];
  y1orig = BCPick[lint, {"Coupling[y1,", "CG[T3Y1CG,"}, 1, {1}];
  y2orig = BCPick[lint, {"Coupling[y2,", "CG[T3Y2CG,"}, 1, {1}];
  y1hc = BCPick[lint, {"Bar[Coupling[y1,", "CG[Bar[T3Y1CG],"}, 1, {1}];
  y2hc = BCPick[lint, {"Bar[Coupling[y2,", "CG[Bar[T3Y2CG],"}, 1, {1}];
  mass = BCPick[lf, {"Coupling[MF,", "NewFermion"}, 1,
    If[id === "A_p0", {-1}, {-1/2}]];
  checks = <|
    "native_assembled_action_validated" -> True,
    "yukawa_original_and_hc_each_plus_one" -> True,
    "dirac_versus_majorana_mass_prefactor" -> True,
    "y1_original_has_GammaCC" -> If[id === "A_p0",
      StringContainsQ[First[y1orig["native_terms_inputform"]], "GammaCC"],
      !StringContainsQ[First[y1orig["native_terms_inputform"]], "GammaCC"]],
    "y1_hc_GammaCC_placement" -> If[id === "A_p0",
      StringContainsQ[First[y1hc["native_terms_inputform"]], "GammaCC"],
      StringContainsQ[First[y1hc["native_terms_inputform"]], "GammaCC"]],
    "y2_original_without_GammaCC" -> !StringContainsQ[
      First[y2orig["native_terms_inputform"]], "GammaCC"],
    "majorana_mass_has_GammaCC" -> If[id === "A_p0", True,
      StringContainsQ[First[mass["native_terms_inputform"]], "GammaCC"]],
    "restricted_flags_zero" -> AllTrue[$bcFlags, Environment[#] === "0" &]
    |>;
  If[!AllTrue[Values[checks], TrueQ], BCFail["Native spinor-operator check failed: " <> id <> " " <> BCText[checks]]];
  cases[id] = <|"native_action_passed" -> True, "mass" -> mass,
    "y1_original" -> y1orig, "y2_original" -> y2orig,
    "y1_hc" -> y1hc, "y2_hc" -> y2hc,
    "projector_strings" -> (BCText /@ Cases[{lf, lint}, HoldPattern[DiracProduct[Proj[_]]], Infinity]),
    "checks" -> checks|>;
 , {entry, specs}];
 (* CRITICAL: GammaCC is normally an abstract Clifford-algebra object,
    NOT an installed concrete four-by-four matrix. Do not infer a matrix
    phase from an empty OwnValues or a textbook choice of C. *)
 sym = GammaCC;
 explicit4 = MatrixQ[sym, NumericQ] && Dimensions[sym] === {4, 4};
 gamma = <|
   "symbol" -> BCText[Unevaluated[GammaCC]],
   "symbol_context" -> Context[GammaCC],
   "own_values_count" -> Length[OwnValues[GammaCC]],
   "down_values_count" -> Length[DownValues[GammaCC]],
   "up_values_count" -> Length[UpValues[GammaCC]],
   "explicit_numeric_4x4_exposed" -> TrueQ[explicit4],
   "installed_absolute_phase_established" -> False,
   "note" -> "Matchete symbolic GammaCC identities are not a numerical representation/phase."|>;
 frc = Names["*`FeynmanRules"];
 info = <|"feynman_rules_symbol_names" -> frc,
   "feynman_rules_called" -> False,
   "reason" -> "Capability is recorded; no undocumented call is treated as phase proof."|>;
 sources = BCSourceHashes[];
 gith = Quiet@Check[RunProcess[{"git", "-C", $bcRoot, "rev-parse", "HEAD"}, "StandardOutput"], "UNAVAILABLE"];
 <|"format" -> "t3-b1-native-spinor-convention-v1",
  "status" -> "LiveNativeB1SpinorOperatorsAudited_AbsoluteGammaCCPhaseUnresolved",
  "models" -> cases, "source_sha256" -> sources, "source_git_head" -> StringTrim[gith],
  "restricted_environment" -> <|"prior" -> $bcPrior,
    "effective" -> AssociationMap[Environment, $bcFlags]|>,
  "installed_gammaCC" -> gamma, "optional_feynman_rules_capability" -> info,
  "matchete_native_relative_sign_from_concrete_spinors_measured" -> False,
  "chosen_Weyl_convention_derived_in_this_script" -> False,
  "physical_source_certified" -> False,
  "matching_or_running_performed" -> False|>
];
$bcResult = CheckAbort[Catch[BCRun[], "BCStop"], $Aborted];
If[!AssociationQ[$bcResult] || Lookup[$bcResult, "status", ""] =!=
  "LiveNativeB1SpinorOperatorsAudited_AbsoluteGammaCCPhaseUnresolved",
 Print["BLOCKED: ", If[AssociationQ[$bcResult], Lookup[$bcResult, "reason", "native failure"], "abort"]]; Exit[5]];
If[!DirectoryQ[DirectoryName[$bcOut]],
 CreateDirectory[DirectoryName[$bcOut], CreateIntermediateDirectories -> True]];
If[Quiet@Check[Export[$bcOut, $bcResult, "RawJSON"], $Failed] === $Failed,
 Print["BLOCKED: JSON export failed"]; Exit[6]];
Print["A_+0, B_-1, C_-1: live native spinor operator and mass-factor checks PASSED"];
Print["GammaCC concrete 4x4 exposed: ", $bcResult["installed_gammaCC", "explicit_numeric_4x4_exposed"]];
Print["Absolute native GammaCC phase: UNRESOLVED (not inferred from chosen Weyl basis)"];
Print["Witness: ", $bcOut];
Print["Physical source-certified: 0/16"];
Exit[0];

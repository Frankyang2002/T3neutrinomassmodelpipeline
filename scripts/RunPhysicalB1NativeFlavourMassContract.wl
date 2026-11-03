(* Read-only native B1 flavour-mass action check, 2026-11-03.
   Usage:
     wolframscript -file scripts/RunPhysicalB1NativeFlavourMassContract.wl OUTPUT.json
   Explicitly switches on ONLY the experimental flavour-diagonal mass flag.
   No Match, EFT, RGBeta, running, or production-source writes. *)
ClearAll["Global`*"];
$fmArgs = Rest[$ScriptCommandLine];
If[Length[$fmArgs] =!= 1,
  Print["Usage: OUTPUT.json"]; Exit[2]];
$fmOut = ExpandFileName[First[$fmArgs]];
If[FileExistsQ[$fmOut],
  Print["BLOCKED: would overwrite existing witness: ", $fmOut]; Exit[2]];
$fmRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$fmFlags = {"T3_EXPERIMENTAL_GENERAL_FERMIONS",
  "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
  "T3_EXPERIMENTAL_GENERAL_QUARTICS",
  "T3_EXPERIMENTAL_FLAVOR_MASSES"};
$fmPrior = AssociationMap[
  Function[x, With[{v = Environment[x]}, If[StringQ[v], v, "UNSET"]]],
  $fmFlags];
SetEnvironment[Thread[$fmFlags -> {"0", "0", "0", "1"}]];
If[(Environment /@ $fmFlags) =!= {"0", "0", "0", "1"},
  Print["BLOCKED: flag isolation failed"]; Exit[3]];
If[CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]]; True, $Failed],
    $Aborted] =!= True,
  Print["BLOCKED: installed Matchete unavailable"]; Exit[4]];
Get[FileNameJoin[{$fmRoot, "Lagrangian", "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{$fmRoot, "Lagrangian", "model", "LagrangianBuilder.wl"}]];
ClearAll[FMText, FMTerms, FMCoef, FMFail, FMHash, FMRun];
FMText[x_] := ToString[x, InputForm];
FMTerms[x_] := Module[{e = Expand[x]},
  If[Head[e] === Plus, List @@ e, {e}]];
FMCoef[x_] := Times @@ Select[
  If[Head[x] === Times, List @@ x, {x}], NumberQ];
FMFail[msg_String] := Throw[
  <|"status" -> "BlockedNativeFlavourMassContract",
    "reason" -> msg, "physical_source_certified" -> False|>, "FMStop"];
FMHash[] := Module[{paths, path, file},
  paths = {"Lagrangian/model/T3ModelCatalog.wl",
    "Lagrangian/model/LagrangianBuilder.wl",
    "Lagrangian/model/T3Fields.wl",
    "Lagrangian/model/SU2Invariants.wl",
    "Lagrangian/interactions/T3Topology.wl",
    "Lagrangian/interactions/ScalarPotential.wl"};
  Association@Table[
    file = FileNameJoin[Prepend[FileNameSplit[path], $fmRoot]];
    If[!FileExistsQ[file], FMFail["Missing source file: " <> path]];
    path -> IntegerString[FileHash[file, "SHA256"], 16, 64],
    {path, paths}]];
FMRun[] := Module[{specs, rows = <||>, id, family, alpha, model, built,
    free, terms, mass, coeff, couplings, ix, fields, fieldFlavours,
    checks, sourceFile, sourceText, hashes, git},
  sourceFile = FileNameJoin[{$fmRoot, "Lagrangian", "model", "T3Fields.wl"}];
  sourceText = Import[sourceFile, "Text"];
  If[!StringContainsQ[sourceText, "T3ExperimentalFlavorMassesQ"] ||
      !StringContainsQ[sourceText, "{kind, mass, {NFlavor}}"],
    FMFail["Source no longer declares flavour-indexed mass specification"]];
  hashes = FMHash[];
  specs = {{"A_p0", "A", 0}, {"B_m1", "B", -1},
    {"C_m1", "C", -1}};
  Do[
    {id, family, alpha} = spec;
    model = T3ModelFromClass[family, alpha];
    If[!AssociationQ[model], FMFail[id <> ": model catalogue failed"]];
    built = CheckAbort[
      Check[UsingFrontEnd[BuildT3Lagrangian[model]], $Failed], $Aborted];
    If[!AssociationQ[built],
      FMFail[id <> ": Matchete did not return an action association"]];
    If[Lookup[built, "Status", ""] =!= "Success" ||
        !TrueQ[Lookup[built, "FullValidation", False]],
      FMFail[id <> ": native action did not validate with flavour masses; status=" <>
        FMText[Lookup[built, "Status", "UNAVAILABLE"]]]];
    free = built["LFree"];
    terms = FMTerms[free];
    mass = Select[terms,
      Function[t, !FreeQ[t, HoldPattern[Coupling[MF, ___]]]]];
    If[Length[mass] =!= 1,
      FMFail[id <> ": expected one canonical MF mass monomial; found " <>
        ToString[Length[mass]]]];
    coeff = FMCoef[First[mass]];
    If[coeff =!= If[id === "A_p0", -1, -1/2],
      FMFail[id <> ": unexpected Dirac/Majorana mass coefficient " <>
        FMText[coeff]]];
    couplings = Cases[First[mass],
      HoldPattern[Coupling[MF, inds_List, ___]] :> inds, Infinity];
    If[Length[couplings] =!= 1 || Length[First[couplings]] =!= 1 ||
        !MatchQ[First[couplings][[1]], HoldPattern[Index[_, NFlavor]]],
      FMFail[id <> ": MF is not an exactly one-index flavour mass; coupling=" <>
        FMText[couplings]]];
    ix = First[couplings][[1]];
    fields = Cases[First[mass],
      HoldPattern[Field[NewFermion, Fermion, inds_List, derivs_List]] :> inds,
      Infinity];
    If[Length[fields] =!= 2,
      FMFail[id <> ": mass term does not contain two fermion field factors"]];
    fieldFlavours = Cases[#, HoldPattern[Index[_, NFlavor]], Infinity] & /@ fields;
    If[!AllTrue[fieldFlavours, # === {ix} &],
      FMFail[id <> ": MF index is not the same as BOTH fermion-flavour indices"]];
    checks = <|
      "native_action_accepted" -> True,
      "exactly_one_MF_mass_monomial" -> True,
      "correct_Dirac_or_Majorana_mass_coefficient" -> True,
      "MF_has_exactly_one_NFlavor_index" -> True,
      "both_fermion_fields_share_MF_flavour_index" -> True,
      "other_three_experimental_flags_disabled" ->
        AllTrue[Take[$fmFlags, 3], Environment[#] === "0" &],
      "flavour_mass_flag_enabled" ->
        Environment[Last[$fmFlags]] === "1"|>;
    If[!AllTrue[Values[checks], TrueQ],
      FMFail[id <> ": executable checks failed"]];
    rows[id] = <|
      "checks" -> checks,
      "mass_prefactor" -> FMText[coeff],
      "mass_term_inputform" -> FMText[First[mass]],
      "mass_coupling_index_inputform" -> FMText[ix],
      "source_MF_is_one_symbol_with_index_not_a_general_flavour_matrix" -> True,
      "diagonal_mass_eigenvalue_interpretation_native" -> True,
      "three_numeric_MF_values_assigned_inside_Matchete" -> False|>;
    Print[id, ": flavour-indexed MF native action VERIFIED (", FMText[coeff], ")"];
    , {spec, specs}];
  git = Quiet@Check[RunProcess[{"git", "-C", $fmRoot, "rev-parse", "HEAD"},
    "StandardOutput"], "UNAVAILABLE"];
  <|
    "format" -> "t3-b1-live-native-flavour-diagonal-mass-contract-v1",
    "status" -> "NativeThreeModelFlavourMassActionAccepted",
    "models" -> rows, "source_sha256" -> hashes,
    "local_git_head" -> StringTrim[git],
    "experimental_flags_prior" -> $fmPrior,
    "experimental_flags_effective" -> AssociationMap[Environment, $fmFlags],
    "general_offdiagonal_MF_matrix_certified" -> False,
    "three_numeric_MF_values_assigned_inside_Matchete" -> False,
    "benchmark_MF_diagonal_extension_supported_by_native_action_form" -> True,
    "installed_GammaCC_numerical_phase_measured" -> False,
    "absolute_cross_package_spinor_phase_certified" -> False,
    "matching_or_running_performed" -> False,
    "physical_source_certified" -> False,
    "certified_models" -> 0|>
];
$fmResult = CheckAbort[Catch[FMRun[], "FMStop"], $Aborted];
If[!AssociationQ[$fmResult],
  $fmResult = <|"status" -> "BlockedNativeFlavourMassContract",
    "reason" -> "Matchete evaluation aborted/non-association",
    "physical_source_certified" -> False|>];
If[!DirectoryQ[DirectoryName[$fmOut]],
  CreateDirectory[DirectoryName[$fmOut], CreateIntermediateDirectories -> True]];
If[Quiet@Check[Export[$fmOut, $fmResult, "RawJSON"], $Failed] === $Failed,
  Print["BLOCKED: unable to export JSON witness"]; Exit[7]];
Print["Witness: ", $fmOut];
Print["Status: ", Lookup[$fmResult, "status", "Blocked"]];
If[Lookup[$fmResult, "status", ""] =!= "NativeThreeModelFlavourMassActionAccepted",
  Print["Reason: ", Lookup[$fmResult, "reason", "unknown"]]; Exit[6]];
Print["Flavour diagonal MF action: 3/3 ACCEPTED; no EFT, no source changes"];
Print["Numerical MF substitution not performed; physical certification: 0/16"];
Exit[0];

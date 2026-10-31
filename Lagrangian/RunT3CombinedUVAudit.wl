(* Workstream A: combined general-fermion and general-scalar UV audit.
   Project path: Lagrangian/RunT3CombinedUVAudit.wl

   wolframscript -file Lagrangian/RunT3CombinedUVAudit.wl \
       OUTPUT.json CLASS alpha MODE

   MODE = quartics: general Yukawas + three independent diagonal fermion
     masses + all extra scalar quartics; no additional scalar bilinears.
   MODE = full: add all allowed extra scalar bilinears.  B/C(-1) are
     expected to be blocked by Matchete's non-diagonal heavy scalar basis.

   This is an assembly/CheckLagrangian audit only: no Match or RGBeta.
   It is not a physical UV-operator-completeness certificate.
*)

ClearAll["Global`*"];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 4,
  Print["Usage: OUTPUT.json CLASS alpha quartics|full"];
  Exit[2]
];

out = ExpandFileName[args[[1]]];
cls = ToUpperCase[args[[2]]];
mode = ToLowerCase[args[[4]]];
T3ParseCombinedAlpha[s_String] := Which[
  StringMatchQ[s, "m" ~~ DigitCharacter ..], -ToExpression[StringDrop[s, 1]],
  StringMatchQ[s, "p" ~~ DigitCharacter ..], ToExpression[StringDrop[s, 1]],
  StringMatchQ[s, ("-" | "") ~~ DigitCharacter ..], ToExpression[s],
  True, $Failed
];
alpha = T3ParseCombinedAlpha[args[[3]]];
If[alpha === $Failed || !MemberQ[{"quartics", "full"}, mode],
  Print["ERROR: invalid alpha or audit mode."];
  Exit[2]
];

(* Independent fixed 16-model structure counts from the Workstream A
   SU(2)-projected catalogue; do not infer expected counts from the builder. *)
T3CombinedExpected = <|
  "A_m4" -> {0, 0}, "A_m2" -> {5, 1}, "A_p0" -> {5, 1},
  "B_m3" -> {1, 0}, "B_m1" -> {7, 1}, "B_p1" -> {1, 0},
  "C_m3" -> {1, 0}, "C_m1" -> {7, 1}, "C_p1" -> {1, 0},
  "D_m2" -> {5, 1}, "D_p0" -> {5, 1}, "D_p2" -> {0, 0},
  "E_m4" -> {0, 0}, "E_m2" -> {6, 1},
  "E_p0" -> {6, 1}, "E_p2" -> {0, 0}
|>;
key = cls <> "_" <> If[alpha < 0, "m", "p"] <> ToString[Abs[alpha]];
If[!KeyExistsQ[T3CombinedExpected, key],
  Print["ERROR: unsupported model ", key];
  Exit[2]
];

expectedQuartics = T3CombinedExpected[key][[1]];
expectedBilinears = If[mode === "full", T3CombinedExpected[key][[2]], 0];
expectedYukawas = If[MemberQ[
  {"A_m2", "A_p0", "D_m2", "D_p0", "E_m2", "E_p0"}, key], 1, 0];
expectedMixedMassBlock = mode === "full" &&
  MemberQ[{"B_m1", "C_m1"}, key];

(* Fail closed if a workstream changes an inherited environment variable. *)
If[Environment["T3_EXPERIMENTAL_FLAVOR_MASSES"] =!= "1" ||
   Environment["T3_EXPERIMENTAL_GENERAL_FERMIONS"] =!= "1" ||
   Environment["T3_EXPERIMENTAL_GENERAL_QUARTICS"] =!= "1" ||
   (Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] === "1") =!=
     (mode === "full"),
  Print["ERROR: set FLAVOR_MASSES=1, GENERAL_FERMIONS=1, ",
    "GENERAL_QUARTICS=1, and GENERAL_QUADRATICS=(1 for full; 0 otherwise)."];
  Exit[2]
];

loaded = CheckAbort[
  Check[UsingFrontEnd[Needs["Matchete`"]; True], $Failed],
  $Aborted
];
If[loaded =!= True,
  Print["ERROR: Matchete not available via the service front end."];
  Exit[3]
];
Get[FileNameJoin[{scriptDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{scriptDir, "model", "LagrangianBuilder.wl"}]];
model = T3ModelFromClass[cls, alpha];
If[model === $Failed, Print["ERROR: model construction failed."]; Exit[4]];

built = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
If[!AssociationQ[built],
  Print["ERROR: Lagrangian builder did not return an association."];
  Exit[5]
];

(* Structural validation of the evaluated free Lagrangian. *)
massIndexLists = DeleteDuplicates @ Cases[
  built["LFree"],
  Coupling[name_, inds_List, ___] /; name === MF :> inds,
  Infinity
];
indexedMassQ = Length[massIndexLists] > 0 &&
  AllTrue[massIndexLists, MatchQ[#, {Index[_, NFlavor]}] &];

actualYukawas = Lookup[built, "AdditionalYukawasIncluded", -1];
validatedYukawas = Lookup[built, "AdditionalYukawasIndividuallyValidated", -1];
actualQuartics = Lookup[built, "AdditionalQuarticsIncluded", -1];
actualBilinears = Lookup[built, "AdditionalQuadraticsIncluded", -1];
blockedBilinears = Lookup[built,
  "AdditionalQuadraticsHeavyMassBasisBlocked", -1];

inventoryOK = indexedMassQ &&
  Lookup[built, "AdditionalYukawasExpected", -1] === expectedYukawas &&
  actualYukawas === expectedYukawas &&
  validatedYukawas === expectedYukawas &&
  TrueQ[Lookup[built, "FermionYukawaInventoryAssembled", False]] &&
  Lookup[built, "AdditionalQuarticsExpected", -1] === expectedQuartics &&
  actualQuartics === expectedQuartics &&
  Lookup[built, "AdditionalQuadraticsExpected", -1] === expectedBilinears &&
  actualBilinears === expectedBilinears;

plainSuccessQ = built["Status"] === "Success" &&
  TrueQ[built["FullValidation"]] && blockedBilinears === 0;
expectedBlockQ = expectedMixedMassBlock &&
  built["Status"] === "RequiresMixedScalarMassDiagonalization" &&
  !TrueQ[built["FullValidation"]] &&
  blockedBilinears === 1 && expectedBilinears === 1;

status = Which[
  !inventoryOK, "OperatorInventoryOrFermionMassFailed",
  expectedMixedMassBlock && expectedBlockQ,
    "ExpectedHeavyScalarMassBasisBlock",
  expectedMixedMassBlock, "UnexpectedBCMixedMassBehaviour",
  plainSuccessQ, "PassedCombinedUVCheckLagrangian",
  True, "CombinedUVCheckFailed"
];

If[!DirectoryQ[DirectoryName[out]],
  CreateDirectory[DirectoryName[out], CreateIntermediateDirectories -> True]
];
freeFile = out <> ".free.inputform.txt";
uvFile = out <> ".uv.inputform.txt";
Export[freeFile, ToString[built["LFree"], InputForm], "Text"];
Export[uvFile, ToString[built["LUV"], InputForm], "Text"];
report = <|
  "format" -> "t3-workstream-a-combined-uv-audit-v1",
  "model" -> key,
  "mode" -> mode,
  "status" -> status,
  "fermion_type" -> If[
    model["Fermion", "Y"] === 0 && OddQ[model["Fermion", "SU2"]],
    "Majorana", "Dirac"],
  "flavour_multiplicity" -> model["Fermion", "Multiplicity"],
  "flavour_indexed_mass_found" -> indexedMassQ,
  "mass_spec_inputform" -> ToString[T3MassSpec["F", MF], InputForm],
  "mass_index_lists_inputform" -> (ToString[#, InputForm] & /@ massIndexLists),
  "extra_yukawas_expected" -> expectedYukawas,
  "extra_yukawas_included" -> actualYukawas,
  "extra_yukawas_validated" -> validatedYukawas,
  "extra_quartics_expected" -> expectedQuartics,
  "extra_quartics_included" -> actualQuartics,
  "extra_bilinears_expected" -> expectedBilinears,
  "extra_bilinears_included" -> actualBilinears,
  "extra_bilinears_heavy_mass_basis_blocked" -> blockedBilinears,
  "operator_inventory_matches" -> inventoryOK,
  "builder_status" -> built["Status"],
  "full_matchete_check_passed" -> TrueQ[built["FullValidation"]],
  "expected_mixed_mass_basis_block" -> expectedMixedMassBlock,
  "observed_expected_mass_basis_block" -> expectedBlockQ,
  "heavy_scalar_matching_mass_basis_usable" -> (plainSuccessQ && inventoryOK),
  "allowed_interactions" -> built["AllowedInteractions"],
  "rejected_interactions" -> built["RejectedInteractions"],
  "free_lagrangian_inputform_file" -> freeFile,
  "uv_lagrangian_inputform_file" -> uvFile,
  "operator_independence_normalization_certified" -> False,
  "explicit_Z2_registered_in_matchete" -> False,
  "matching_performed" -> False,
  "UV_physical_completeness_certified" -> False
|>;
Export[out, report, "RawJSON"];
Print[key, "/", mode, ": ", status,
  "; quartics ", actualQuartics, "/", expectedQuartics,
  "; bilinears ", actualBilinears, "/", expectedBilinears,
  "; indexed mass ", indexedMassQ,
  "; builder status ", built["Status"]];
If[!MemberQ[{"PassedCombinedUVCheckLagrangian",
  "ExpectedHeavyScalarMassBasisBlock"}, status], Exit[6]];

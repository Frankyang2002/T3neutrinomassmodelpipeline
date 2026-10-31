(* Workstream A: Matchete-backed extra-fermion Yukawa smoke test.
   Usage from project root:
     wolframscript -file Lagrangian/RunGeneralFermionAudit.wl OUTPUT.json CLASS alpha
   T3_EXPERIMENTAL_GENERAL_FERMIONS=1 must be set in the environment.
   The script never runs Match, RGBeta or numerical code.
   A successful check is NOT a full-flavour mass or symmetry certificate. *)
ClearAll["Global`*"];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 3,
  Print["Usage: OUTPUT.json CLASS alpha (e.g. A m2)"]; Exit[2]];
out = ExpandFileName[args[[1]]];
cls = ToUpperCase[args[[2]]];
T3ParseAuditAlpha[s_String] := Which[
  StringMatchQ[s, "m" ~~ DigitCharacter ..], -ToExpression[StringDrop[s, 1]],
  StringMatchQ[s, "p" ~~ DigitCharacter ..], ToExpression[StringDrop[s, 1]],
  StringMatchQ[s, ("-" | "") ~~ DigitCharacter ..], ToExpression[s],
  True, $Failed
];
alpha = T3ParseAuditAlpha[args[[3]]];
If[alpha === $Failed ||
   Environment["T3_EXPERIMENTAL_GENERAL_FERMIONS"] =!= "1",
  Print["Extra-fermion audit mode not enabled or invalid alpha."]; Exit[2]];
loaded = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"] ; True],
  $Failed], $Aborted];
If[loaded =!= True, Print["Matchete unavailable"]; Exit[3]];
Get[FileNameJoin[{scriptDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{scriptDir, "model", "LagrangianBuilder.wl"}]];
model = T3ModelFromClass[cls, alpha];
If[model === $Failed, Print["Unsupported T3 model"]; Exit[2]];
result = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
If[!AssociationQ[result], Print["Lagrangian build failed"]; Exit[5]];
expected = Lookup[result, "AdditionalYukawasExpected", -1];
included = Lookup[result, "AdditionalYukawasIncluded", -1];
validated = Lookup[result, "AdditionalYukawasIndividuallyValidated", -1];
report = <|
  "format" -> "t3-extra-fermion-matchete-check-v1",
  "model" -> cls <> "_" <> If[alpha < 0, "-", "+"] <> ToString[Abs[alpha]],
  "status" -> result["Status"],
  "extra_yukawa_expected" -> expected,
  "extra_yukawa_included" -> included,
  "extra_yukawa_individually_validated" -> validated,
  "fermion_yukawa_inventory_assembled" ->
    TrueQ[result["FermionYukawaInventoryAssembled"]],
  "full_matchete_check_passed" -> TrueQ[result["FullValidation"]],
  "allowed_interactions" -> result["AllowedInteractions"],
  "rejected_interactions" -> result["RejectedInteractions"],
  "fermion_flavor_mass_generality_certified" -> False,
  "discrete_z2_registered_in_matchete" -> False,
  "physical_symmetry_completeness_certified" -> False,
  "matching_performed" -> False
|>;
If[!DirectoryQ[DirectoryName[out]],
  CreateDirectory[DirectoryName[out], CreateIntermediateDirectories -> True]];
interactionFile = out <> ".interaction.inputform.txt";
uvFile = out <> ".uv.inputform.txt";
Export[interactionFile, ToString[result["LInt"], InputForm], "Text"];
Export[uvFile, ToString[result["LUV"], InputForm], "Text"];
report = Association[report,
  "interaction_inputform_file" -> interactionFile,
  "uv_lagrangian_inputform_file" -> uvFile
];
Export[out, report, "RawJSON"];
Print[report["model"], ": ", report["status"],
  ", extra Yukawas ", included, "/", expected, " -> ", out];
If[!TrueQ[report["fermion_yukawa_inventory_assembled"]] ||
   validated =!= expected || included =!= expected ||
   report["status"] =!= "Success" ||
   !TrueQ[report["full_matchete_check_passed"]], Exit[6]];

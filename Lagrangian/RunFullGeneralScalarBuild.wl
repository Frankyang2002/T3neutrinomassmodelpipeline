(* Symbolic UV Lagrangian build with all permitted scalar quartics and masses.
   Usage: wolframscript -file Lagrangian/RunFullGeneralScalarBuild.wl OUT.json A m2
   Requires environment T3_EXPERIMENTAL_GENERAL_QUADRATICS=1 (set by Python
   driver).  Writes the assembled LInt as full InputForm for later inspection.
   Does NOT call Match or claim a scalar mass-eigenstate calculation. *)
ClearAll["Global`*"];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 3,
  Print["Usage: OUTPUT.json CLASS alpha (e.g. A m2)"]; Exit[2]];
out = ExpandFileName[args[[1]]];
cls = ToUpperCase[args[[2]]];
T3ParseGeneralAlpha[s_String] := Which[
  StringMatchQ[s, "m" ~~ DigitCharacter ..], -ToExpression[StringDrop[s, 1]],
  StringMatchQ[s, "p" ~~ DigitCharacter ..], ToExpression[StringDrop[s, 1]],
  StringMatchQ[s, ("-" | "") ~~ DigitCharacter ..], ToExpression[s],
  True, $Failed
];
alpha = T3ParseGeneralAlpha[args[[3]]];
If[alpha === $Failed ||
   Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] =!= "1",
  Print["Quadratic mode not enabled or invalid alpha"]; Exit[2]];
loaded = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"] ; True],
  $Failed], $Aborted];
If[loaded =!= True, Print["Matchete unavailable"]; Exit[3]];
Get[FileNameJoin[{scriptDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{scriptDir, "model", "LagrangianBuilder.wl"}]];
model = T3ModelFromClass[cls, alpha];
If[model === $Failed, Print["Unsupported T3 model"]; Exit[2]];
result = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
If[!AssociationQ[result], Print["Lagrangian build failed"]; Exit[5]];
expected = Length[T3ExtendedExtraQuadraticStructures[model]];
count = Lookup[result, "AdditionalQuadraticsIncluded", -1];
If[!DirectoryQ[DirectoryName[out]],
  CreateDirectory[DirectoryName[out], CreateIntermediateDirectories -> True]];
expressionFile = out <> ".interaction.inputform.txt";
uvFile = out <> ".uv.inputform.txt";
Export[expressionFile, ToString[result["LInt"], InputForm], "Text"];
Export[uvFile, ToString[result["LUV"], InputForm], "Text"];
report = <|
  "format" -> "t3-full-general-scalar-lagrangian-candidate-v1",
  "model" -> cls <> "_" <> If[alpha < 0, "-", "+"] <> ToString[Abs[alpha]],
  "status" -> result["Status"],
  "quadratic_expected" -> expected,
  "quadratic_included" -> count,
  "quadratic_individually_validated" -> result["AdditionalQuadraticsIndividuallyValidated"],
  "quadratic_mass_basis_blocked" -> result["AdditionalQuadraticsHeavyMassBasisBlocked"],
  "matching_mass_basis_ready" -> result["MatchingMassBasisReady"],
  "full_matchete_check_passed" -> TrueQ[result["FullValidation"]],
  "interactions" -> result["AllowedInteractions"],
  "rejected_interactions" -> result["RejectedInteractions"],
  "interaction_inputform_file" -> expressionFile,
  "uv_lagrangian_inputform_file" -> uvFile,
  "physical_matchete_rgbeta_mapping_certified" -> False,
  "scalar_mass_diagonalization_performed" -> False,
  "matching_performed" -> False
|>;
Export[out, report, "RawJSON"];
If[count =!= expected, Print["Missing symmetry-allowed mass operator"]; Exit[6]];
If[!MemberQ[{"Success", "RequiresMixedScalarMassDiagonalization"}, result["Status"]],
  Print["General scalar Lagrangian invalid: ", result["Status"]]; Exit[7]];
Print[report["model"], ": ", report["status"],
  ", bilinears included=", count, " -> ", out];

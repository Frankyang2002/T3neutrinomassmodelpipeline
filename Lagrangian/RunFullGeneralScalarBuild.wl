(* Build the general, symmetry-allowed ordinary-T3 BSM Lagrangian.
   This is an operator-inventory/Matchete-validation export, NOT matching.
   The invoking driver enables both general scalar quadratics and extra
   fermionic Yukawas.  A mixed scalar mass basis may block matching. *)
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
   Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] =!= "1" ||
   Environment["T3_EXPERIMENTAL_GENERAL_FERMIONS"] =!= "1",
  Print["Complete scalar + fermion symbolic mode not enabled or invalid alpha"];
  Exit[2]];
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
yExpected = Lookup[result, "AdditionalYukawasExpected", -1];
yIncluded = Lookup[result, "AdditionalYukawasIncluded", -1];
qExpected = Lookup[result, "AdditionalQuarticsExpected", -1];
qIncluded = Lookup[result, "AdditionalQuarticsIncluded", -1];
If[!DirectoryQ[DirectoryName[out]],
  CreateDirectory[DirectoryName[out], CreateIntermediateDirectories -> True]];
expressionFile = out <> ".interaction.inputform.txt";
uvFile = out <> ".uv.inputform.txt";
bsmFile = out <> ".bsm.inputform.txt";
Export[expressionFile, ToString[result["LInt"], InputForm], "Text"];
Export[uvFile, ToString[result["LUV"], InputForm], "Text"];
Export[bsmFile, ToString[result["LBSM"], InputForm], "Text"];
report = <|
  "format" -> "t3-full-general-scalar-lagrangian-candidate-v1",
  "model" -> cls <> "_" <> If[alpha < 0, "-", "+"] <> ToString[Abs[alpha]],
  "status" -> result["Status"],
  "quadratic_expected" -> expected,
  "quadratic_included" -> count,
  "quadratic_individually_validated" -> result["AdditionalQuadraticsIndividuallyValidated"],
  "quadratic_mass_basis_blocked" -> result["AdditionalQuadraticsHeavyMassBasisBlocked"],
  "quartic_expected" -> qExpected,
  "quartic_included" -> qIncluded,
  "extra_yukawa_expected" -> yExpected,
  "extra_yukawa_included" -> yIncluded,
  "extra_yukawas_requested" -> TrueQ[result["AdditionalYukawasRequested"]],
  "general_scalar_quadratics_requested" -> TrueQ[result["AdditionalQuadraticsRequested"]],
  "matching_mass_basis_ready" -> result["MatchingMassBasisReady"],
  "full_matchete_check_passed" -> TrueQ[result["FullValidation"]],
  "interactions" -> result["AllowedInteractions"],
  "rejected_interactions" -> result["RejectedInteractions"],
  "interaction_inputform_file" -> expressionFile,
  "uv_lagrangian_inputform_file" -> uvFile,
  "bsm_lagrangian_inputform_file" -> bsmFile,
  "physical_matchete_rgbeta_mapping_certified" -> False,
  "scalar_mass_diagonalization_performed" -> False,
  "matching_performed" -> False
|>;
Export[out, report, "RawJSON"];
If[count =!= expected || qIncluded =!= qExpected ||
    yIncluded =!= yExpected || yExpected < 0 || qExpected < 0,
  Print["Incomplete general scalar/Yukawa interaction inventory"]; Exit[6]];
If[!MemberQ[{"Success", "RequiresMixedScalarMassDiagonalization"}, result["Status"]],
  Print["General scalar+fermion Lagrangian invalid: ", result["Status"]]; Exit[7]];
Print[report["model"], ": ", report["status"],
  ", extra quartics=", qIncluded,
  ", extra bilinears=", count,
  ", extra Yukawas=", yIncluded,
  " -> ", out];

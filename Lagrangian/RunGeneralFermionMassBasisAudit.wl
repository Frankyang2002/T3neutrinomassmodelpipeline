(* Workstream A: opt-in, three-flavour diagonal fermion mass validation.
   Usage from the project root:
     wolframscript -file Lagrangian/RunGeneralFermionMassBasisAudit.wl \
       OUTPUT.json CLASS alpha
   Required environment:
     T3_EXPERIMENTAL_FLAVOR_MASSES=1
     T3_EXPERIMENTAL_GENERAL_FERMIONS=1
     T3_EXPERIMENTAL_GENERAL_QUARTICS=0
     T3_EXPERIMENTAL_GENERAL_QUADRATICS=0
   Uses Matchete CheckLagrangian; never invokes Match, RGBeta or numerics.
   A diagonal 3-flavour mass basis is not a UV operator-completeness proof. *)

ClearAll["Global`*"];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 3,
  Print["Usage: OUTPUT.json CLASS alpha (example: A m2)"];
  Exit[2]
];

out = ExpandFileName[args[[1]]];
cls = ToUpperCase[args[[2]]];
T3ParseMassAuditAlpha[s_String] := Which[
  StringMatchQ[s, "m" ~~ DigitCharacter ..], -ToExpression[StringDrop[s, 1]],
  StringMatchQ[s, "p" ~~ DigitCharacter ..], ToExpression[StringDrop[s, 1]],
  StringMatchQ[s, ("-" | "") ~~ DigitCharacter ..], ToExpression[s],
  True, $Failed
];
alpha = T3ParseMassAuditAlpha[args[[3]]];
If[alpha === $Failed,
  Print["ERROR: invalid alpha label: ", args[[3]]]; Exit[2]
];

If[Environment["T3_EXPERIMENTAL_FLAVOR_MASSES"] =!= "1" ||
   Environment["T3_EXPERIMENTAL_GENERAL_FERMIONS"] =!= "1",
  Print["ERROR: both experimental flavour-mass and general-fermion flags must be 1."];
  Exit[2]
];
If[Environment["T3_EXPERIMENTAL_GENERAL_QUARTICS"] === "1" ||
   Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] === "1",
  Print["ERROR: disable general scalar extensions for this isolated mass audit."];
  Exit[2]
];

loaded = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"] ; True],
  $Failed], $Aborted];
If[loaded =!= True,
  Print["ERROR: Matchete unavailable through the service front end."];
  Exit[3]
];

Get[FileNameJoin[{scriptDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{scriptDir, "model", "LagrangianBuilder.wl"}]];
model = T3ModelFromClass[cls, alpha];
If[model === $Failed,
  Print["ERROR: unsupported T3 model."];
  Exit[2]
];

result = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
If[!AssociationQ[result],
  Print["ERROR: Lagrangian build did not return an association: ",
    ToString[result, InputForm]];
  Exit[5]
];

(* Probe the evaluated FreeLag, not the field-option syntax alone.  In the
   successful standalone probe the mass term is:
     Coupling[MF, {Index[d$$1, NFlavor]}, 0]
   There must be no unindexed Coupling[MF, {}, 0] remaining. *)
massIndexLists = DeleteDuplicates @ Cases[
  result["LFree"],
  Coupling[name_, indices_List, ___] /; name === MF :> indices,
  Infinity
];
indexedMassQ = Length[massIndexLists] > 0 &&
  AllTrue[massIndexLists, MatchQ[#, {Index[_, NFlavor]}] &];

extraExpected = result["AdditionalYukawasExpected"];
extraIncluded = result["AdditionalYukawasIncluded"];
extraValidated = result["AdditionalYukawasIndividuallyValidated"];
extraOK = TrueQ[result["FermionYukawaInventoryAssembled"]] &&
  extraIncluded === extraExpected && extraValidated === extraExpected;
uvOK = result["Status"] === "Success" &&
  TrueQ[result["FullValidation"]];

status = Which[
  !indexedMassQ, "MissingOrNonDiagonalFlavorMassIndex",
  !extraOK, "GeneralFermionYukawaValidationFailed",
  !uvOK, "UVCheckLagrangianFailed",
  True, "PassedDiagonalFermionMassBasis"
];

If[!DirectoryQ[DirectoryName[out]],
  CreateDirectory[DirectoryName[out], CreateIntermediateDirectories -> True]
];
freeFile = out <> ".free.inputform.txt";
uvFile = out <> ".uv.inputform.txt";
Export[freeFile, ToString[result["LFree"], InputForm], "Text"];
Export[uvFile, ToString[result["LUV"], InputForm], "Text"];

report = <|
  "format" -> "t3-flavour-diagonal-mass-matchete-audit-v1",
  "model" -> cls <> "_" <> If[alpha < 0, "-", "+"] <>
    ToString[Abs[alpha]],
  "status" -> status,
  "fermion_type" -> If[
    model["Fermion", "Y"] === 0 && OddQ[model["Fermion", "SU2"]],
    "Majorana", "Dirac"
  ],
  "fermion_flavour_multiplicity" -> model["Fermion", "Multiplicity"],
  "mass_spec_requested_inputform" ->
    ToString[T3MassSpec["F", MF], InputForm],
  "flavour_mass_index_detected" -> indexedMassQ,
  "fermion_mass_coupling_instances" -> Length[massIndexLists],
  "fermion_mass_index_lists_inputform" ->
    (ToString[#, InputForm] & /@ massIndexLists),
  "general_fermion_yukawas_expected" -> extraExpected,
  "general_fermion_yukawas_included" -> extraIncluded,
  "general_fermion_yukawas_validated" -> extraValidated,
  "full_matchete_check_passed" -> TrueQ[result["FullValidation"]],
  "lagrangian_status" -> result["Status"],
  "free_lagrangian_inputform_file" -> freeFile,
  "uv_lagrangian_inputform_file" -> uvFile,
  "general_off_diagonal_mass_matrix_tested" -> False,
  "matching_performed" -> False,
  "scalar_extensions_tested" -> False,
  "physical_symmetry_completeness_certified" -> False
|>;
Export[out, report, "RawJSON"];
Print[report["model"], ": ", status,
  "; MF index lists=", report["fermion_mass_index_lists_inputform"],
  "; Matchete UV status=", result["Status"]];
If[status =!= "PassedDiagonalFermionMassBasis", Exit[6]];

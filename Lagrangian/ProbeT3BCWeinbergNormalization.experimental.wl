(* Optional installed-model probe. The Python v9.1 source-only fallback does
   not require this Wolfram script. Never run Match or modify production code. *)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[!MemberQ[{1, 2}, Length[args]],
  Print["Usage: wolframscript -file Lagrangian/ProbeT3BCWeinbergNormalization.experimental.wl OUTPUT.json [WARSAW_MODEL_PATH]"];
  Exit[2]];
outPath = ExpandFileName[args[[1]]];
If[!DirectoryQ[DirectoryName[outPath]],
  CreateDirectory[DirectoryName[outPath], CreateIntermediateDirectories -> True]];
format = "t3-bc-matchete-warsaw-dim5-probe-v1";
reportFailure[stage_, code_, models_] := (
  Export[outPath, <|
    "format" -> format, "status" -> stage,
    "Mathematica_version" -> ToString[$Version],
    "available_models" -> models,
    "requested_model_path" -> If[Length[args] == 2, args[[2]], Null],
    "physical_C5" -> Null,
    "absolute_Matchete_to_Warsaw_factor_certified" -> False
  |>, "RawJSON"];
  Print["T3_WARSAW: ", stage, "; report=", outPath];
  Exit[code]
);
Print["T3_WARSAW_STAGE: matchete-load-start"];
loaded = CheckAbort[Quiet@Check[UsingFrontEnd[Needs["Matchete`"]]; True, $Failed], $Aborted];
If[loaded =!= True, reportFailure["MatcheteLoadFailed", 3, {}]];
modelsRaw = Quiet@Check[GetModels[], {}];
modelNames = If[ListQ[modelsRaw], ToString[#, InputForm] & /@ modelsRaw, {}];
Print["T3_WARSAW: available models: ", modelNames];
modelTarget = If[Length[args] == 2, ExpandFileName[args[[2]]], "SMEFT_Warsaw"];
If[Length[args] == 2 && !FileExistsQ[modelTarget],
  reportFailure["WarsawModelUnavailable", 4, modelNames]];
Print["T3_WARSAW_STAGE: load-reference-model ", modelTarget];
lag = CheckAbort[Quiet@Check[LoadModel[modelTarget], $Failed], $Aborted];
If[lag === $Failed || lag === $Aborted ||
   StringContainsQ[ToString[lag, InputForm], "LoadModel["],
  reportFailure["WarsawModelUnavailable", 4, modelNames]];
all = Expand[lag];
allTerms = If[Head[all] === Plus, List@@all, {all}];
strings = ToString[#, InputForm] & /@ allTerms;
selectQ[s_String] :=
  StringCount[s, "Field[H, Scalar"] == 2 &&
  StringCount[s, "Field[l, Fermion"] == 2 &&
  StringContainsQ[s, "GammaCC"];
weinbergTerms = Select[strings, selectQ];
termData = Map[Function[s, <|
    "inputform" -> s,
    "H_fields" -> StringCount[s, "Field[H, Scalar"],
    "barred_H_fields" -> StringCount[s, "Bar[Field[H, Scalar"],
    "lepton_fields" -> StringCount[s, "Field[l, Fermion"],
    "barred_lepton_fields" -> StringCount[s, "Bar[Field[l, Fermion"],
    "has_GammaCC" -> StringContainsQ[s, "GammaCC"],
    "has_PL" -> StringContainsQ[s, "Proj[-1]"],
    "has_PR" -> StringContainsQ[s, "Proj[1]"],
    "coupling_labels" -> DeleteDuplicates[StringCases[s,
      RegularExpression["Coupling\\[([A-Za-z0-9$]+),"] -> "$1"]]
  |>], weinbergTerms];
status = If[Length[termData] > 0, "WarsawDim5StructureExported",
            "WarsawDim5NotRecognized"];
Export[outPath, <|
  "format" -> format, "status" -> status,
  "Mathematica_version" -> ToString[$Version],
  "reference_model" -> modelTarget,
  "available_models" -> modelNames,
  "model_term_count" -> Length[allTerms],
  "weinberg_source_terms_found" -> Length[termData],
  "weinberg_terms" -> termData,
  "definition_from_public_Matchete_reference" ->
    "Q_llHH^pr = eps_in eps_jm H^i H^j [(ell_p^n)^T C ell_r^m]",
  "absolute_Matchete_to_Warsaw_factor_certified" -> False,
  "physical_C5" -> Null, "production_code_modified" -> False
|>, "RawJSON"];
Print["T3_WARSAW: status=", status, "; terms=", Length[termData]];
If[Length[termData] == 0, Exit[5]];
Exit[0];

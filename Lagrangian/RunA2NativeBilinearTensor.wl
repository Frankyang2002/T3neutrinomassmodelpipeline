(* A(-2): extract the ACTUAL Matchete quadratic invariant tensor.
   No guessed Clebsch-Gordan normalization, no Match[] or RGE changes.
   Usage: wolframscript -file Lagrangian/RunA2NativeBilinearTensor.wl OUT.json
   Requires the opt-in full-general-scalar builder installed in this project.
*)
ClearAll["Global`*"];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 1, Print["Usage: OUT.json"]; Exit[2]];
outPath = ExpandFileName[First[args]];
If[Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] =!= "1" ||
   Environment["T3_EXPERIMENTAL_GENERAL_QUARTICS"] =!= "1",
  Print["Full scalar candidate flags not set"]; Exit[2]];
loaded = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"] ; True], $Failed], $Aborted];
If[loaded =!= True, Print["Matchete unavailable"]; Exit[3]];
Get[FileNameJoin[{scriptDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{scriptDir, "model", "LagrangianBuilder.wl"}]];
model = T3ModelFromClass["A", -2];
If[model === $Failed, Exit[4]];
report = CheckAbort[UsingFrontEnd[Module[
    {built, structures, counts, names, basis, metric, k, residual,
     candidateNames, metricComponents, result},
    built = BuildT3Lagrangian[model];
    If[!AssociationQ[built] || built["Status"] =!= "Success", Return[$Failed]];
    structures = T3ExtendedExtraQuadraticStructures[model];
    If[Length[structures] =!= 1 ||
       Lookup[built, "AdditionalQuadraticsIncluded", -1] =!= 1,
      Return[$Failed]];
    counts = First[structures];
    names = T3ExtendedSpeciesList[counts];
    If[names =!= {"S2dag", "S2dag"}, Return[$Failed]];
    candidateNames = Select[Lookup[built, "AllowedInteractions", {}],
      StringStartsQ[#, "T3ExtraM_"] &];
    If[Length[candidateNames] =!= 1, Return[$Failed]];
    basis = T3ExtendedBosonicTensorBasis[model, names];
    If[!AssociationQ[basis] || basis["PhysicalCount"] =!= 1,
      Return[$Failed]];
    metric = Normal[First[basis["Basis"]]];
    If[Dimensions[metric] =!= {3,3}, Return[$Failed]];
    k = FullSimplify[Tr[metric]/3];
    If[TrueQ[PossibleZeroQ[k]], Return[$Failed]];
    residual = FullSimplify[metric - k IdentityMatrix[3]];
    If[!TrueQ[And @@ (TrueQ[PossibleZeroQ[#]] & /@ Flatten[residual])],
      Print["Matchete triplet quadratic tensor is NOT proportional to identity: ",
        ToString[metric, InputForm]];
      Return[$Failed]];
    metricComponents = Map[ToString[#, InputForm] &, metric, {2}];
    result = <|
      "format" -> "t3-a-minus2-native-bilinear-tensor-v1",
      "model" -> "A_-2", "status" -> "NativeMatcheteBilinearMetricExtracted",
      "species" -> {"S2dag", "S2dag"},
      "native_triplet_dimension" -> 3,
      "native_metric_inputform" -> metricComponents,
      "native_metric_factor_inputform" -> ToString[k, InputForm],
      "native_metric_is_factor_times_identity" -> True,
      "mass_operator_in_full_lagrangian" -> True,
      "mass_operator_names" -> candidateNames,
      "matchete_lagrangian_convention" -> "+ muMat * (S2dag)^T (k I3) S2dag + h.c.",
      "rgbeta_mass_tensor_checked" -> False,
      "matching_performed" -> False
    |>;
    result
  ]], $Aborted];
If[!AssociationQ[report], Print["Native mass tensor extraction failed"]; Exit[5]];
If[!DirectoryQ[DirectoryName[outPath]],
   CreateDirectory[DirectoryName[outPath], CreateIntermediateDirectories -> True]];
Export[outPath, report, "RawJSON"];
Print["A_-2 native Matchete metric coefficient = ",
  report["native_metric_factor_inputform"], " -> ", outPath];

(* Read-only native SU(2) CG evidence for A_-4 F-only exchange.
   Project path: output/general_observable/lo_fermion_boundary/native_cg_a_m4_probe.wl
   Never touches the 16-model Wilson boundary or the earlier native FR report.
   Uses same installed Matchete source and A_-4 builder as the successful
   native_feynman_rules_audit.wl, and reports raw CG phases explicitly.
*)
ClearAll["Global`*"];
probeFile = ExpandFileName[$InputFileName];
repoRoot = Nest[DirectoryName, probeFile, 4];
probeDir = DirectoryName[probeFile];
outDir = FileNameJoin[{probeDir, "normalization_probes"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];
outPath = FileNameJoin[{outDir, "native_cg_a_m4_" <>
  StringReplace[CreateUUID[], "-" -> ""] <> ".json"}];
stop[reason_String] := (Print["NATIVE_CG_BLOCKED=", reason]; Exit[2]);
Print["NATIVE_CG_STAGE=load"];
loaded = CheckAbort[UsingFrontEnd[Needs["Matchete`"]; True], $Aborted];
If[loaded =!= True, stop["Matchete frontend loading failed"]];
sourceDir = FileNameJoin[{repoRoot, "Lagrangian"}];
Get[FileNameJoin[{sourceDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{sourceDir, "model", "LagrangianBuilder.wl"}]];
If[!And @@ (Environment[#] === "1" & /@ {
     "T3_EXPERIMENTAL_GENERAL_FERMIONS",
     "T3_EXPERIMENTAL_GENERAL_QUARTICS",
     "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
     "T3_EXPERIMENTAL_FLAVOR_MASSES"}),
   stop["All four experimental general flags must equal 1"]];
model = T3ModelFromClass["A", -4, 3];
If[!AssociationQ[model] || SetT3HeavyFields[{"F"}] =!= True,
   stop["A_-4 F-first model initialization failed"]];
Print["NATIVE_CG_STAGE=BuildT3Lagrangian"];
built = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
If[!AssociationQ[built] || Lookup[built, "Status", ""] =!= "Success" ||
    Lookup[built, "FullValidation", False] =!= True,
  stop["Native full-general UV validation failed"]];
ClearAll[cgEvidence];
cgEvidence[name_Symbol] := Module[{tensor, reps, dense, dims, values},
  tensor = Quiet@Check[
    Matchete`CGManipulations`PackagePrivate`$CGtensors[name], $Failed];
  reps = Quiet@Check[Matchete`PackageScope`$CGproperties[name, Indices], $Failed];
  If[tensor === $Failed || reps === $Failed || MissingQ[tensor] || MissingQ[reps],
    Return[$Failed]];
  dense = Normal[tensor];
  dims = Dimensions[dense];
  values = Flatten[dense];
  If[!VectorQ[dims, IntegerQ] || Times @@ dims =!= Length[values],
    Return[$Failed]];
  <|
    "name" -> ToString[Unevaluated[name], InputForm],
    "representation_inputform" -> ToString[InputForm[reps]],
    "tensor_inputform" -> ToString[InputForm[tensor]],
    "dimensions" -> dims,
    "real_flat" -> (N[Re[#]] & /@ values),
    "imag_flat" -> (N[Im[#]] & /@ values)
   |>
];
Print["NATIVE_CG_STAGE=read-native-CG-registry"];
y1 = cgEvidence[T3Y1CG];
y2 = cgEvidence[T3Y2CG];
mix = cgEvidence[T3MixCG];
If[!And @@ (AssociationQ /@ {y1, y2, mix}),
  stop["One or more native CG tensors could not be read"]];
If[y1["dimensions"] =!= {2, 2} || y2["dimensions"] =!= {2, 2, 3} ||
   mix["dimensions"] =!= {2, 2, 3},
  stop["A_-4 CG dimensions inconsistent with the 2x2 and 2x2x3 vertices"]];
res = <|
  "format" -> "t3-native-a-m4-yukawa-mixing-cg-evidence-v1",
  "model" -> "A_-4",
  "status" -> "NativeCGTensorsExtracted_NotPhaseCertified",
  "native_yukawa1_cg" -> y1,
  "native_yukawa2_cg" -> y2,
  "native_matched_mix_cg" -> mix,
  "holomorphic_native_source_has_conjugated_yukawa_CGs" -> True,
  "cg_product_identity_computed" -> False,
  "master_weyl_phase_certified" -> False,
  "physical_16_model_wilson_normalization_certified" -> False,
  "previous_release_changed" -> False
|>;
Export[outPath, res, "RawJSON"];
Print["NATIVE_CG_REPORT=", outPath];
Print["NATIVE_CG_STATUS=", res["status"]];
Exit[0];

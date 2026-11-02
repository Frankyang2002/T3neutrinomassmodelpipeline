(* WS4: live Matchete restricted UV operator inventory for A(0), B(-1).
   This is an additive native witness, NOT a component Feynman-rule or
   cross-package normalisation certificate. No Match[], RGBeta or EFT runs.
   Run: wolframscript -file scripts/RunLODirectUVNativeOperatorExport.wl OUTPUT.json
   Output path must not already exist. *)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 1, Print["Usage: OUTPUT.json"]; Exit[2]];
out = ExpandFileName[First[args]];
If[FileExistsQ[out], Print["Refusing to overwrite ", out]; Exit[2]];
If[AnyTrue[{"T3_EXPERIMENTAL_GENERAL_FERMIONS", "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
             "T3_EXPERIMENTAL_GENERAL_QUARTICS", "T3_EXPERIMENTAL_FLAVOR_MASSES"},
            Environment[#] === "1" &],
   Print["BLOCKED: this witness must use the original restricted default Lagrangian"]; Exit[2]];
repo = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
source = FileNameJoin[{repo, "Lagrangian", "model"}];
loading = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]; True], $Failed], $Aborted];
If[loading =!= True, Print["Matchete unavailable"]; Exit[3]];
Get[FileNameJoin[{source,"T3ModelCatalog.wl"}]];
Get[FileNameJoin[{source,"LagrangianBuilder.wl"}]];

ClearAll[ProbeModel];
ProbeModel[class_String, alpha_Integer] := Module[
  {model,built,expr,allow,forms,scope},
  model = T3ModelFromClass[class, alpha];
  If[!AssociationQ[model], Return[$Failed]];
  built = CheckAbort[Check[UsingFrontEnd[BuildT3Lagrangian[model]], $Failed], $Aborted];
  If[!AssociationQ[built], Return[$Failed]];
  If[Lookup[built,"Status",""] =!= "Success" ||
     !TrueQ[Lookup[built,"FullValidation",False]], Return[$Failed]];
  expr = ToString[built["LInt"], InputForm];
  allow = Lookup[built,"AllowedInteractions",{}];
  forms = <|
    "unindexed_lambdaH1_function_in_LInt" -> StringContainsQ[expr,"lambdaH1["],
    "unindexed_lambdaH2_function_in_LInt" -> StringContainsQ[expr,"lambdaH2["],
    "indexed_lambdaH1Inv_in_LInt" -> StringContainsQ[expr,"lambdaH1Inv"],
    "indexed_lambdaH2Inv_in_LInt" -> StringContainsQ[expr,"lambdaH2Inv"],
    "topology_y1_in_LInt" -> StringContainsQ[expr,"y1["],
    "topology_y2_in_LInt" -> StringContainsQ[expr,"y2["],
    "T3_mixing_in_LInt" -> StringContainsQ[expr,"lambdaT3["]
  |>;
  scope = <|
    "model" -> class<>"_"<>If[alpha>=0,"+",""]<>ToString[alpha],
    "native_build_status" -> built["Status"],
    "native_full_lagrangian_validation" -> True,
    "allowed_operator_names" -> allow,
    "native_restricted_LInt_function_names" -> forms,
    "native_LInt_string_bytes" -> StringLength[expr],
    "native_full_LInt_inputform" -> expr,
    "physical_component_spinor_mass_and_portal_normalisation_certified" -> False
  |>;
  scope
];
models = {ProbeModel["A",0],ProbeModel["B",-1]};
If[!AllTrue[models,AssociationQ],
   Print["BLOCKED: native Matchete build failed for A_+0 or B_-1"];Exit[5]];
res = <|
 "format" -> "t3-ws4-native-matchete-restricted-LInt-operator-export-v1",
 "status" -> "native_restricted_operator_inventory_extracted_not_component_physical",
 "models" -> models,
 "source_files" -> <|
   "LagrangianBuilder" -> FileHash[FileNameJoin[{source,"LagrangianBuilder.wl"}], "SHA256"],
   "ScalarPotential" -> FileHash[FileNameJoin[{repo,"Lagrangian","interactions","ScalarPotential.wl"}], "SHA256"],
   "T3Topology" -> FileHash[FileNameJoin[{repo,"Lagrangian","interactions","T3Topology.wl"}], "SHA256"]|>,
 "physical_cross_package_normalisation_certified" -> False,
 "missing_native_objects" -> {"Concrete GammaCC numeric components in declared Weyl axes",
   "Both Yukawa Feynman rule coefficients with canonical kinetic terms",
   "RGBeta epsilon/triplet basis in same component axes",
   "Scalar Hessian derivative of isolated lambdaHiInvN operator"},
 "one_loop_or_EFT_matching_performed" -> False
|>;
If[!DirectoryQ[DirectoryName[out]],
   CreateDirectory[DirectoryName[out],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[out,res,"RawJSON"],$Failed]===$Failed,
   Print["BLOCKED: native report export failed"];Exit[6]];
Print["Live Matchete restricted operator inventory exported: ",out];
Print["Physical Matchete/RGBeta component normalisation: NOT CERTIFIED"];
Exit[0];

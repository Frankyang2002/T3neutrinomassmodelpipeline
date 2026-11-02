(* Workstream 4: read-only registered Matchete lambdaT3 mixing tensor for A-E.
   wolframscript -file scripts/RunLODirectUVNativeMixingTensorExport.wl OUTPUT.json
   Uses ORIGINAL restricted UV Lagrangian. This is not a full physical
   Matchete/RGBeta-to-Weyl certification or an EFT matching calculation. *)
ClearAll["Global`*"];
$MXArgs = Rest[$ScriptCommandLine];
If[Length[$MXArgs] =!= 1, Print["Usage: OUTPUT.json"]; Exit[2]];
$MXOut = ExpandFileName[First[$MXArgs]];
If[FileExistsQ[$MXOut], Print["Refusing to overwrite ", $MXOut]; Exit[2]];
$MXRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$MXFlags = {"T3_EXPERIMENTAL_GENERAL_FERMIONS",
  "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
  "T3_EXPERIMENTAL_GENERAL_QUARTICS",
  "T3_EXPERIMENTAL_FLAVOR_MASSES"};
$MXInherited = AssociationMap[
  Function[n, With[{v = Environment[n]}, If[StringQ[v], v, "UNSET"]]],
  $MXFlags];
$MXOverride = Check[SetEnvironment[Thread[$MXFlags -> ConstantArray["0",Length[$MXFlags]]]]; True,$Failed];
If[$MXOverride =!= True || !AllTrue[$MXFlags, Environment[#] === "0" &],
  Print["BLOCKED: restricted environment could not be established"]; Exit[3]];
$MXEffective = AssociationMap[Environment,$MXFlags];
$MXLoad = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"];True],$Failed],$Aborted];
If[$MXLoad =!= True, Print["BLOCKED: Matchete unavailable"]; Exit[4]];
Get[FileNameJoin[{$MXRoot,"Lagrangian","model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{$MXRoot,"Lagrangian","model","LagrangianBuilder.wl"}]];
ClearAll[MXFail, MXText, MXModel];
MXFail[msg_String] := Throw[<|"status"->"blocked","error"->msg|>,"MXAbort"];
MXText[x_] := ToString[FullSimplify[x],InputForm];
MXModel[class_String, alpha_Integer] := Module[
  {short,model,built,d1,d2,df,expr,allowed,form,raw,arr,dims,expected,flat,symHc},
  short = class<>"_"<>If[alpha>=0,"+",""]<>ToString[alpha];
  model = T3ModelFromClass[class,alpha];
  If[!AssociationQ[model],MXFail[short<>" catalog model missing"]];
  {d1,d2,df} = {model["Scalar1","SU2"],model["Scalar2","SU2"],model["Fermion","SU2"]};
  built = CheckAbort[Check[UsingFrontEnd[BuildT3Lagrangian[model]],$Failed],$Aborted];
  If[!AssociationQ[built] || Lookup[built,"Status",""] =!= "Success" ||
     !TrueQ[Lookup[built,"FullValidation",False]],
     MXFail[short<>" restricted native Lagrangian build failed"]];
  allowed = Lookup[built,"AllowedInteractions",{}];
  If[Count[allowed,"T3ScalarMix_GenericCG"] =!= 1,
    MXFail[short<>" mixing operator not selected exactly once"]];
  (* Test the ACTUAL expanded Matchete expression rather than a Module-local
     HoldPattern[CG[sym,...]], which previously produced a false negative. *)
  expr = ToString[built["LInt"],InputForm];
  If[!StringContainsQ[expr,"CG[T3MixCG,"] ||
     !StringContainsQ[expr,"Coupling[lambdaT3,"] ||
     !StringContainsQ[expr,"CG[Bar[T3MixCG],"] ||
     !StringContainsQ[expr,"Bar[Coupling[lambdaT3,"],
     MXFail[short<>" assembled LInt lacks both native mixing orientations"]];
  raw = CheckAbort[Quiet@Check[GetCGTensor[T3MixCG],$Failed],$Aborted];
  If[raw === $Failed || raw === $Aborted,
     MXFail[short<>" Matchete GetCGTensor[T3MixCG] unavailable"]];
  arr = If[Head[raw] === SparseArray,Normal[raw],raw];
  dims = Dimensions[arr];
  expected = Select[{2,2,d1,d2},#>1&];
  If[dims =!= expected,
    MXFail[short<>" native mixing tensor dimensions "<>ToString[dims,InputForm]<>
     " not "<>ToString[expected,InputForm]]];
  flat = Flatten[arr];
  If[Length[flat] =!= Times@@expected ||
     !AllTrue[flat,TrueQ[NumericQ[FullSimplify[#]]] &] ||
     !AnyTrue[flat,!TrueQ[PossibleZeroQ[#]] &],
    MXFail[short<>" native mixing tensor is incomplete or zero"]];
  If[!AllTrue[Flatten[arr-Transpose[arr,Join[{2,1},Range[3,Length[dims]]]]],
       TrueQ[PossibleZeroQ[FullSimplify[#]]] &],
    MXFail[short<>" native tensor not symmetric in identical Higgs axes"]];
  <|"class"->class,"model"->short,"alpha"->alpha,
    "dS1"->d1,"dS2"->d2,"dF"->df,
    "native_build_status"->"Success","native_full_lagrangian_validation"->True,
    "mixing_operator_name"->"T3ScalarMix_GenericCG",
    "mixing_coupling"->"lambdaT3",
    "native_cg"-><|"symbol"->"T3MixCG","tensor_accessor"->"Matchete`GetCGTensor",
      "dimensions"->dims,"full_dimensions_including_singlets"->{2,2,d1,d2},
      "components_inputform"->(MXText/@flat),
      "components_real"->(N[Re[#],MachinePrecision]&/@flat),
      "components_imag"->(N[Im[#],MachinePrecision]&/@flat),
      "registered_in_assembled_interaction"->True|>,
    "native_nonhc_operator_present"->True,
    "native_hc_operator_present"->True,
    "Higgs_fields_in_nonhc_operator"->"H H",
    "scalar_fields_in_nonhc_operator"->"S1 Bar[S2]",
    "native_nonhc_LInt_coefficient"->"+lambdaT3 * T3MixCG (no 1/2)",
    "physical_source_to_Weyl_certified"->False|>
];
$MXSources = {"Lagrangian/model/T3ModelCatalog.wl",
  "Lagrangian/model/LagrangianBuilder.wl",
  "Lagrangian/model/SU2Invariants.wl",
  "Lagrangian/model/T3Fields.wl",
  "Lagrangian/interactions/T3Topology.wl",
  "scripts/BuildT3DirectOneLoopInputs.py"};
$MXDigest = Association@Table[With[{f=FileNameJoin[{$MXRoot,path}]},
  If[!FileExistsQ[f],MXFail["missing native source "<>path]];
  path->IntegerString[FileHash[f,"SHA256"],16,64]],{path,$MXSources}];
$MXCases = {{"A",0},{"B",-1},{"C",-1},{"D",0},{"E",0}};
$MXResult = CheckAbort[Catch[Table[MXModel@@spec,{spec,$MXCases}],"MXAbort"],$Aborted];
If[!ListQ[$MXResult] || Length[$MXResult] =!= 5 ||
    !AllTrue[$MXResult,AssociationQ],
  Print["BLOCKED: ",If[AssociationQ[$MXResult],Lookup[$MXResult,"error","unknown"],
    "one or more Matchete native builds aborted"]]; Exit[5]];
$MXReport = <|
  "format"->"t3-ws4-native-mixing-tensor-five-classes-v1",
  "status"->"five_classes_native_registered_mixing_tensors_extracted",
  "native_source_sha256"->$MXDigest,
  "restricted_environment"-><|"inherited"->$MXInherited,"effective"->$MXEffective,
    "isolated_to_wolframscript_process"->True|>,
  "models"->$MXResult,
  "component_axes"->"Original registered Matchete index axes, not generator weight basis",
  "H0_field_orientation_source"->"H[2] in generated weight basis (subject to generator rotation)",
  "vev_declaration"->"H0=v/Sqrt[2]; LInt has +lambdaT3 H H S1 Bar[S2]+HC",
  "native_rgbeta_mixing_coupling_dictionary_certified"->False,
  "native_spinor_and_mass_dictionary_certified"->False,
  "matching_or_running_performed"->False|>;
If[!DirectoryQ[DirectoryName[$MXOut]],
  CreateDirectory[DirectoryName[$MXOut],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[$MXOut,$MXReport,"RawJSON"],$Failed] === $Failed,
  Print["BLOCKED: native JSON export failed"];Exit[6]];
Print["Registered Matchete lambdaT3 mixing tensors: 5/5 classes"];
Print["Higgs-axis symmetry and both LInt orientations: 5/5"];
Print["Physical Matchete/RGBeta-to-Weyl normalisation: NOT CERTIFIED"];
Print["Output: ",$MXOut];
Exit[0];

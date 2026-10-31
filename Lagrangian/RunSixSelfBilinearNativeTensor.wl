(* Six ordinary A/D/E self-bilinear models: independent native Matchete
   operator and CG extraction, including the gauge-singlet case.
   Usage: wolframscript -file Lagrangian/RunSixSelfBilinearNativeTensor.wl OUTPUT.json A 0-2
   This never invokes Match[] and never modifies production output.

   WolframScript may consume a literal negative positional argument ("-2") as
   a command-line option before constructing $ScriptCommandLine. The existing
   Python batch passes -2 literally. Recover that value only when the complete
   command line ends with -2, or when the two surviving arguments identify the
   batch's exact <family>_m2 output directory. This preserves the original
   Python CLI, output format, native extraction and physical calculations.
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] === 2 &&
   MemberQ[{"A", "D", "E"}, ToUpperCase[args[[2]]]] &&
   (TrueQ[Length[$CommandLine] > 0 && Last[$CommandLine] === "-2"] ||
    FileNameTake[DirectoryName[ExpandFileName[args[[1]]]]] ===
      ToUpperCase[args[[2]]] <> "_m2"),
  args = Append[args, "-2"];
  Print["Recovered alpha=-2 positional argument for six-self native tensor batch"];
];
If[Length[args] =!= 3, Print["Usage: OUTPUT.json A|D|E -2|0"]; Exit[2]];
out = ExpandFileName[args[[1]]];
class = ToUpperCase[args[[2]]];
alpha = Quiet @ Check[ToExpression[args[[3]]], $Failed];
If[!MemberQ[{"A", "D", "E"}, class] || !MemberQ[{-2, 0}, alpha],
  Print["Only A/D/E alpha=-2,0 are supported"]; Exit[2]];
If[Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] =!= "1" ||
   Environment["T3_EXPERIMENTAL_GENERAL_QUARTICS"] =!= "1",
   Print["General scalar interaction flags must both be 1"]; Exit[3]];
load = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]; True], $Failed], $Aborted];
If[load =!= True, Print["Matchete unavailable"]; Exit[3]];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];
Get[FileNameJoin[{scriptDir,"model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{scriptDir,"model","LagrangianBuilder.wl"}]];
model = T3ModelFromClass[class,alpha];
If[model === $Failed, Print["Invalid T3 model"]; Exit[4]];
res = CheckAbort[UsingFrontEnd[Module[
  {build, struct, names, data, tensor, dim, factor, residual, operators, identity,
   result, species},
  build = BuildT3Lagrangian[model];
  If[!AssociationQ[build] || build["Status"] =!= "Success" ||
     Lookup[build,"AdditionalQuadraticsIncluded",-1] =!= 1, Return[$Failed]];
  struct = T3ExtendedExtraQuadraticStructures[model];
  If[Length[struct] =!= 1, Return[$Failed]];
  names = T3ExtendedSpeciesList[First[struct]];
  species = If[alpha === -2, "S2dag", "S1dag"];
  If[names =!= {species,species}, Return[$Failed]];
  dim = model[If[alpha === -2,"Scalar2","Scalar1"],"SU2"];
  If[!MemberQ[{1,3},dim], Return[$Failed]];
  data = T3ExtendedBosonicTensorBasis[model,names];
  If[!AssociationQ[data] || data["PhysicalCount"] =!= 1, Return[$Failed]];
  tensor = Normal[First[data["Basis"]]];
  If[dim === 1,
    (* A(0) and D(-2) have a trivial SU(2) singlet with no tensor index. *)
    factor = If[Head[tensor] === List, First[Flatten[tensor]], tensor];
    residual = 0;
    identity = {{factor}},
    If[Dimensions[tensor] =!= {3,3}, Return[$Failed]];
    factor = FullSimplify[Tr[tensor]/3];
    residual = FullSimplify[tensor-factor IdentityMatrix[3]];
    identity = tensor
  ];
  If[TrueQ[PossibleZeroQ[factor]] ||
     !And @@ (TrueQ[PossibleZeroQ[#]] & /@ Flatten[{residual}]),
     Print["Native self metric not proportional to identity: ",InputForm[tensor]];
     Return[$Failed]];
  operators = Select[Lookup[build,"AllowedInteractions",{}],
    StringStartsQ[#,"T3ExtraM_"] &];
  If[Length[operators] =!= 1, Return[$Failed]];
  result = <|
   "format"->"t3-six-self-native-bilinear-tensor-v1",
   "status"->"NativeMassTensorExtracted",
   "model"->(class<>"_"<>If[alpha<0,"-","+"]<>ToString[Abs[alpha]]),
   "alpha"->alpha,"species"->names,"dimension"->dim,
   "native_metric_inputform"->Map[ToString[#,InputForm]&,identity,{2}],
   "native_metric_factor_inputform"->ToString[factor,InputForm],
   "native_metric_proportional_identity"->True,
   "mass_operator_in_full_lagrangian"->True,
   "operator_names"->operators,
   "rgbeta_matchete_vertex_map_certified"->False,
   "mass_eigenstate_matching_performed"->False
  |>;
  result
]],$Aborted];
If[!AssociationQ[res], Print["Native extraction failed"]; Exit[5]];
If[!DirectoryQ[DirectoryName[out]],CreateDirectory[DirectoryName[out],CreateIntermediateDirectories->True]];
Export[out,res,"RawJSON"];
Print[res["model"],": native metric = ",res["native_metric_factor_inputform"]];

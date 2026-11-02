(* Workstream 4: registered native Matchete Higgs-portal CG tensors, five SU(2)
   representation classes. Completely read-only; no Match[], running, or EFT.

   wolframscript -file scripts/RunLODirectUVNativePortalTensorExport.wl OUTPUT.json

   This child wolframscript enforces the ORIGINAL restricted UV Lagrangian.
   Actual registered tensors come from Matchete`GetCGTensor, not regenerated
   guesses. All components are exported in Matchete's original component axes.
*)
ClearAll["Global`*"];
$WPArgs = Rest[$ScriptCommandLine];
If[Length[$WPArgs] =!= 1, Print["Usage: OUTPUT.json"]; Exit[2]];
$WPOut = ExpandFileName[First[$WPArgs]];
If[FileExistsQ[$WPOut], Print["Refusing to overwrite existing file: ", $WPOut]; Exit[2]];
$WPRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$WPFlags = {"T3_EXPERIMENTAL_GENERAL_FERMIONS",
    "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
    "T3_EXPERIMENTAL_GENERAL_QUARTICS",
    "T3_EXPERIMENTAL_FLAVOR_MASSES"};
$WPInherited = AssociationMap[
  Function[name, With[{v=Environment[name]},If[StringQ[v],v,"UNSET"]]],
  $WPFlags];
$WPOverride = Check[SetEnvironment[Thread[$WPFlags -> ConstantArray["0",Length[$WPFlags]]]]; True, $Failed];
If[$WPOverride =!= True || !AllTrue[$WPFlags,Environment[#] === "0" &],
 Print["BLOCKED: restricted-process environment could not be set"]; Exit[3]];
$WPEffective = AssociationMap[Environment,$WPFlags];
$WPLoad = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"];True],$Failed],$Aborted];
If[$WPLoad =!= True,Print["BLOCKED: Matchete unavailable"];Exit[4]];
Get[FileNameJoin[{$WPRoot,"Lagrangian","model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{$WPRoot,"Lagrangian","model","LagrangianBuilder.wl"}]];

ClearAll[WPError, WPText, WPPortalAssemblyPresentQ, WPRegistered, WPModel];
WPError[message_String] := Throw[<|"status"->"blocked","error"->message|>,"WPAbort"];
WPText[x_] := ToString[FullSimplify[x],InputForm];
(* The previous guard FreeQ[expr,HoldPattern[CG[sym,_List]]] was erroneous:
   HoldPattern holds the local Module variable `sym` rather than evaluating
   it to the actual registered CG symbol. A prior live LInt InputForm export
   confirms the representation CG[T3HiggsPortal...CGInvN, {...}] and
   Coupling[lambdaH...InvN, {}, 0]. Match those exact operator tokens in
   the assembled, checked expression, then query GetCGTensor on its symbol.
   This check does not rely on a generated or guessed tensor component. *)
WPPortalAssemblyPresentQ[expr_,cgName_String,couplingName_String] := Module[
  {form},
  form = ToString[expr,InputForm];
  StringContainsQ[form,"CG["<>cgName<>","] &&
    StringContainsQ[form,"Coupling["<>couplingName<>","]
];
WPRegistered[modelName_String,allowed_List,expr_,which_Integer,d_Integer] := Module[
  {names,results},
  names = Select[allowed, StringMatchQ[#,
    ("HiggsPortal"<>ToString[which]<>"Inv")~~DigitCharacter..] &];
  If[Length[names] < 1 || Length[names] > 2,
    WPError[modelName<>" unexpected portal count for scalar "<>ToString[which]]];
  results = Table[Module[{index,sym,raw,arr,shape,expected,flat},
    index = ToExpression[StringDrop[name,StringLength["HiggsPortal"<>ToString[which]<>"Inv"]]];
    sym = Symbol["T3HiggsPortal"<>ToString[which]<>"CGInv"<>ToString[index]];
    If[!WPPortalAssemblyPresentQ[expr,SymbolName[sym],
        "lambdaH"<>ToString[which]<>"Inv"<>ToString[index]],
      WPError[modelName<>" assembled LInt missing native portal CG and/or coupling "<>
        SymbolName[sym]]];
    raw = CheckAbort[Quiet@Check[GetCGTensor[sym],$Failed],$Aborted];
    If[raw === $Failed || raw === $Aborted,
      WPError[modelName<>" registered GetCGTensor failed: "<>SymbolName[sym]]];
    arr = If[Head[raw] === SparseArray,Normal[raw],raw];
    shape = Dimensions[arr];
    expected = Select[{2,2,d,d},#>1&];
    If[shape =!= expected,
      WPError[modelName<>" native portal tensor shape mismatch for "<>SymbolName[sym]<>
        ": "<>ToString[shape,InputForm]]];
    flat = Flatten[arr];
    If[Length[flat] =!= Times@@expected ||
       !AllTrue[flat,TrueQ[NumericQ[FullSimplify[#]]] &] ||
       !AnyTrue[flat,!TrueQ[PossibleZeroQ[#]] &],
      WPError[modelName<>" incomplete, nonnumeric or zero registered portal "<>SymbolName[sym]]];
    <|"name"->name,"coupling"->("lambdaH"<>ToString[which]<>"Inv"<>ToString[index]),
      "scalar"->which,"invariant"->index,"cg_symbol"->SymbolName[sym],
      "dimensions"->expected,"components_inputform"->(WPText/@flat),
      "tensor_accessor"->"Matchete`GetCGTensor",
      "registered_in_assembled_interaction"->True|>
  ],{name,names}];
  If[Sort[Lookup[results,"invariant"]] =!= Range[Length[results]],
     WPError[modelName<>" portal invariant numbering gap"]];
  results
];
WPModel[class_String,alpha_Integer] := Module[
  {model,built,d1,d2,dF,short,allowed,r1,r2},
  short = class<>"_"<>If[alpha>=0,"+",""]<>ToString[alpha];
  model = T3ModelFromClass[class,alpha];
  If[!AssociationQ[model],WPError[short<>" invalid catalog model"]];
  {d1,d2,dF} = {model["Scalar1","SU2"],model["Scalar2","SU2"],model["Fermion","SU2"]};
  built = CheckAbort[Check[UsingFrontEnd[BuildT3Lagrangian[model]],$Failed],$Aborted];
  If[!AssociationQ[built] || Lookup[built,"Status",""] =!= "Success" ||
     !TrueQ[Lookup[built,"FullValidation",False]],
    WPError[short<>" restricted native Lagrangian did not pass CheckLagrangian"]];
  allowed = Lookup[built,"AllowedInteractions",{}];
  r1 = WPRegistered[short,allowed,built["LInt"],1,d1];
  r2 = WPRegistered[short,allowed,built["LInt"],2,d2];
  <|"model"->short,"class"->class,"alpha"->alpha,
    "dS1"->d1,"dS2"->d2,"dF"->dF,
    "native_build_status"->"Success","native_full_lagrangian_validation"->True,
    "portals"->Join[r1,r2]|>
];
$WPSources = {
  "Lagrangian/model/T3ModelCatalog.wl",
  "Lagrangian/model/LagrangianBuilder.wl",
  "Lagrangian/model/SU2Invariants.wl",
  "Lagrangian/model/T3Fields.wl",
  "Lagrangian/interactions/T3Topology.wl",
  "Lagrangian/interactions/ScalarPotential.wl"};
$WPDigest = Association@Table[With[{file=FileNameJoin[{$WPRoot,path}]},
 If[!FileExistsQ[file],WPError["Missing native source "<>path]];
 path->IntegerString[FileHash[file,"SHA256"],16,64]],{path,$WPSources}];
$WPCases = {{"A",0},{"B",-1},{"C",-1},{"D",0},{"E",0}};
$WPAttempt = CheckAbort[Catch[Table[WPModel@@spec,{spec,$WPCases}],"WPAbort"],$Aborted];
If[!ListQ[$WPAttempt] || Length[$WPAttempt] != 5 ||
   !AllTrue[$WPAttempt,AssociationQ],
  Print["BLOCKED: ",If[AssociationQ[$WPAttempt],Lookup[$WPAttempt,"error","unknown"],
    "one of five native builds failed or was aborted"]]; Exit[5]];
$WPReport = <|
 "format"->"t3-ws4-native-portal-registered-tensors-five-classes-v1",
 "status"->"five_classes_native_registered_portal_tensors_extracted",
 "native_source_sha256"->$WPDigest,
 "restricted_environment"-><|
   "inherited"->$WPInherited,"effective"->$WPEffective,
   "isolated_to_wolframscript_process"->True|>,
 "models"->$WPAttempt,
 "registered_component_axes"->"Matchete native, not T3 weight basis",
 "native_rgbeta_portal_map_certified"->False,
 "spinor_and_T3_mixing_normalisation_certified"->False,
 "matching_or_running_performed"->False|>;
If[!DirectoryQ[DirectoryName[$WPOut]],
  CreateDirectory[DirectoryName[$WPOut],CreateIntermediateDirectories->True]];
If[Check[Export[$WPOut,$WPReport,"RawJSON"],$Failed] === $Failed,
 Print["BLOCKED: report export failed"];Exit[6]];
Print["Native registered portal CG tensors: five classes / ",
  Total[Length[Lookup[#,"portals",{}]]&/@$WPAttempt]," operators"];
Print["Output: ",$WPOut];
Print["Full physical source-to-Weyl normalisation: NOT CERTIFIED"];
Exit[0];

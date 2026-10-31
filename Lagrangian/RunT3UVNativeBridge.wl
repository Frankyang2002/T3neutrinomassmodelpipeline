(* Workstream A: live Matchete baseline topology and Yukawa CG witness.
   Command: wolframscript -file Lagrangian/RunT3UVNativeBridge.wl OUTPUT.json A|B|C|D|E m4|m3|m2|m1|p0|p1|p2
   This extracts registered tensors via GetCGTensor, not regenerated candidates.
   No Match[], RGBeta state mutation, or physical common-basis claim.
*)
ClearAll["Global`*"];
$T3NBArgs=Rest[$ScriptCommandLine];
If[Length[$T3NBArgs]=!=3,Print["Usage: OUTPUT.json FAMILY alphaTag"];Exit[2]];
$T3NBOut=ExpandFileName[$T3NBArgs[[1]]];
$T3NBClass=ToUpperCase[$T3NBArgs[[2]]];
$T3NBTag=$T3NBArgs[[3]];
If[!MemberQ[{"A","B","C","D","E"},$T3NBClass] ||
   !StringMatchQ[$T3NBTag,("m"|"p")~~DigitCharacter..],
  Print["Bad family or alpha tag"];Exit[2]];
$T3NBAlpha=If[StringStartsQ[$T3NBTag,"m"],-1,1]*
  ToExpression[StringDrop[$T3NBTag,1]];
$T3NBModel=$T3NBClass<>"_"<>$T3NBTag;
$T3NBDir=DirectoryName@ExpandFileName[$InputFileName];
$T3NBFlags={"T3_EXPERIMENTAL_GENERAL_FERMIONS",
    "T3_EXPERIMENTAL_FLAVOR_MASSES","T3_EXPERIMENTAL_GENERAL_QUARTICS"};
If[!AllTrue[$T3NBFlags,Environment[#]==="1"&] ||
    Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"]==="1",
  Print["Require general fermion/flavour/quartic flags=1 and quadratics=0"];Exit[2]];
(* The package must be loaded before public Matchete symbols are parsed. *)
$T3NBLoad=CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"];True],$Failed],$Aborted];
If[$T3NBLoad=!=True,Print["Matchete load failed"];Exit[3]];
Get[FileNameJoin[{$T3NBDir,"model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{$T3NBDir,"model","LagrangianBuilder.wl"}]];

ClearAll[T3NBError,T3NBString,T3NBTensor,T3NBExecute];
T3NBError[msg_String]:=Throw[<|"status"->"NativeRegisteredCGExtractionFailed",
   "error"->msg|>,"T3NBError"];
T3NBString[x_]:=ToString[x,InputForm];
T3NBTensor[name_Symbol,expected_List,interaction_] := Module[
  {raw,arr,shape,flat,examples},
  If[FreeQ[interaction,HoldPattern[CG[name,_List]]],
    T3NBError[SymbolName[name]<>" absent from assembled interaction expression"]];
  raw=CheckAbort[Quiet@Check[GetCGTensor[name],$Failed],$Aborted];
  If[raw===$Failed || raw===$Aborted,
    T3NBError["GetCGTensor unavailable for "<>SymbolName[name]]];
  arr=If[Head[raw]===SparseArray,Normal[raw],raw];
  shape=Dimensions[arr];
  If[shape=!=expected,
    T3NBError["Unexpected "<>SymbolName[name]<>" tensor dimensions "<>
      T3NBString[shape]<>" expected "<>T3NBString[expected]]];
  flat=Flatten[arr];
  If[Length[flat]=!=Times@@expected ||
     !AllTrue[flat,TrueQ[NumericQ[FullSimplify[#]]]&],
    T3NBError["Non-numeric or incomplete "<>SymbolName[name]<>" tensor"]];
  If[!AnyTrue[flat,!TrueQ[PossibleZeroQ[#]]&],
    T3NBError["Identically zero "<>SymbolName[name]<>" tensor"]];
  <|"symbol"->SymbolName[name],"dimensions"->shape,
    "tensor_accessor"->"Matchete`GetCGTensor",
    "source"->"live registered Matchete CG",
    "components_inputform"->(T3NBString[FullSimplify[#]]&/@flat),
    "component_count"->Length[flat],
    "registered_in_assembled_interaction"->True|>
];
T3NBExecute[]:=Module[
  {model,d1,d2,dF,built,syms,names,dims,rows,indexed,flavourMass,counts},
  model=T3ModelFromClass[$T3NBClass,$T3NBAlpha];
  If[!AssociationQ[model],T3NBError["Unknown T3 model"]];
  {d1,d2,dF}={model["Scalar1","SU2"],model["Scalar2","SU2"],
     model["Fermion","SU2"]};
  built=BuildT3Lagrangian[model];
  If[!AssociationQ[built] || Lookup[built,"Status",""]=!="Success" ||
     !TrueQ[Lookup[built,"FullValidation",False]],
    T3NBError["Matchete Lagrangian builder did not return Success"]];
  indexed=DeleteDuplicates@Cases[built["LFree"],
    HoldPattern[Coupling[n_,ids_List,___]]/;n===MF:>ids,Infinity];
  flavourMass=indexed=!={} &&
    AllTrue[indexed,MatchQ[#,{Index[_,NFlavor]}]&];
  If[!flavourMass,T3NBError["Flavour indexed MF absent from free sector"]];
  names={"Yukawa1","Yukawa2","TopologyHHSS"};
  syms={T3Y1CG,T3Y2CG,T3MixCG};
  dims={Select[{2,dF,d1},#>1&],Select[{2,dF,d2},#>1&],
    Select[{2,2,d1,d2},#>1&]};
  rows=AssociationThread[names,MapThread[T3NBTensor[#1,#2,built["LInt"]]&, 
      {syms,dims}]];
  If[!AllTrue[Values[rows],AssociationQ],T3NBError["Incomplete CG records"]];
  <|"status"->"LiveRegisteredBaselineCGsCaptured",
    "model"->$T3NBModel,
    "dimensions"->{d1,d2,dF},"alpha"->$T3NBAlpha,
    "native_tensors"->rows,
    "full_matchete_check_passed"->True,
    "indexed_flavour_mf"->True,
    "quartics_included"->Lookup[built,"AdditionalQuarticsIncluded",-1],
    "extra_yukawas_included"->Lookup[built,"AdditionalYukawasIncluded",-1],
    "rgbeta_component_basis_aligned"->False,
    "signed_physical_coupling_map_certified"->False,
    "matching_performed"->False|>
];
$T3NBRun=CheckAbort[
  Check[UsingFrontEnd[Catch[T3NBExecute[],"T3NBError"]],$Failed],
  $Aborted];
If[!AssociationQ[$T3NBRun],$T3NBRun=<|"status"->"WolframKernelFailed"|>];
If[!DirectoryQ[DirectoryName[$T3NBOut]],
 CreateDirectory[DirectoryName[$T3NBOut],CreateIntermediateDirectories->True]];
$T3NBReport=Join[<|"format"->"t3-workstream-a-live-baseline-cg-v1",
  "model"->$T3NBModel|>,$T3NBRun];
$T3NBWrite=Quiet@Check[Export[$T3NBOut,$T3NBReport,"RawJSON"],$Failed];
If[$T3NBWrite===$Failed,Print["Could not write report"];Exit[7]];
Print[$T3NBModel,": ",Lookup[$T3NBReport,"status","Missing"]];
If[Lookup[$T3NBReport,"status",""]=!="LiveRegisteredBaselineCGsCaptured",Exit[6]];
Exit[0];

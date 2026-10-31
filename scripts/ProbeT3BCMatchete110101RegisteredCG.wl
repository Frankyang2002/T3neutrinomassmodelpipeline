(* B/C(-1) 110101: extract TWO actually registered Matchete CG arrays.
   GetCGTensor is Matchete's public tensor accessor; it is NOT the
   T3ExtendedBosonicTensorBasis source-array used by the Stage22 exporter.
   All native calls must remain inside UsingFrontEnd. Fail closed and
   write a diagnostic JSON even on a blocked/missing native CG.

   wolframscript -file scripts/ProbeT3BCMatchete110101RegisteredCG.wl \
      output/audits/t3_bc_matchete_registered_cg/live_registered_cg_v1.json
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: ProbeT3BCMatchete110101RegisteredCG.wl OUTPUT.json"];Exit[64]];
outFile = ExpandFileName[$ScriptCommandLine[[2]]];
If[FileExistsQ[outFile], Print["Refusing to overwrite: ",outFile];Exit[64]];
root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
load = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"];True],$Failed],$Aborted];
If[load =!= True, Print["Matchete/FrontEnd unavailable"];Exit[3]];
Get[FileNameJoin[{root,"Lagrangian","model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{root,"Lagrangian","model","LagrangianBuilder.wl"}]];

ClearAll[T3BCGFail,T3BCGInput,T3BCGOne,T3BCGRun,T3BCGComponent];
T3BCGFail[model_String,msg_String] := Throw[
  <|"model"->model,"status"->"LiveRegisteredMatcheteCGBlocked",
    "error"->msg,"live_registered_cg_components_verified"->False,
    "direct_live_cross_tool_action_identity_verified"->False,
    "physical_matching_certified"->False|>,"T3BCRegisteredCG"];
T3BCGInput[x_] := ToString[x,InputForm];
T3BCGComponent[x_] := T3BCGInput[FullSimplify[x]];

T3BCGOne[label_String,cls_String] := Module[
  {model,sm,baseCG,all,selected,expectedNames,cgNames={},term,
   cg,raw,norm,shape,flat,items={},k,expressionSeen},
  Print["Extracting live Matchete registered 110101 CG for ",label," ..."];
  model=T3ModelFromClass[cls,-1];
  If[!AssociationQ[model] || model["Alpha"]=!=-1 ||
     model["Scalar1","SU2"]=!=2 || model["Scalar2","SU2"]=!=2,
    T3BCGFail[label,"B/C model identity changed"]];
  ResetAll[];
  sm=CheckAbort[Check[LoadModel["SM"],$Failed],$Aborted];
  If[sm===$Failed,T3BCGFail[label,"LoadModel[SM] failed"]];
  If[DefineT3Fields[model]===$Failed || DefineT3Couplings[]=!=True,
    T3BCGFail[label,"field or coupling registration failed"]];
  baseCG=CheckAbort[Check[DefineT3InvariantCGs[model],$Failed],$Aborted];
  If[baseCG===$Failed,T3BCGFail[label,"baseline CG registration failed"]];
  all=CheckAbort[Check[T3ExtendedBuildExtraQuarticCandidates[model],$Failed],$Aborted];
  If[!ListQ[all] || Length[all]=!=7,
    T3BCGFail[label,"incorrect extra quartic inventory"]];
  selected=Select[all,Lookup[#,"Counts",{}]==={1,1,0,1,0,1}&];
  expectedNames={"T3ExtraQ_1_1_0_1_0_1_1","T3ExtraQ_1_1_0_1_0_1_2"};
  If[Length[selected]=!=2 || Lookup[selected,"Name",{}]=!=expectedNames,
    T3BCGFail[label,"incorrect 110101 term selection"]];
  Do[
    cg=Symbol["T3ExtraCGQ110101I"<>ToString[k]];
    term=selected[[k]]["Expression"];
    expressionSeen=!FreeQ[term,CG[cg,_List]];
    If[!TrueQ[expressionSeen],
      T3BCGFail[label,"CG symbol absent from selected live interaction: "<>SymbolName[cg]]];
    (* GetCGTensor reads the CURRENT Matchete registration, not the
       source basis. No fallback to recomputing InvariantTensors is allowed. *)
    raw=CheckAbort[Quiet@Check[GetCGTensor[cg],$Failed],$Aborted];
    If[raw===$Failed || raw===$Aborted,
      T3BCGFail[label,"GetCGTensor failed for "<>SymbolName[cg]]];
    norm=If[Head[raw]===SparseArray,Normal[raw],raw];
    shape=Dimensions[norm];
    If[!ListQ[norm] || shape=!={2,2,2,2},
      T3BCGFail[label,"GetCGTensor returned unexpected head/dimensions for "<>
        SymbolName[cg]<>"; head="<>T3BCGInput[Head[raw]]<>
        "; dimensions="<>T3BCGInput[shape]<>
        "; preview="<>StringTake[T3BCGInput[raw],UpTo[350]]]];
    flat=Flatten[norm];
    If[Length[flat]=!=16 || !AllTrue[flat,NumericQ[FullSimplify[#]]&],
      T3BCGFail[label,"registered CG has non-numeric or incomplete components: "<>
        SymbolName[cg]]];
    AppendTo[items,<|
      "symbol"->SymbolName[cg],
      "matchete_coupling"->("lambdaT3ExtraQ110101I"<>ToString[k]),
      "term_name"->expectedNames[[k]],
      "cg_present_in_native_plus_hc_interaction"->True,
      "tensor_accessor"->"Matchete`GetCGTensor",
      "returned_head"->T3BCGInput[Head[raw]],
      "dimensions"->shape,
      "components"->(T3BCGComponent /@ flat),
      "component_count"->Length[flat]|>],
    {k,1,2}];
  Print[label,": registered CG tensors 2/2, components 32/32 extracted"];
  <|"model"->label,"status"->"LiveRegisteredMatcheteCGComponentsCaptured",
    "extra_quartic_inventory"->Length[all],"registered_tensors"->items,
    "registered_cg_components_extracted"->32,
    "live_registered_cg_components_verified"->False,
    "direct_live_cross_tool_action_identity_verified"->False,
    "physical_matching_certified"->False|>
];

T3BCGRun[] := <|
  "format"->"t3-bc-live-matchete-110101-registered-cg-v1",
  "status"->"LiveRegisteredCGExtractionPendingIndependentComparison",
  "source"->"installed Matchete GetCGTensor called after T3ExtendedBuildExtraQuarticCandidates",
  "models"-><|
    "B_-1"->Catch[T3BCGOne["B_-1","B"],"T3BCRegisteredCG"],
    "C_-1"->Catch[T3BCGOne["C_-1","C"],"T3BCRegisteredCG"]|>,
  "source_basis_tensors_substituted"->False,
  "direct_live_cross_tool_action_identity_verified"->False,
  "physical_coupling_mapping_certified"->False,
  "mass_basis_matching_certified"->False,
  "sequential_eft_matching_certified"->False,
  "production_modified"->False|>;

(* Keep front end active across ResetAll, DefineCG and GetCGTensor. *)
result=CheckAbort[Check[UsingFrontEnd[T3BCGRun[]],$Failed],$Aborted];
If[!AssociationQ[result],
  Print["Live registered CG probe failed before report construction"];Exit[5]];
If[!DirectoryQ[DirectoryName[outFile]],
  CreateDirectory[DirectoryName[outFile],CreateIntermediateDirectories->True]];
write=Quiet@Check[Export[outFile,result,"RawJSON"],$Failed];
If[write===$Failed || !FileExistsQ[outFile],
  Print["Failed to write live CG diagnostic JSON"];Exit[6]];
Print["Registered Matchete CG diagnostic: ",outFile];
If[!AllTrue[Values[result["models"]],
     Lookup[#,"status",""]==="LiveRegisteredMatcheteCGComponentsCaptured"&],
  Print["Live registered CG extraction BLOCKED; diagnostic JSON saved"];Exit[7]];
Print["Live registered Matchete CGs extracted: 64/64; independent comparison pending"];

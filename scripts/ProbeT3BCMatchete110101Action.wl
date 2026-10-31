(* Live Matchete construction/validation of the two B/C(-1) 110101
   extra-quartic LAGRANGIAN terms. Uses the existing T3 builder; no new
   CG, potential, threshold or mass-basis convention is introduced.

   This is deliberately not yet a cross-tool polynomial or matching
   certificate. It tests the Matchete half of that comparison.
   v2: hold the entire model build and CheckLagrangian operations inside
   UsingFrontEnd, not just the initial Matchete package load.

   wolframscript -file scripts/ProbeT3BCMatchete110101Action.wl \
     output/audits/t3_bc_matchete_action/matchete_110101_v1.json
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: ProbeT3BCMatchete110101Action.wl OUTPUT.json"];Exit[64]];
outFile=ExpandFileName[$ScriptCommandLine[[2]]];
If[FileExistsQ[outFile],Print["Refusing to overwrite: ",outFile];Exit[64]];
root=DirectoryName@DirectoryName@ExpandFileName[$InputFileName];
load=CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"];True],$Failed],$Aborted];
If[load=!=True,Print["Matchete/FrontEnd failed to load"];Exit[3]];
Get[FileNameJoin[{root,"Lagrangian","model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{root,"Lagrangian","model","LagrangianBuilder.wl"}]];

ClearAll[T3BCMFail,T3BCMOne,T3BCMText];
T3BCMFail[msg_String] := Throw[<|
  "status"->"Matchete110101ActionBlocked", "error"->msg,
  "matchete_action_validated"->False,
  "cross_tool_action_equality_verified"->False,
  "physical_coupling_mapping_certified"->False,
  "physical_matching_certified"->False|>,"T3BCMatchete"];
T3BCMText[x_,max_Integer:24000] := Module[{str=ToString[x,InputForm]},
  If[StringLength[str]>max,
    T3BCMFail["Matchete expression length exceeds diagnostic limit"]];str];
T3BCMOne[label_String,cls_String] := Module[
  {model,lsm,lfree,cg,all,selected,expected,expressions,individual={},total,
   valid,sourceNames,sourceCouplings,foundSymbols,one,check,status},
  Print["Checking live Matchete 110101 Lagrangian for ",label," ..."];
  model=T3ModelFromClass[cls,-1];
  If[!AssociationQ[model] || model["Scalar1","SU2"]=!=2 ||
     model["Scalar2","SU2"]=!=2 ||
     model["Alpha"]=!=-1,T3BCMFail[label<>" model identity mismatch"]];
  ResetAll[];
  lsm=CheckAbort[Check[LoadModel["SM"],$Failed],$Aborted];
  If[lsm===$Failed,T3BCMFail[label<>" SM load failed"]];
  If[DefineT3Fields[model]===$Failed,T3BCMFail[label<>" field registration failed"]];
  If[CheckAbort[Check[DefineT3Couplings[],$Failed],$Aborted]=!=True,
    T3BCMFail[label<>" baseline coupling registration failed"]];
  cg=CheckAbort[Check[DefineT3InvariantCGs[model],$Failed],$Aborted];
  If[cg===$Failed,T3BCMFail[label<>" baseline CG registration failed"]];
  lfree=CheckAbort[Check[FreeLag[NewFermion,NewScalar1,NewScalar2]//RelabelIndices,$Failed],$Aborted];
  If[lfree===$Failed,T3BCMFail[label<>" free Lagrangian failed"]];
  all=CheckAbort[Check[T3ExtendedBuildExtraQuarticCandidates[model],$Failed],$Aborted];
  If[!ListQ[all] || Length[all]=!=7,T3BCMFail[label<>" extra-quartic inventory is not 7"]];
  selected=Select[all,Lookup[#,"Counts",{}]==={1,1,0,1,0,1}&];
  expected={"T3ExtraQ_1_1_0_1_0_1_1","T3ExtraQ_1_1_0_1_0_1_2"};
  sourceNames=Lookup[selected,"Name",{}];
  If[Length[selected]=!=2 || sourceNames=!=expected ||
     !AllTrue[selected,Lookup[#,"Class",""]==="ExtraQuartic"&],
    T3BCMFail[label<>" incorrect 110101 Matchete interaction selection"]];
  sourceCouplings=Symbol/@{"lambdaT3ExtraQ110101I1","lambdaT3ExtraQ110101I2"};
  Do[
    one=selected[[k]];
    If[FreeQ[one["Expression"],sourceCouplings[[k]]],
      T3BCMFail[label<>" expected coupling absent in native Matchete expression "<>ToString[k]]];
    check=CheckAbort[Check[
      CheckLagrangian[(lsm+lfree+one["Expression"])//RelabelIndices],$Failed],$Aborted];
    AppendTo[individual,<|
      "name"->one["Name"],
      "coupling"->SymbolName[sourceCouplings[[k]]],
      "coupling_present_in_native_expression"->True,
      "check_lagrangian_return_inputform"->T3BCMText[check,500],
      "check_lagrangian_passed"->TrueQ[check],
      "native_expression_inputform"->T3BCMText[one["Expression"]]|>],
    {k,1,2}];
  total=Total[Lookup[selected,"Expression"]]//RelabelIndices;
  valid=CheckAbort[Check[CheckLagrangian[(lsm+lfree+total)//RelabelIndices],$Failed],$Aborted];
  status=If[AllTrue[individual,TrueQ[Lookup[#,"check_lagrangian_passed",False]]&]&&TrueQ[valid],
    "LiveMatchete110101ActionTermsValidated", "Matchete110101ActionValidationFailed"];
  Print[label,": Matchete terms=",Length[selected],
    "; individual CheckLagrangian=",Lookup[individual,"check_lagrangian_passed"],
    "; combined=",TrueQ[valid]];
  <|"status"->status,"model"->label,
    "native_model_built"->True,"scalar_extra_quartics_count"->Length[all],
    "selected_110101_couplings_count"->2,"selected_terms"->individual,
    "combined_check_lagrangian_passed"->TrueQ[valid],
    "combined_check_lagrangian_return_inputform"->T3BCMText[valid,500],
    "source_operator_uses_plus_hc"->"verified separately from GeneralScalarOperators.wl",
    "matchete_action_validated"->(status==="LiveMatchete110101ActionTermsValidated"),
    "cross_tool_action_equality_verified"->False,
    "physical_coupling_mapping_certified"->False,
    "physical_matching_certified"->False|>
];
(* Matchete requires the FrontEnd throughout its model construction and
   CheckLagrangian calls. Wrapping only Needs["Matchete`"] allows the
   FrontEnd session to end before T3BCMOne invokes native operations.
   Follow the existing ExportGeneralT3ScalarTensorMatchete.wl pattern. *)
result=CheckAbort[Check[
  UsingFrontEnd[
  <|"format"->"t3-bc-live-matchete-110101-action-terms-v1",
    "status"->"NativeMatcheteInteractionConstructionDiagnostic",
    "models"-><|
      "B_-1"->Catch[T3BCMOne["B_-1","B"],"T3BCMatchete"],
      "C_-1"->Catch[T3BCMOne["C_-1","C"],"T3BCMatchete"]|>,
    "source"->"Matchete live T3ExtendedBuildExtraQuarticCandidates and CheckLagrangian",
    "matchete_scalar_off_diagonal_mass_included"->False,
    "cross_tool_action_equality_verified"->False,
    "physical_coupling_mapping_certified"->False,
    "physical_matching_certified"->False|>
  ],$Failed],$Aborted];
If[!AssociationQ[result],Print["Could not construct Matchete action diagnostic"];Exit[5]];
If[!DirectoryQ[DirectoryName[outFile]],
  CreateDirectory[DirectoryName[outFile],CreateIntermediateDirectories->True]];
If[Check[Export[outFile,result,"RawJSON"],$Failed]===$Failed ||
   !FileExistsQ[outFile],Print["Cannot write Matchete action diagnostic"];Exit[6]];
Print["Matchete action diagnostic: ",outFile];
If[!AllTrue[Values[result["models"]],
    Lookup[#,"status",""]==="LiveMatchete110101ActionTermsValidated"&],
  Print["Matchete action diagnostic FAILED/BLOCKED; JSON saved. ",outFile];Exit[7]];
Print["Both Matchete 110101 terms validated separately and together."];
Print["Cross-tool potential/action coefficient equality and physical matching: NOT CERTIFIED"];

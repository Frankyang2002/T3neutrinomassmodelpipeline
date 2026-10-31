(* Workstream B: read-only LIVE B/C(-1) mixed quadratic CG extraction.
   Does not run Match[], does not claim a common Matchete/RGBeta field basis.
   wolframscript -file scripts/ProbeT3BCMixedBilinearMatchete.wl OUT.json
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
 Print["Usage: ProbeT3BCMixedBilinearMatchete.wl NEW_OUTPUT.json"];Exit[64]];
out = ExpandFileName[$ScriptCommandLine[[2]]];
If[FileExistsQ[out],Print["Refusing to overwrite ",out];Exit[64]];
root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
loaded = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"];True],$Failed],$Aborted];
If[loaded =!= True,Print["Matchete FrontEnd not available"];Exit[3]];
Get[FileNameJoin[{root,"Lagrangian","model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{root,"Lagrangian","model","LagrangianBuilder.wl"}]];

T3BMFail[model_,message_] := Throw[<|
 "model"->model,"status"->"Blocked","error"->message,
 "physical_mapping_certified"->False|>,"T3BMMixed"];
T3BMOne[label_,cls_] := Module[
 {model,sm,basisCG,structures,candidates,entry,coupling,cg,raw,tensor,
  original,componentTexts,nonZero},
 Print["Live Matchete mixed quadratic for ",label];
 model=T3ModelFromClass[cls,-1];
 If[!AssociationQ[model] || model["Alpha"]=!=-1 ||
   model["Scalar1","SU2"]=!=2 || model["Scalar2","SU2"]=!=2,
   T3BMFail[label,"Unexpected model representation"]];
 ResetAll[];
 sm=CheckAbort[Check[LoadModel["SM"],$Failed],$Aborted];
 If[sm===$Failed,T3BMFail[label,"SM initialization failed"]];
 If[DefineT3Fields[model]===$Failed ||
    CheckAbort[Check[DefineT3Couplings[],$Failed],$Aborted]=!=True,
    T3BMFail[label,"Field/coupling registration failed"]];
 basisCG=CheckAbort[Check[DefineT3InvariantCGs[model],$Failed],$Aborted];
 If[basisCG===$Failed,T3BMFail[label,"Baseline invariant registration failed"]];
 structures=T3ExtendedExtraQuadraticStructures[model];
 If[structures=!={{0,0,0,1,0,1}},
    T3BMFail[label,"Unexpected extra quadratic inventory"]];
 candidates=CheckAbort[Check[
   T3ExtendedBuildCandidates[model,2,structures],$Failed],$Aborted];
 If[!ListQ[candidates] || Length[candidates]=!=1 ||
    (First[candidates])["Counts"]=!={0,0,0,1,0,1} ||
    (First[candidates])["Class"]=!="ExtraQuadratic",
    T3BMFail[label,"Candidate bilinear registration failed"]];
 entry=First[candidates];
 coupling=Symbol["lambdaT3ExtraM000101I1"];
 cg=Symbol["T3ExtraCGM000101I1"];
 If[FreeQ[entry["Expression"],coupling,{0,Infinity}] ||
    FreeQ[entry["Expression"],cg,{0,Infinity}],
    T3BMFail[label,"Registered coupling or CG not in native interaction"]];
 raw=CheckAbort[Quiet@Check[GetCGTensor[cg],$Failed],$Aborted];
 If[raw===$Failed || raw===$Aborted,
    T3BMFail[label,"Live GetCGTensor failed"]];
 tensor=If[Head[raw]===SparseArray,Normal[raw],raw];
 If[Dimensions[tensor]=!={2,2} ||
    !And@@(NumericQ[FullSimplify[#]]& /@ Flatten[tensor]),
    T3BMFail[label,"Live registered CG is not a numerical 2x2 tensor"]];
 nonZero=Count[Flatten[tensor], x_ /; !TrueQ[PossibleZeroQ[x]]];
 If[nonZero=!=2,T3BMFail[label,"Unexpected mixed bilinear tensor rank/structure"]];
 componentTexts=Map[ToString[FullSimplify[#],InputForm]&,tensor,{2}];
 <|"model"->label,"status"->"LiveMatcheteMixedBilinearCGCaptured",
   "source"->"T3ExtendedBuildCandidates degree 2; GetCGTensor on registered CG",
   "counts"->{0,0,0,1,0,1},
   "matchete_coupling"->SymbolName[coupling],
   "matchete_cg"->SymbolName[cg],
   "cg_present_in_actual_plus_hc_expression"->True,
   "tensor_accessor"->"Matchete`GetCGTensor",
   "cg_components_inputform"->componentTexts,
   "matchete_lagrangian_prefactor"->"+1 * muM * CG * Bar[S1] Bar[S2] + h.c.",
   "matchete_source_interaction_inputform"->ToString[entry["Expression"],InputForm],
   "heavy_mass_basis_matching_performed"->False,
   "physical_mapping_certified"->False|>
];
report=CheckAbort[Check[UsingFrontEnd[
 <|"format"->"t3-bc-live-mixed-matchete-bilinear-v1",
  "status"->"LiveNativeCGExtractionDiagnostic",
  "source_sha256"-><|"Lagrangian/interactions/GeneralScalarOperators.wl"->
   IntegerString[FileHash[FileNameJoin[{root,"Lagrangian","interactions","GeneralScalarOperators.wl"}],"SHA256"],16,64]|>,
  "models"-><|
    "B_-1"->Catch[T3BMOne["B_-1","B"],"T3BMMixed"],
    "C_-1"->Catch[T3BMOne["C_-1","C"],"T3BMMixed"]|>,
  "physical_mapping_certified"->False,
  "mass_basis_matching_certified"->False|>],$Failed],$Aborted];
If[!AssociationQ[report],Print["Matchete mixed bilinear probe blocked"];Exit[5]];
If[!DirectoryQ[DirectoryName[out]],
 CreateDirectory[DirectoryName[out],CreateIntermediateDirectories->True]];
If[Check[Export[out,report,"RawJSON"],$Failed]===$Failed || !FileExistsQ[out],
 Print["Cannot write ",out];Exit[6]];
Print["Matchete mixed quadratic diagnostic: ",out];
If[!AllTrue[Values[report["models"]],
  Lookup[#,"status",""]==="LiveMatcheteMixedBilinearCGCaptured"&],Exit[7]];

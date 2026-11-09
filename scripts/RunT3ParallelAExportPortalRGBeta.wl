(* A v5: inspect LIVE RGBeta basic and adjoint Higgs-portal registrations.
   Source is symbolic tGen, not numerical SU2 component generators.
   wolframscript -file scripts/RunT3ParallelAExportPortalRGBeta.wl OUTPUT.json
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: RunT3ParallelAExportPortalRGBeta.wl OUTPUT.json"];Exit[64]];
out=ExpandFileName[$ScriptCommandLine[[2]]];
If[FileExistsQ[out],Print["Refusing to overwrite ",out];Exit[6]];
root=DirectoryName@DirectoryName@ExpandFileName[$InputFileName];
sourceFiles={"RGE/running/rgbeta/T3RGBetaModel.wl",
             "scripts/RGBetaNativeAccess.wl"};
If[!AllTrue[sourceFiles,FileExistsQ[FileNameJoin[{root,#}]]&],
  Print["Required native RGBeta source unavailable"];Exit[2]];
loaded=CheckAbort[Check[Needs["RGBeta`"];True,$Failed],$Aborted];
If[loaded=!=True,Print["RGBeta unavailable"];Exit[3]];
Get[FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
Get[FileNameJoin[{root,"scripts","RGBetaNativeAccess.wl"}]];
T3A5RFail[s_] := (Print["A v5 RGBeta ERROR: ",s];Throw[$Failed,"T3APortalRGB5"]);
T3A5RText[x_] := ToString[x,InputForm];
T3A5RCases=Join[({"A",#}&/@{-4,-2,0}),({"B",#}&/@{-3,-1,1}),
  ({"C",#}&/@{-3,-1,1}),({"D",#}&/@{-2,0,2}),
  ({"E",#}&/@{-4,-2,0,2})];
T3A5RDims[cls_] := Switch[cls,"A",{1,3,2},"B",{2,2,1},
  "C",{2,2,3},"D",{3,1,2},"E",{3,3,2},_,{}];
T3A5RLabel[cl_,a_] := cl<>"_"<>If[a<0,"-","+"]<>ToString[Abs[a]];
T3A5RPortal[which_,dim_,kind_] := Module[
 {symbol,info,f,inv,coup,perms,sc,expected},
 symbol=Symbol["lambdaH"<>ToString[which]<>If[kind==="adjoint","Adj",""]];
 info=T3Stage22NativeInfo[symbol];
 f=T3Stage22NativeGet[info,"Fields"];
 inv=T3Stage22NativeGet[info,"Invariant"];
 coup=T3Stage22NativeGet[info,"Coupling"];
 perms=T3Stage22NativeGet[info,"UniqueArrangements"];
 sc=T3Stage22NativeGet[info,"SelfConjugate"];
 expected={Bar[H],H,Bar[Symbol["S"<>ToString[which]]],Symbol["S"<>ToString[which]]};
 If[!TrueQ[f===expected] || sc=!=True ||
    !ListQ[perms] || Length[perms]<1 || Length[perms]>24 ||
    !AllTrue[perms,Sort[#]==={1,2,3,4}&] ||
    Length[DeleteDuplicates[perms]]=!=Length[perms],
   T3A5RFail["Invalid native RGBeta portal registration "<>SymbolName[symbol]]];
 If[!TrueQ[coup@@{} === 6 symbol],
   T3A5RFail["Nonstandard 6*Matrix coupling convention "<>SymbolName[symbol]]];
 If[kind==="adjoint" && FreeQ[inv,tGen,Infinity],
   T3A5RFail["Adjoint native invariant missing tGen "<>SymbolName[symbol]]];
 <|"symbol"->SymbolName[symbol],"which"->which,"kind"->kind,
   "dimensions"->{2,2,dim,dim},
   "fields_inputform"->T3A5RText[f],
   "native_invariant_inputform"->T3A5RText[inv],
   "native_coupling_inputform"->T3A5RText[coup],
   "native_self_conjugate"->True,
   "native_permutation_count"->Length[perms],
   "native_numerical_generator_components_exported"->False|>
];
T3A5RModel[cl_,a_] := Module[
 {tag,dims,d1,d2,df,build,portals={},which,d,available},
 tag=T3A5RLabel[cl,a];dims=T3A5RDims[cl];
 If[Length[dims]=!=3,T3A5RFail["Unknown model class "<>cl]];
 {d1,d2,df}=dims;
 Print["A v5 RGBeta registered portals: ",tag];
 build=CheckAbort[Check[T3RGBetaBuildEFT1[d1,d2,df,a],$Failed],$Aborted];
 If[!AssociationQ[build],T3A5RFail["Native EFT1 RGBeta build failed "<>tag]];
 Do[
   d=If[which===1,d1,d2];
   AppendTo[portals,T3A5RPortal[which,d,"singlet"]];
   If[d>1,AppendTo[portals,T3A5RPortal[which,d,"adjoint"]]],
   {which,1,2}];
 <|"dimensions"->dims,"alpha"->a,"portals"->portals|>
];
run=CheckAbort[Check[Catch[Catch[Association@Table[
  With[{cl=First[c],a=Last[c]},T3A5RLabel[cl,a]->T3A5RModel[cl,a]],
 {c,T3A5RCases}],"T3Stage22Native"],"T3APortalRGB5"],$Failed],$Aborted];
If[!AssociationQ[run] || Length[Keys[run]]=!=16 ||
   Total[Length[Lookup[#,"portals",{}]]&/@Values[run]]=!=58,
  Print["A v5 RGBeta all-model native portal registrations BLOCKED"];Exit[5]];
sha=Association@Table[
 s->IntegerString[FileHash[FileNameJoin[{root,s}],"SHA256"],16,64],
 {s,sourceFiles}];
report=<|"format"->"t3-parallel-a-native-rgbeta-portal-registrations-v5",
 "status"->"NativeSymbolicPortalRegistrations_ComponentsNotNumerical",
 "model_count"->16,"portal_coupling_slots"->58,
 "source_sha256"->sha,"native_numeric_tGen_components_certified"->False,
 "physical_coupling_map_certified"->False,"models"->run|>;
If[!DirectoryQ[DirectoryName[out]],
 CreateDirectory[DirectoryName[out],CreateIntermediateDirectories->True]];
write=Quiet@Check[Export[out,report,"RawJSON"],$Failed];
If[write===$Failed || !FileExistsQ[out],
 Print["A v5 RGBeta JSON write failed"];Exit[7]];
Print["58/58 symbolic native RGBeta portal registrations (no numeric tGen) -> ",out];

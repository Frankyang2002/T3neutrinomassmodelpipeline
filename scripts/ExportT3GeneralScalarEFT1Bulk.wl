(* Bulk 16-model conditional renormalisable general-scalar EFT1 one-loop
   RGBeta export. Keeps the existing production model unchanged.
   Usage: wolframscript -file scripts/ExportT3GeneralScalarEFT1Bulk.wl MODEL OUTPUT.json
   No Wilson coefficients, scalar quadratic extensions, or certified
   Matchete<->RGBeta coupling normalization are asserted. *)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 3,
  Print["Usage: ... MODEL OUTPUT.json"]; Exit[64]];
model = $ScriptCommandLine[[2]];
out = ExpandFileName[$ScriptCommandLine[[3]]];
root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
If[FileExistsQ[out], DeleteFile[out]];
If[CheckAbort[Quiet @ Check[Needs["RGBeta`"] ; True, $Failed], $Aborted] =!= True,
  Print["Bulk EFT1: Cannot load RGBeta"]; Exit[2]];
T3BulkFail[msg_String] := (Print["Bulk EFT1 ERROR: ",msg];Throw[$Failed,"T3Bulk"]);
T3BulkExport[] := Module[
 {registry,entry,dims,alpha,newRows,reusedRows,names,reused,extraFile,modelFile,
 commonFile,build,baseline,baseNames,qNames,projected,messages={},actual,
 nonQ,baseQ,extraQ,betas,invalid,converted,meta},
 registry=Quiet @ Check[Import[FileNameJoin[{root,"RGE","general","generated",
     "general_rgbeta_quartic_registry.json"}],"RawJSON"],$Failed];
 If[!AssociationQ[registry] || Length[Lookup[registry,"models",<||>]] =!= 16 ||
   !KeyExistsQ[registry["models"],model],T3BulkFail["Invalid 16-model registry/model"]];
 entry=registry["models"][model];
 dims=Lookup[entry,"dimensions",{}];alpha=Lookup[entry,"alpha",Null];
 newRows=Lookup[entry,"new_rgbeta_symbols",{}];
 reusedRows=Lookup[entry,"preexisting_rgbeta_symbols",{}];
 If[!ListQ[dims] || Length[dims]=!=3 || !AllTrue[dims,IntegerQ] ||
    !IntegerQ[alpha] || !ListQ[newRows] || !ListQ[reusedRows],
    T3BulkFail["Invalid model registry representation"]];
 names=Lookup[#,"name",Missing["name"]]& /@ newRows;
 reused=Lookup[#,"existing_name",Missing["existing_name"]]& /@ reusedRows;
 If[!AllTrue[Join[names,reused],StringQ] ||
    Length[Join[names,reused]] =!= Lookup[entry,"total_additional_quartics",-1] ||
    Length[DeleteDuplicates[Join[names,reused]]]=!=Length[Join[names,reused]],
    T3BulkFail["Invalid extra quartic symbol inventory"]];
 extraFile=FileNameJoin[{root,"RGE","general","generated",model<>"_extra_quartics.wl"}];
 modelFile=FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaModel.wl"}];
 commonFile=FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaRunnerCommon.wl"}];
 If[!And@@(FileExistsQ/@{extraFile,modelFile,commonFile}),
    T3BulkFail["Missing existing generated/native source"]];
 Get[commonFile]; Get[modelFile];
 build=Quiet@Check[Apply[T3RGBetaBuildEFT1,Append[dims,alpha]],$Failed];
 If[!AssociationQ[build] ||
    {Lookup[build,"dS1"],Lookup[build,"dS2"],Lookup[build,"dF"],Lookup[build,"alpha"]}
       =!= Append[dims,alpha] || Lookup[build,"ActiveBSMFields",{}]=!={"S1","S2"},
    T3BulkFail["Incorrect native EFT1 representation"]];
 If[CheckAbort[Check[Get[extraFile];True,$Failed],$Aborted]=!=True ||
    !ValueQ[T3GeneralQuarticAddedSymbols] ||
    Sort[SymbolName/@T3GeneralQuarticAddedSymbols]=!=Sort[names] ||
    T3GeneralQuarticExpectedCount=!=Length[names]+Length[reused] ||
    T3GeneralQuarticReusedCount=!=Length[reused],
    T3BulkFail["Generated native quartic vertex registration failed"]];
 baseline=Quiet@Check[T3RGBetaEFT1OneLoopBetas[], $Failed];
 If[!AssociationQ[baseline] || Length[baseline]<15,
    T3BulkFail["Missing scalar EFT1 baseline beta functions"]];
 baseNames=Keys[baseline];
 If[!SubsetQ[baseNames,reused] || Length[Intersection[baseNames,names]]>0,
    T3BulkFail["New/reused RGBeta quartic inventory inconsistent"]];
 Print["Bulk EFT1 ",model,": native full Quartic-class projection, new=",Length[names]];
 projected=CheckAbort[Block[{$MessageList={}},
    Module[{result},result=BetaTerm[Quartic,1];
      messages=DeleteDuplicates[ToString[#,InputForm]&/@$MessageList];result]],$Aborted];
 If[!AssociationQ[projected] || AnyTrue[messages,StringContainsQ[#,"singular"]&],
   Print["Native warnings: ",ToString[Take[messages,UpTo[15]],InputForm]];
   T3BulkFail["Native quartic-class projector failed"]];
 qNames=Select[baseNames,StringStartsQ[#,"lambda"]&];
 actual=Sort[SymbolName/@Keys[projected]];
 If[actual=!=Sort[Join[qNames,names]],
    T3BulkFail["Native projected quartic inventory incomplete: "<>ToString[actual,InputForm]]];
 nonQ=KeySelect[baseline,!StringStartsQ[#,"lambda"]&];
 baseQ=Association@Table[name->Lookup[projected,Symbol[name],Missing["base",name]],{name,qNames}];
 extraQ=Association@Table[name->Lookup[projected,Symbol[name],Missing["new",name]],{name,names}];
 betas=Join[nonQ,baseQ,extraQ];
 invalid[x_]:=x===$Failed||x===$Aborted||x===Indeterminate||MissingQ[x]||
      !FreeQ[HoldComplete[x],_BetaTerm|_QuarticTensors|_CheckProjection];
 If[Sort[Keys[betas]]=!=Sort[Join[baseNames,names]]||
    AnyTrue[Values[betas],invalid],T3BulkFail["Incomplete native coupled one-loop equations"]];
 converted=T3RGBetaConventionalReportBetas[betas];
 meta=<|"SharedScalar"->False,"dS1"->dims[[1]],"dS2"->dims[[2]],
  "dF"->dims[[3]],"alpha"->alpha,"IntegratedField"->"F",
  "ActiveBSMFields"->{"S1","S2"},
  "ReportBetaConvention"->"16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"|>;
 <|"format"->"t3-bulk-general-scalar-eft1-rgbeta-v1","status"->"Success",
   "model"->model,"metadata"->meta,
   "base_beta_symbols"->baseNames,"new_quartic_symbols"->names,
   "reused_quartic_symbols"->reused,
   "physical_coupling_mapping_certified"->False,
   "general_quadratic_terms_evolved"->False,
   "wilson_matching_performed"->False,
   "report_betas"->T3RGBetaSerializeAssociation[converted]|>
];
result=CheckAbort[Catch[T3BulkExport[],"T3Bulk"],$Aborted];
If[!AssociationQ[result] || Lookup[result,"status",""]=!="Success",
 Print["Bulk EFT1 export failed for ",model];Exit[2]];
If[!DirectoryQ[DirectoryName[out]],CreateDirectory[DirectoryName[out],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[out,result,"RawJSON"],$Failed]===$Failed||!FileExistsQ[out],
 Print["Bulk EFT1 JSON export failed"];Exit[3]];
Print["Bulk EFT1 ",model,": ",Length[result["base_beta_symbols"]],
 " baseline + ",Length[result["new_quartic_symbols"]]," new betas EXPORTED"];
Exit[0];

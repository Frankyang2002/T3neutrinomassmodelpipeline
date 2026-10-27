(* Stage 31: ONE ordinary T3 model per Wolfram kernel, 15 models excluding A_+0.
   Full native RGBeta quartic-class projector; also valid when there are zero
   newly generated quartics. Does not touch production model or GitHub.
   wolframscript -file scripts/ExportGeneralT3Stage31OneLoop.wl A_-2 output/.../stage31_A_m2_full_uv_betas.json
*)
ClearAll["Global`*"];
(* Must load RGBeta BEFORE evaluating definitions referring to BetaTerm. *)
If[CheckAbort[Quiet@Check[Needs["RGBeta`"] ; True, $Failed], $Aborted] =!= True,
  Print["Stage31 EXPORT ERROR: Cannot load RGBeta"]; Exit[2]];
If[Length[$ScriptCommandLine] =!= 3,
  Print["Usage: ExportGeneralT3Stage31OneLoop.wl MODEL OUTPUT.json"]; Exit[64]];
model = $ScriptCommandLine[[2]];
out = ExpandFileName[$ScriptCommandLine[[3]]];
root = DirectoryName@DirectoryName@ExpandFileName[$InputFileName];
If[FileExistsQ[out], DeleteFile[out]];
T3S31Fail[msg_String] := (Print["Stage31 EXPORT ERROR: ", msg]; Throw[$Failed,"T3S31"]);
T3S31LoadJSON[path_String] := Module[{obj=Quiet@Check[Import[path,"RawJSON"],$Failed]},
  If[!AssociationQ[obj],T3S31Fail["Missing/invalid input: "<>path]];obj];

T3S31Export[] := Module[
 {reg,all,entry,newRows,reusedRows,names,reused,dims,alpha,ext,modelPath,commonPath,accessPath,
  build,beta,baseKeys,quartic,expectedQ,actualQ,projectedBase,extended,report,
  symbol,info,invalidQ,raw,meta,nativeMessages,nativeResult,expression},
 reg=T3S31LoadJSON[FileNameJoin[{root,"RGE","general","generated",
   "general_rgbeta_quartic_registry.json"}]];
 all=Lookup[reg,"models",<||>];
 If[!AssociationQ[all] || Length[all]=!=16 || model==="A_+0" ||
   !KeyExistsQ[all,model] || !StringMatchQ[model,
      ("A"|"B"|"C"|"D"|"E")~~"_"~~("+"|"-")~~DigitCharacter..],
   T3S31Fail["Expected one of the other 15 registered ordinary T3 models"]];
 entry=all[model];
 dims=Lookup[entry,"dimensions",{}];alpha=Lookup[entry,"alpha",Missing["alpha"]];
 (* Wolfram Lookup[{}, key] returns Missing, not {}. Extract each
    record via Map to preserve empty model-specific symbol lists. *)
 newRows=Lookup[entry,"new_rgbeta_symbols",Missing["new_rgbeta_symbols"]];
 reusedRows=Lookup[entry,"preexisting_rgbeta_symbols",Missing["preexisting_rgbeta_symbols"]];
 If[!ListQ[newRows]||!ListQ[reusedRows]||
    !AllTrue[Join[newRows,reusedRows],AssociationQ],
    T3S31Fail["Invalid registry symbol-record lists for "<>model]];
 names=(Lookup[#,"name",Missing["name"]]& /@ newRows);
 reused=(Lookup[#,"existing_name",Missing["existing_name"]]& /@ reusedRows);
 If[!ListQ[dims]||Length[dims]=!=3||!AllTrue[dims,IntegerQ]||!IntegerQ[alpha]||
   !AllTrue[Join[names,reused],StringQ]||Length[names]+Length[reused]=!=
      Lookup[entry,"total_additional_quartics",-1] ||
   Length[DeleteDuplicates[Join[names,reused]]]=!=Length[names]+Length[reused]||
   !AllTrue[names,StringMatchQ[#,"lambdaT3GeneralQ"~~DigitCharacter..~~"I"~~DigitCharacter..]&],
   T3S31Fail["Invalid model/extension registry inventory"]];
 modelPath=FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaModel.wl"}];
 commonPath=FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaRunnerCommon.wl"}];
 accessPath=FileNameJoin[{root,"scripts","Stage22RGBetaNativeAccess.wl"}];
 ext=FileNameJoin[{root,"RGE","general","generated",model<>"_extra_quartics.wl"}];
 If[!And@@(FileExistsQ/@{modelPath,commonPath,accessPath,ext}),
   T3S31Fail["Missing existing source or generated extension"]];
 Get[commonPath];Get[modelPath];Get[accessPath];
 build=Quiet@Check[Apply[T3RGBetaBuild,Append[dims,alpha]],$Failed];
 If[!AssociationQ[build] || Lookup[build,"dS1",-1]=!=dims[[1]] ||
   Lookup[build,"dS2",-1]=!=dims[[2]] ||
   Lookup[build,"dF",-1]=!=dims[[3]] ||
   Lookup[build,"alpha",Null]=!=alpha,
   T3S31Fail["Model build failed or representation changed"]];
 If[CheckAbort[Check[Get[ext];True,$Failed],$Aborted]=!=True,
   T3S31Fail["Generated quartic extension failed to register"]];
 If[!ValueQ[T3GeneralQuarticAddedSymbols] ||
   Sort[SymbolName/@T3GeneralQuarticAddedSymbols]=!=Sort[names] ||
   T3GeneralQuarticExpectedCount=!=Length[names]+Length[reused] ||
   T3GeneralQuarticReusedCount=!=Length[reused],
   T3S31Fail["Registered extension inventory differs from registry"]];
 Do[symbol=Symbol[name];info=T3Stage22NativeInfo[symbol];
   If[!AssociationQ[info]||T3Stage22NativeGet[info,"SelfConjugate"]=!=False||
      Length[T3Stage22NativeGet[info,"Fields"]]=!=4,
      T3S31Fail["Unexpected native generated quartic: "<>name]],{name,names}];
 Print["Stage31 ",model,": exporting full UV betas; new=",Length[names],
       " reused=",Length[reused]];
 beta=Quiet@Check[T3RGBetaOneLoopBetas[],$Failed];
 If[!AssociationQ[beta]||!And@@(KeyExistsQ[beta,#]&/@
     {"gY","g2","g3","yu","yd","ye","y1","y2","MF","mS1Sq",
      "mS2Sq","lambdaH","lambdaS1","lambdaS2","lambdaH1","lambdaH2",
      "lambda12","lambdaT3"}),
      T3S31Fail["Base model UV beta inventory incomplete"]];
 baseKeys=Keys[beta];
 If[!SubsetQ[baseKeys,reused]||Length[Intersection[names,baseKeys]]>0,
   T3S31Fail["Reused quartics not present in native base beta inventory"]];
 Print["Stage31 ",model,": projecting native full quartic class"];
 nativeMessages={};
 quartic=CheckAbort[Block[{$MessageList={}},
    nativeResult=BetaTerm[Quartic,1];
    nativeMessages=DeleteDuplicates[ToString[#,InputForm]&/@$MessageList];
    nativeResult],$Aborted];
 Print["Stage31 ",model,": quartic head=",Head[quartic],
   "; native messages=",Length[nativeMessages]];
 If[!AssociationQ[quartic],
   If[nativeMessages=!={},Print["Stage31 native message names: ",
      StringRiffle[Take[nativeMessages,UpTo[20]],", "]]];
   T3S31Fail["Class-level quartic projector returned "<>
      ToString[Head[quartic],InputForm]<>"; do not bypass singular projection"]];
 expectedQ=Join[Select[baseKeys,StringStartsQ[#,"lambda"]&],names];
 actualQ=Sort[SymbolName/@Keys[quartic]];
 If[Sort[expectedQ]=!=actualQ,
   T3S31Fail["Projected quartic inventory mismatch expected="<>
      ToString[Sort[expectedQ],InputForm]<>" actual="<>
      ToString[actualQ,InputForm]]];
 projectedBase=Association@Table[name->Lookup[quartic,Symbol[name],
      Missing["MissingBaseQuarticBeta",name]],
      {name,Select[baseKeys,StringStartsQ[#,"lambda"]&]}];
 extended=Association@Table[name->Lookup[quartic,Symbol[name],
      Missing["MissingNewQuarticBeta",name]],{name,names}];
 invalidQ[x_]:=x===$Failed||x===$Aborted||x===Indeterminate||MissingQ[x]||
   !FreeQ[HoldComplete[x],_BetaTerm|_QuarticTensors|_CheckProjection];
 Do[expression=Lookup[Join[projectedBase,extended],name,Missing["beta"]];
    If[invalidQ[expression],T3S31Fail["Missing/unevaluated beta: "<>name]],
    {name,expectedQ}];
 beta=Join[beta,projectedBase];
 report=T3RGBetaConventionalReportBetas[Join[beta,extended]];
 If[Keys[report]=!=Join[baseKeys,names],
   T3S31Fail["Combined beta inventory is inconsistent"]];
 meta=<|"SharedScalar"->False,"dS1"->dims[[1]],"dS2"->dims[[2]],
     "dF"->dims[[3]],"alpha"->alpha,
     "ReportBetaConvention"->"16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"|>;
 raw=<|"format"->"stage31-multi-model-coupled-uv-rgbeta-betas-v1",
    "status"->"Success","model"->model,"metadata"->meta,
    "new_quartic_symbols"->names,"reused_quartic_symbols"->reused,
    "base_beta_symbols"->baseKeys,
    "source"->"Live RGBeta BetaTerm[Quartic,1] projected after model-specific generated quartic registration",
    "physical_coupling_mapping_certified"->False,
    "report_betas"->T3RGBetaSerializeAssociation[report]|>;
 raw
];
(* No surrounding Check[]: RGBeta's Monitor emits harmless FrontEndObject
   warnings in headless wolframscript. Explicit failures remain fail-closed. *)
result=CheckAbort[Catch[T3S31Export[],"T3S31"],$Aborted];
If[!AssociationQ[result]||Lookup[result,"status",""]=!="Success",
  Print["Stage31 full UV beta export FAILED"];Exit[2]];
If[!DirectoryQ[DirectoryName[out]],
  CreateDirectory[DirectoryName[out],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[out,result,"RawJSON"],$Failed]===$Failed||!FileExistsQ[out],
  Print["Stage31 JSON export FAILED"];Exit[3]];
Print["Stage31 ",model,": ",Length[result["base_beta_symbols"]],
  " base + ",Length[result["new_quartic_symbols"]]," new one-loop betas EXPORTED"];
Print["Physical cross-package coupling map: UNCERTIFIED"];
Print["Report: ",out];

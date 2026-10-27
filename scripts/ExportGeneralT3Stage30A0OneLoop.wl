(* Stage 30: full coupled UV one-loop beta export for A_+0 ONLY.
   Does not modify the production RGBeta model. Uses the pre-existing
   generated 5-quartic extension AFTER T3RGBetaBuild, so every beta function
   includes the new interactions. Native class-level quartic projection is
   essential where several invariants have the same field content.
   Numerical basis remains conventional.

   wolframscript -file scripts/ExportGeneralT3Stage30A0OneLoop.wl \
     output/general_t3_stage22/stage30_A_p0_full_uv_betas.json
*)
ClearAll["Global`*"];
(* RGBeta must be loaded before the delayed exporter definition is parsed.
   Otherwise the unqualified BetaTerm in that definition interns a
   Global`BetaTerm symbol, which shadows RGBeta`BetaTerm at Needs[] time.
   Check[] then regards the shadowing warning as a load failure. *)
If[CheckAbort[Quiet@Check[Needs["RGBeta`"] ; True, $Failed], $Aborted] =!= True,
  Print["Stage30 EXPORT ERROR: Cannot load RGBeta package"]; Exit[2]];
If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: ExportGeneralT3Stage30A0OneLoop.wl OUTPUT.json"]; Exit[64]];
out = ExpandFileName[$ScriptCommandLine[[2]]];
root = DirectoryName@DirectoryName@ExpandFileName[$InputFileName];
If[FileExistsQ[out], DeleteFile[out]];
T3S30Fail[msg_String] := (Print["Stage30 EXPORT ERROR: ",msg]; Throw[$Failed,"T3S30"]);
T3S30LoadJSON[path_String] := Module[{obj=Quiet@Check[Import[path,"RawJSON"],$Failed]},
  If[!AssociationQ[obj],T3S30Fail["Missing/invalid input: "<>path]];obj];

T3S30Export[] := Module[
  {reg,entry,names,ext,modelPath,commonPath,accessPath,build,beta,extended,
   baseKeys,extraKeys,report,symbol,info,meta,raw,quartic,expectedQ,actualQ,
   projectedBase,invalidQ,expression,firstInvalid,nativeMessages,nativeResult},
  reg=T3S30LoadJSON[FileNameJoin[{root,"RGE","general","generated",
    "general_rgbeta_quartic_registry.json"}]];
  entry=Lookup[Lookup[reg,"models",<||>],"A_+0",Missing["model"]];
  If[!AssociationQ[entry] || Lookup[entry,"total_additional_quartics",-1]=!=5 ||
      Length[Lookup[entry,"new_rgbeta_symbols",{}]]=!=5 ||
      Length[Lookup[entry,"preexisting_rgbeta_symbols",{}]]=!=0,
    T3S30Fail["A_+0 registry must specify 5 new and 0 reused quartics"]];
  names=Lookup[entry["new_rgbeta_symbols"],"name"];
  If[!ListQ[names]||Length[DeleteDuplicates[names]]=!=5||
     !AllTrue[names,StringMatchQ[#,"lambdaT3GeneralQ"~~DigitCharacter..~~"I"~~DigitCharacter..]&],
     T3S30Fail["Unexpected new-quartic symbol inventory"]];
  modelPath=FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaModel.wl"}];
  commonPath=FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaRunnerCommon.wl"}];
  accessPath=FileNameJoin[{root,"scripts","Stage22RGBetaNativeAccess.wl"}];
  If[!And@@(FileExistsQ/@{modelPath,commonPath,accessPath}),
    T3S30Fail["Missing current RGBeta model/common/native-access source"]];
  Get[commonPath];Get[modelPath];Get[accessPath];
  build=Quiet@Check[T3RGBetaBuild[1,3,2,0],$Failed];
  If[!AssociationQ[build]||Lookup[build,"dS1",0]=!=1||
      Lookup[build,"dS2",0]=!=3||Lookup[build,"dF",0]=!=2||
      Lookup[build,"alpha",1]=!=0,T3S30Fail["A_+0 model build failed"]];
  ext=FileNameJoin[{root,"RGE","general","generated","A_+0_extra_quartics.wl"}];
  If[!FileExistsQ[ext],T3S30Fail["Missing generated A_+0 extension"]];
  If[CheckAbort[Check[Get[ext];True,$Failed],$Aborted]=!=True,
    T3S30Fail["A_+0 generated quartics failed to register"]];
  If[!ValueQ[T3GeneralQuarticAddedSymbols]||
     Sort[SymbolName/@T3GeneralQuarticAddedSymbols]=!=Sort[names]||
     T3GeneralQuarticExpectedCount=!=5||T3GeneralQuarticReusedCount=!=0,
     T3S30Fail["Registered quartic symbols differ from registry"]];
  Do[
    symbol=Symbol[name];info=T3Stage22NativeInfo[symbol];
    If[!AssociationQ[info] || T3Stage22NativeGet[info,"SelfConjugate"]=!=False ||
       Length[T3Stage22NativeGet[info,"Fields"]]=!=4,
       T3S30Fail["Native quartic declaration missing or altered: "<>name]],
    {name,names}];
  Print["Stage30: exporting full base UV beta system with added five quartics"];
  beta=Quiet@Check[T3RGBetaOneLoopBetas[],$Failed];
  If[!AssociationQ[beta]||Length[beta]<15,
    T3S30Fail["Production UV beta export failed"]];
  baseKeys=Keys[beta];
  If[!And@@(KeyExistsQ[beta,#]& /@ {"gY","g2","g3","yu","yd","ye","y1","y2","MF",
      "mS1Sq","mS2Sq","lambdaH","lambdaS1","lambdaS2","lambdaH1",
      "lambdaH2","lambda12","lambdaT3","lambdaH2Adj","lambdaS2Adj"}),
      T3S30Fail["Core A_+0 UV beta inventory incomplete"]];
  (* RGBeta's single-coupling BetaTerm can mix several quartic invariants.
     RGBeta Package/Betas.m explicitly directs the caller to the type-level
     BetaTerm[Quartic,1], which diagonalizes the full projection matrix.
     NEVER use five independent BetaTerm[symbol,1] outputs as physical
     beta functions when mixing is possible. *)
  Print["Stage30: projecting full native RGBeta quartic beta system"];
  (* DO NOT wrap a complete class-level beta calculation in Check[].
     Wolfram Check[] returns $Failed on ANY message, including warnings, even
     when BetaTerm successfully returns a valid, projected Association.
     Leave native messages visible to distinguish harmless warnings from
     RGBeta's actual ProjectionToUnmixedBetas::singular / Abort[]. *)
  nativeMessages={};
  quartic=CheckAbort[
    Block[{$MessageList={}},
      nativeResult=BetaTerm[Quartic,1];
      nativeMessages=DeleteDuplicates[ToString[#,InputForm]& /@ $MessageList];
      nativeResult
    ],
    $Aborted
  ];
  Print["Stage30: class-level RGBeta result head: ",Head[quartic],
    "; native message count: ",Length[nativeMessages]];
  If[!AssociationQ[quartic],
    If[nativeMessages=!={},
      Print["Stage30: native message names: ",
        StringRiffle[Take[nativeMessages,UpTo[20]],", "]]];
    T3S30Fail["RGBeta class-level BetaTerm[Quartic,1] returned "<>
      ToString[Head[quartic],InputForm]<>
      " rather than an Association. Check preceding native messages; "<>
      "a projector singularity cannot be bypassed"]];
  (* The class-level function keys are native coupling Symbols, unlike the
     production base association whose keys are strings. *)
  expectedQ=Join[Select[baseKeys,StringStartsQ[#,"lambda"]&],names];
  actualQ=Sort[SymbolName /@ Keys[quartic]];
  If[Sort[expectedQ]=!=actualQ,
    T3S30Fail["Native class quartic inventory differs from UV model: "<>
      "expected="<>ToString[Sort[expectedQ],InputForm]<>
      " actual="<>ToString[actualQ,InputForm]]];
  (* All base quartics are REPLACED by the same projected class-level betas;
     retaining the unprojected base quartics would break coupled consistency. *)
  projectedBase=Association@Table[
    name->Lookup[quartic,Symbol[name],Missing["AbsentQuarticBeta",name]],
    {name,Select[baseKeys,StringStartsQ[#,"lambda"]&]}];
  extended=Association@Table[
    name->Lookup[quartic,Symbol[name],Missing["AbsentQuarticBeta",name]],
    {name,names}];
  invalidQ[x_]:=x===$Failed || x===$Aborted || x===Indeterminate ||
    MissingQ[x] || !FreeQ[HoldComplete[x],_BetaTerm|_QuarticTensors|_CheckProjection];
  Do[
    expression=extended[name];
    If[invalidQ[expression],
      Print["Stage30: invalid projected new quartic beta: ",name,
        " (",ToString[Short[Unevaluated[expression],3],InputForm],")"];
      T3S30Fail["Projected new quartic beta missing or unevaluated: "<>name]],
    {name,names}];
  Do[
    expression=projectedBase[name];
    If[invalidQ[expression],
      T3S30Fail["Projected base quartic beta missing or unevaluated: "<>name]],
    {name,Keys[projectedBase]}];
  beta=Join[beta,projectedBase];
  extraKeys=Keys[extended];
  If[Length[Intersection[baseKeys,extraKeys]]>0 || Length[extraKeys]=!=5,
    T3S30Fail["New and existing quartic beta symbols overlap"]];
  report=T3RGBetaConventionalReportBetas[Join[beta,extended]];
  If[Keys[report]=!=Join[baseKeys,extraKeys]||
     !AllTrue[Values[report],#=!= $Failed && #=!= $Aborted&],
    T3S30Fail["Conventional beta conversion failed"]];
  meta=<|"SharedScalar"->False,"dS1"->1,"dS2"->3,"dF"->2,"alpha"->0,
    "ReportBetaConvention"->"16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"|>;
  raw=<|"format"->"stage30-A-p0-coupled-uv-rgbeta-betas-v1",
    "status"->"Success","model"->"A_+0","metadata"->meta,
    "new_quartic_symbols"->names,"base_beta_symbols"->baseKeys,
    "source"->"Live RGBeta BetaTerm[Quartic,1] projected after A_+0 quartic registration; ordinary nonquartic BetaTerm",
    "triplet_component_basis"->"Stage27 provisional reference",
    "physical_coupling_mapping_certified"->False,
    "report_betas"->T3RGBetaSerializeAssociation[report]|>;
  raw
];
(* Do NOT wrap the exporter in Check[]. RGBeta's native
   ProjectionToUnmixedBetas uses Monitor[], which emits
   FrontEndObject::notavail under headless wolframscript even when the
   projected beta system is a valid Association. A surrounding Check[]
   would discard that valid result. Real failures still propagate via
   T3S30Fail/Catch, CheckAbort, and the strict result validation below. *)
result=CheckAbort[Catch[T3S30Export[],"T3S30"],$Aborted];
If[!AssociationQ[result]||Lookup[result,"status",""]=!="Success",
  Print["Stage30 full UV beta export FAILED"];Exit[2]];
If[!DirectoryQ[DirectoryName[out]],
  CreateDirectory[DirectoryName[out],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[out,result,"RawJSON"],$Failed]===$Failed || !FileExistsQ[out],
  Print["Stage30 JSON export FAILED"];Exit[3]];
Print["Stage30 full coupled RGBeta beta export: ",Length[result["base_beta_symbols"]],
  " base + 5 new beta functions (native quartic-class projected)"];
Print["Physical cross-package coupling map: UNCERTIFIED"];
Print["Report: ",out];

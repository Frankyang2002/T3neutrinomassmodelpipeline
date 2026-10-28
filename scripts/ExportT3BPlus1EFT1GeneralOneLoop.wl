(* Batch 20: isolated B(+1) post-fermion renormalisable EFT1 native RGBeta
   beta export, including lambdaT3GeneralQ022000I1.

   wolframscript -file scripts/ExportT3BPlus1EFT1GeneralOneLoop.wl OUTPUT.json

   This is NOT Wilson matching or running, nor an independent physical
   RGBeta/Matchete coupling normalisation certificate. *)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: ExportT3BPlus1EFT1GeneralOneLoop.wl OUTPUT.json"]; Exit[64]];
out = ExpandFileName[$ScriptCommandLine[[2]]];
root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
If[FileExistsQ[out], DeleteFile[out]];
If[CheckAbort[Quiet @ Check[Needs["RGBeta`"]; True, $Failed], $Aborted] =!= True,
  Print["Batch20: Cannot load RGBeta"]; Exit[2]];
T3B20Fail[msg_String] := (Print["Batch20 native EFT1 EXPORT ERROR: ", msg];
  Throw[$Failed, "T3B20"]);

T3B20Export[] := Module[
  {reg, entry, names, modelFile, commonFile, extraFile, build, baseline,
   baseNames, qNames, projected, msgs={}, invalid, qActual, expected,
   nonQ, baseQ, extraQ, betas, converted, meta},
  reg = Quiet @ Check[Import[FileNameJoin[{root,"RGE","general","generated",
      "general_rgbeta_quartic_registry.json"}], "RawJSON"], $Failed];
  If[!AssociationQ[reg] || Length[Keys[Lookup[reg,"models",<||>]]] =!= 16,
    T3B20Fail["Missing 16-model quartic registry"]];
  entry = Lookup[reg["models"], "B_+1", Missing["B_+1"]];
  If[!AssociationQ[entry] || Lookup[entry,"dimensions",{}] =!= {2,2,1} ||
      Lookup[entry,"alpha",Null] =!= 1 ||
      Lookup[entry,"total_additional_quartics",-1] =!= 1 ||
      Lookup[entry,"preexisting_rgbeta_symbols",Missing["reused"]] =!= {},
    T3B20Fail["B(+1) registry representation/extra-quartic count changed"]];
  names = Lookup[Lookup[entry,"new_rgbeta_symbols",{}],"name",{}];
  If[names =!= {"lambdaT3GeneralQ022000I1"},
    T3B20Fail["B(+1) expected extra quartic changed"]];
  modelFile = FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaModel.wl"}];
  commonFile = FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaRunnerCommon.wl"}];
  extraFile = FileNameJoin[{root,"RGE","general","generated","B_+1_extra_quartics.wl"}];
  If[!And @@ (FileExistsQ /@ {modelFile,commonFile,extraFile}),
    T3B20Fail["Missing RGBeta model/common/generated vertex source"]];
  Get[commonFile]; Get[modelFile];
  build = Quiet @ Check[T3RGBetaBuildEFT1[2,2,1,1], $Failed];
  If[!AssociationQ[build] || Lookup[build,"dS1",0] =!= 2 ||
      Lookup[build,"dS2",0] =!= 2 || Lookup[build,"dF",0] =!= 1 ||
      Lookup[build,"alpha",0] =!= 1 ||
      Lookup[build,"ActiveBSMFields",{}] =!= {"S1","S2"},
    T3B20Fail["B(+1) scalar-only EFT1 build failed"]];
  (* Register precisely the same native quartic vertex used in Stage31 UV,
     this time AFTER building the F-absent EFT1 model. *)
  If[CheckAbort[Check[Get[extraFile]; True, $Failed],$Aborted] =!= True ||
      !ValueQ[T3GeneralQuarticAddedSymbols] ||
      T3GeneralQuarticAddedSymbols =!= {lambdaT3GeneralQ022000I1} ||
      T3GeneralQuarticExpectedCount =!= 1 || T3GeneralQuarticReusedCount =!= 0,
    T3B20Fail["Generated B(+1) vertex not registered exactly once"]];
  baseline = Quiet @ Check[T3RGBetaEFT1OneLoopBetas[], $Failed];
  If[!AssociationQ[baseline] || Length[baseline] =!= 18,
    T3B20Fail["Native EFT1 baseline 18 beta slots absent"]];
  baseNames = Keys[baseline];
  expected = {"gY","g2","g3","yu","yd","ye","mS1Sq","mS2Sq",
    "lambdaH","lambdaS1","lambdaS2","lambdaH1","lambdaH2",
    "lambda12","lambdaT3","lambdaH1Adj","lambdaH2Adj","lambda12Adj"};
  If[Sort[baseNames] =!= Sort[expected],
    T3B20Fail["Unexpected baseline EFT1 beta symbols"]];
  Print["Batch20: native B(+1) EFT1 Quartic class projection started"];
  projected = CheckAbort[Block[{$MessageList={}},
    Module[{raw},raw=BetaTerm[Quartic,1];
      msgs=DeleteDuplicates[ToString[#,InputForm]& /@ $MessageList]; raw]], $Aborted];
  If[!AssociationQ[projected] ||
      AnyTrue[msgs,StringContainsQ[#,"singular"]&],
    Print["Native messages: ",ToString[Take[msgs,UpTo[12]],InputForm]];
    T3B20Fail["Native full quartic-class projection failed"]];
  qNames = Select[baseNames,StringStartsQ[#,"lambda"]&];
  qActual = Sort[SymbolName /@ Keys[projected]];
  If[qActual =!= Sort[Join[qNames,names]],
    T3B20Fail["Projected quartics do not cover exactly 10 base + 1 extra symbols: " <>
      ToString[qActual,InputForm]]];
  nonQ = KeySelect[baseline,!StringStartsQ[#,"lambda"]&];
  baseQ = Association @ Table[name -> Lookup[projected,Symbol[name],
      Missing["MissingQuartic",name]],{name,qNames}];
  extraQ = <|First[names] -> Lookup[projected,lambdaT3GeneralQ022000I1,
      Missing["MissingExtra"]]|>;
  betas = Join[nonQ,baseQ,extraQ];
  invalid[expr_] := expr === $Failed || expr === $Aborted ||
      expr === Indeterminate || MissingQ[expr] ||
      !FreeQ[HoldComplete[expr],_BetaTerm|_QuarticTensors|_CheckProjection];
  If[Sort[Keys[betas]] =!= Sort[Join[baseNames,names]] ||
      AnyTrue[Values[betas],invalid],
    T3B20Fail["One or more projected coupled EFT1 beta equations missing"]];
  converted = T3RGBetaConventionalReportBetas[betas];
  meta = <|"SharedScalar"->False,"dS1"->2,"dS2"->2,"dF"->1,
    "alpha"->1,"IntegratedField"->"F","ActiveBSMFields"->{"S1","S2"},
    "ReportBetaConvention"->"16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"|>;
  <|"format"->"batch20-b-p1-coupled-eft1-rgbeta-betas-v1",
    "status"->"Success","model"->"B_+1","metadata"->meta,
    "base_beta_symbols"->baseNames,"new_quartic_symbols"->names,
    "physical_coupling_mapping_certified"->False,
    "wilson_matching_performed"->False,
    "stage29_quartic_evolved_in_this_export"->False,
    "source"->"Native EFT1 RGBeta BetaTerm[Quartic,1] with Stage31 B(+1) generated vertex",
    "report_betas"->T3RGBetaSerializeAssociation[converted]|>
];
result = CheckAbort[Catch[T3B20Export[],"T3B20"],$Aborted];
If[!AssociationQ[result] || Lookup[result,"status",""] =!= "Success",
  Print["Batch20 EFT1 native export FAILED"];Exit[2]];
If[!DirectoryQ[DirectoryName[out]],
  CreateDirectory[DirectoryName[out],CreateIntermediateDirectories->True]];
If[Quiet @ Check[Export[out,result,"RawJSON"],$Failed] === $Failed || !FileExistsQ[out],
  Print["Batch20 JSON export failed"];Exit[3]];
Print["Batch20 B(+1) native EFT1 export PASSED: 18 baseline + 1 extra betas -> ",out];
Print["Matchete/RGBeta physical mapping UNCERTIFIED; Wilson matching NOT PERFORMED"];
Exit[0];

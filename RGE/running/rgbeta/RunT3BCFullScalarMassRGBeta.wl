(* Workstream B: LIVE one-loop UV/EFT1 Higgs + B/C(-1) mixed scalar masses.
   This is an independent, opt-in experiment. It does not change
   T3RGBetaModel.wl or the production RG state.

   wolframscript -file RGE/running/rgbeta/RunT3BCFullScalarMassRGBeta.wl \
       EFT1 B output/audits/workstream_b_bc_full_mass/B_m1_eft1_full_mass_v1.json
   The four complex parameter betas represent five real mass equations:
   mHSq, mS1Sq, mS2Sq, Re[muT3GeneralS1S2], Im[muT3GeneralS1S2].
   All scalar quartics are actually registered before evaluating the betas.
   No mass-basis threshold matching or full Matchete/RGBeta phase certificate.

   v2: RGBeta's beta evaluation may emit FrontEndObject::notavail when
   WolframScript is called without a front-end context. Running with
   UsingFrontEnd follows the locally successful Matchete invocation pattern.
   DO NOT wrap BetaTerm in Check: unrelated FrontEnd messages can cause Check
   to discard an otherwise valid one-loop expression. Explicitly validate
   the actual symbolic result instead, and write a blocked JSON diagnostic
   if a beta remains unresolved. All previous success metadata stays intact.
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 4,
 Print["Usage: ... UV|EFT1 B|C NEW_OUTPUT.json"];Exit[64]];
stage = ToUpperCase[$ScriptCommandLine[[2]]];
cls = ToUpperCase[$ScriptCommandLine[[3]]];
out = ExpandFileName[$ScriptCommandLine[[4]]];
(* WolframScript -file counts the script itself as argument 1, so four total
   entries means script, STAGE, CLASS, OUT. *)
If[!MemberQ[{"UV", "EFT1"},stage] || !MemberQ[{"B","C"},cls],
 Print["Unsupported stage/model: ",stage,"/",cls];Exit[64]];
If[FileExistsQ[out],Print["Refusing to overwrite ",out];Exit[64]];
root = DirectoryName @ DirectoryName @ DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];
load = CheckAbort[Check[Needs["RGBeta`"];True,$Failed],$Aborted];
If[load =!= True,Print["RGBeta unavailable"];Exit[3]];
Get[FileNameJoin[{scriptDir,"T3RGBetaRunnerCommon.wl"}]];
Get[FileNameJoin[{scriptDir,"T3RGBetaModel.wl"}]];
T3BHMFail[msg_, details_:<||>] := (
 Print["B/C Higgs mass RGE: ",msg];
 Throw[Join[<|
  "format"->"workstream-b-bc-full-scalar-mass-live-rgbeta-v1",
  "status"->"Blocked",
  "model"->cls<>"_-1",
  "stage"->stage,
  "error"->msg,
  "physical_mapping_certified"->False,
  "mass_eigenstate_threshold_matching_performed"->False
 |>, details],"T3BHM"]
);

T3BHMRun[] := Module[
 {dF,label,model,extension,reg,entry,registered,massRecord,
  massSymbols,raw,reported,metadata,invalid,expectedNew,needed,
  fileSource,baseSource,sourceHashes,failedBetas},
 dF = If[cls === "B",1,3];
 label = cls<>"_-1";
 Print["Live B/C full scalar mass beta: ",label," ",stage];
 model = CheckAbort[Check[
   If[stage === "UV",T3RGBetaBuild[2,2,dF,-1],
                      T3RGBetaBuildEFT1[2,2,dF,-1]],$Failed],$Aborted];
 If[!AssociationQ[model],T3BHMFail["RGBeta baseline build failed"]];
 If[KeyExistsQ[$scalarMasses,mHSq],
   T3BHMFail["mHSq was already registered in the model; refuse duplicate"]];

 (* The SM Higgs is one complex SU(2) doublet; its real mass parameter
    multiplies H^dag H in the RGBeta scalar-potential convention.
    Registration must precede BetaTerm, and occurs after ResetModel[]. *)
 registered = CheckAbort[Check[
   AddScalarMass[mHSq,{Bar@H,H},
     GroupInvariant -> (del[SU2L@fund,#1,#2]&),
     SelfConjugate -> True];True,$Failed],$Aborted];
 If[registered =!= True || !KeyExistsQ[$scalarMasses,mHSq],
   T3BHMFail["Higgs quadratic mass registration failed"]];
 massRecord = $scalarMasses[mHSq];
 If[!AssociationQ[massRecord],T3BHMFail["Higgs native mass record missing"]];

 reg = Quiet@Check[Import[FileNameJoin[{root,"RGE","general","generated",
                 "general_rgbeta_quartic_registry.json"}],"RawJSON"],$Failed];
 If[!AssociationQ[reg] || !KeyExistsQ[Lookup[reg,"models",<||>],label],
   T3BHMFail["General scalar registry missing B/C model"]];
 entry = reg["models",label];
 If[Lookup[entry,"dimensions",{}] =!= {2,2,dF} ||
    Lookup[entry,"alpha",Null] =!= -1 ||
    Lookup[entry,"total_additional_quartics",0] =!= 7,
   T3BHMFail["B/C model registry changed"]];
 extension = FileNameJoin[{root,"RGE","general","generated",
               label<>"_extra_quartics.wl"}];
 T3GeneralQuarticAddedSymbols = {};
 registered = CheckAbort[Check[Get[extension];True,$Failed],$Aborted];
 expectedNew = {lambdaT3GeneralQ110101I1};
 If[registered =!= True || T3GeneralQuarticAddedSymbols =!= expectedNew ||
    T3GeneralQuarticReusedCount =!= 6 ||
    !KeyExistsQ[$quartics,First[expectedNew]],
   T3BHMFail["Seven-quartic RGBeta extension incomplete"]];

 registered = CheckAbort[Check[
   AddScalarMass[muT3GeneralS1S2,{Bar@S1,Bar@S2},
     GroupInvariant -> (eps[SU2L@fund,#1,#2]&),
     SelfConjugate -> False];True,$Failed],$Aborted];
 If[registered =!= True || !KeyExistsQ[$scalarMasses,muT3GeneralS1S2],
   T3BHMFail["B/C mixed mass registration failed"]];

 massSymbols = {mHSq,mS1Sq,mS2Sq,muT3GeneralS1S2};
 (* RGBeta produces FrontEnd messages during beta-function evaluation in
    some WolframScript sessions. Check[BetaTerm[...],$Failed] is too broad:
    it returns $Failed if ANY message is emitted, even when the expression
    itself is valid. Quiet suppresses diagnostic display only; any actual
    failure value, unresolved BetaTerm or nonfinite expression is rejected. *)
 raw = Association @ Table[
   With[{sym = massSymbols[[i]]},
     Print["Computing one-loop mass beta: ",SymbolName[sym]];
     SymbolName[sym] -> CheckAbort[Quiet[BetaTerm[sym,1]],$Aborted]],
   {i,Length[massSymbols]}];
 failedBetas = Select[Keys[raw], Function[name,
   With[{value=raw[name]},
     value === $Failed || value === $Aborted ||
     !FreeQ[value,BetaTerm|_Missing|Indeterminate|ComplexInfinity]
   ]]];
 If[Sort[Keys[raw]] =!= Sort[SymbolName /@ massSymbols] ||
    failedBetas =!= {},
   T3BHMFail["A required one-loop mass beta is unresolved",
     <|"beta_diagnostics" -> Association @ Table[
       name -> <|
          "head" -> ToString[Head[raw[name]],InputForm],
          "inputform_excerpt" -> StringTake[
             ToString[raw[name],InputForm],UpTo[1200]],
          "unresolved" -> MemberQ[failedBetas,name]
       |>,{name,Keys[raw]}]|>]];
 Print["All four mass beta expressions returned concrete results."];
 (* This convention is identical to the project's full-quartic RGE runner.
    No independent individual-quartic beta projections are used here. *)
 reported = T3RGBetaConventionalReportBetas[raw];
 metadata = Join[
   T3RGBetaCommonMetadata[<|
      "SharedMode"->False,"dS1"->2,"dS2"->2,"dF"->dF,
      "alpha"->-1|>,model],
   <|"GeneralScalarModel"->label,"GeneralScalarStage"->stage,
     "GeneralScalarQuadraticMassMixing"->"IncludedCandidate",
     "RegisteredQuarticCount"->7,
     "HiggsScalarMassSymbol"->"mHSq",
     "MassBetaParameters"->(SymbolName/@massSymbols),
     "QuadraticPotentialConvention"->
       "V2=mHSq Bar[H]H+mS1Sq Bar[S1]S1+mS2Sq Bar[S2]S2+mu Bar[S1]Bar[S2]eps+h.c.; no extra half",
     "ReportBetaConvention"->
       "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g",
     "PhysicalMatcheteRGBetaMassBasisCertified"->False,
     "FullRenormalizableMassSectorClosed"->"PendingIndependentVerification",
     "DimensionlessCouplingsRunInThisReport"->False|>];
 <|"format"->"workstream-b-bc-full-scalar-mass-live-rgbeta-v1",
   "status"->"LiveHiggsAndMixedScalarMassBetasExtractedNotPhysicalCertificate",
   "model"->label,"stage"->stage,
   "metadata"->metadata,
   "registered_mass_symbols"->SymbolName/@massSymbols,
   "report_betas"->T3RGBetaSerializeAssociation[reported],
   "mass_eigenstate_threshold_matching_performed"->False,
   "physical_mapping_certified"->False|>
];
(* A working Mathematica FrontEnd was already used by the project's live
   Matchete probes on Windows. RGBeta may invoke front-end progress monitors.
   We avoid the outer Check[...] here for the same message-versus-result reason.
   Failures are checked by inspecting the returned Association. *)
result = CheckAbort[Catch[UsingFrontEnd[T3BHMRun[]],"T3BHM"],$Aborted];
If[!AssociationQ[result],
 Print["Live B/C mass RGE blocked before producing a report"];Exit[5]];
If[!DirectoryQ[DirectoryName[out]],
 CreateDirectory[DirectoryName[out],CreateIntermediateDirectories->True]];
write=Quiet@Check[Export[out,result,"RawJSON"],$Failed];
If[write === $Failed || !FileExistsQ[out],
 Print["Failed writing B/C scalar mass JSON"];Exit[6]];
If[Lookup[result,"status",""] === "Blocked",
 Print["Live B/C mass RGE remains BLOCKED; see beta_diagnostics: ",out];
 Exit[5]];
Print["Five real mass equations captured: ",result["model"]," ",stage];
Print["Report: ",out];
Print["Physical mass-basis thresholds: NOT CERTIFIED"];
Exit[0];

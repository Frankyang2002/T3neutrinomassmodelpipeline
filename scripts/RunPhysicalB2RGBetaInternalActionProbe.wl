(* B2: READ-ONLY fresh-kernel RGBeta internal coupling normalization witness.
   It probes AddYukawa/AddFermionMass/AddQuartic registry functions, NOT the
   physical Lagrangian density, Feynman rules, Matchete gamma phase, or numeric CG.
   Usage: wolframscript -file scripts/RunPhysicalB2RGBetaInternalActionProbe.wl OUTPUT.json
   OUTPUT must not exist. All T3_EXPERIMENTAL_GENERAL_* flags must be off.
   Invoke through scripts/RunPhysicalB2RGBetaInternalActionNative.py so
   Python supplies source-verified Git metadata via environment. Wolfram
   must not execute git or infer Git provenance itself.
*)
ClearAll["Global`*"];
If[Length[Rest[$ScriptCommandLine]] =!= 1,
 Print["Usage: OUTPUT.json"];Exit[2]];
$B4Out=ExpandFileName[Last[$ScriptCommandLine]];
If[FileExistsQ[$B4Out],Print["Refusing existing output: ",$B4Out];Exit[2]];
$B4Env={"T3_EXPERIMENTAL_GENERAL_FERMIONS",
 "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
 "T3_EXPERIMENTAL_GENERAL_QUARTICS",
 "T3_EXPERIMENTAL_FLAVOR_MASSES"};
If[AnyTrue[$B4Env,Environment[#]==="1"&],
 Print["BLOCKED: nonrestricted environment"];Exit[3]];
$B4Root=DirectoryName@DirectoryName@ExpandFileName[$InputFileName];
If[CheckAbort[Check[Needs["RGBeta`"];True,$Failed],$Aborted]=!=True,
 Print["BLOCKED: RGBeta package unavailable"];Exit[4]];
Get[FileNameJoin[{$B4Root,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
$B4Models=<|"A"->{1,3,2,0},"B"->{2,2,1,-1},
 "C"->{2,2,3,-1},"D"->{3,1,2,0},"E"->{3,3,2,0}|>;
$B4Sources={"RGE/running/rgbeta/T3RGBetaModel.wl",
 "Lagrangian/interactions/T3Topology.wl",
 "Lagrangian/model/SU2Invariants.wl"};
(* Python verifies tracked source blobs, then passes HEAD and clean-filter
   hashes to Wolfram. Missing provenance blocks before native evaluation. *)
$B4HeadFromPython=Environment["T3_B2_PROVENANCE_HEAD"];
$B4BlobsFromPython=<|
  "RGE/running/rgbeta/T3RGBetaModel.wl"->Environment["T3_B2_PROVENANCE_RGBETA"],
  "Lagrangian/interactions/T3Topology.wl"->Environment["T3_B2_PROVENANCE_TOPOLOGY"],
  "Lagrangian/model/SU2Invariants.wl"->Environment["T3_B2_PROVENANCE_SU2"]|>;
B4SHA1Q[x_] := StringQ[x] && StringMatchQ[x,RegularExpression["[0-9a-f]{40}"]];
If[!B4SHA1Q[$B4HeadFromPython] ||
   !AllTrue[Values[$B4BlobsFromPython],B4SHA1Q],
  Print["BLOCKED: missing Python-verified Git provenance. Use the Python native launcher."];
  Exit[3]];
ClearAll[B4String,B4Key,B4Read,B4CouplingRow,B4Family,
 B4GroupInvDownValuesString,B4JSONLeafQ,B4JSONUnexpected,
 B4JSONSafe];
B4String[x_]:=ToString[x,InputForm];
B4Key[x_]:=Last@StringSplit[B4String[x],"`"];
B4Read[rec_,name_String]:=Module[{k},
 If[!AssociationQ[rec],Return[Missing["NoAssociation",name]]];
 k=Select[Keys[rec],B4Key[#]===name&];
 If[Length[k]=!=1,Return[Missing["MissingOrAmbiguous",name]]];
 Lookup[rec,First[k]]];

(* DownValues has HoldAll: DownValues[ToExpression[name]] tries to inspect
   the unevaluated ToExpression[...] expression, raising DownValues::sym.
   With inserts the resolved symbol before DownValues evaluates. *)
B4GroupInvDownValuesString[names_List] := If[Length[names]===1,
 With[{symbol=ToExpression[First[names]]}, B4String[DownValues[symbol]]],
 "unresolved"];
B4CouplingRow[registry_,sym_Symbol,kind_String,expected_,expectedFields_List]:=Module[
 {rec,fields,inv,chir,sc,c,cb,expr,exprbar,base,basebar,ratio,ratiobar,
   ok,barSymbol,keys,unique,params,arity},
 rec=If[AssociationQ[registry],Lookup[registry,sym,Missing["Absent"]],Missing["NoRegistry"]];
 If[!AssociationQ[rec],Return[<|"coupling"->SymbolName[sym],"status"->"MissingRegistry"|>]];
 fields=B4Read[rec,"Fields"];inv=B4Read[rec,"Invariant"];
 chir=B4Read[rec,"Chirality"];sc=B4Read[rec,"SelfConjugate"];
 c=B4Read[rec,"Coupling"];cb=B4Read[rec,"CouplingBar"];
 unique=B4Read[rec,"UniqueArrangements"];
 params=If[kind==="quartic",{}, {b4p,b4r}];
 barSymbol=If[kind==="quartic",sym,Bar[sym]];
 base=If[params==={},Matrix[barSymbol][],Matrix[barSymbol][b4p,b4r]];
 basebar=If[params==={},Matrix[Bar[barSymbol]][],Matrix[Bar[barSymbol]][b4p,b4r]];
 (* Quartics use unbarred lambda for Coupling, conjugate for CouplingBar;
    y1,y2,MF have Chirality Right and reverse this orientation. *)
 expr=If[Head[c]===Missing,$Failed,Quiet@Check[Apply[c,params],$Failed]];
 exprbar=If[Head[cb]===Missing,$Failed,Quiet@Check[Apply[cb,params],$Failed]];
 ratio=If[expr===$Failed,Missing["NoCoupling"],Quiet@Check[FullSimplify[expr/base],$Failed]];
 ratiobar=If[exprbar===$Failed,Missing["NoConjugate"],Quiet@Check[FullSimplify[exprbar/basebar],$Failed]];
 ok=SameQ[fields,expectedFields] && Head[inv]=!=Missing &&
   ratio=!=$Failed && ratiobar=!=$Failed &&
   TrueQ[FullSimplify[ratio==expected]] &&
   TrueQ[FullSimplify[ratiobar==expected]] &&
   If[kind==="quartic",TrueQ[sc=== (sym=!=lambdaT3)],chir===Right];
 <|"coupling"->SymbolName[sym],"kind"->kind,
   "status"->If[ok,"NativeRegistryRatioMeasured","BlockedUnexpectedRegistry"],
   "registered_fields"->B4String[fields],
   "actual_expected_fields_match"->SameQ[fields,expectedFields],
   "chirality"->If[kind==="quartic","not_applicable",B4String[chir]],
   "self_conjugate"->If[kind==="quartic",sc,Null],
   "sample_Coupling_inputform"->B4String[expr],
   "sample_CouplingBar_inputform"->B4String[exprbar],
   "sample_Coupling_over_Matrix_inputform"->B4String[ratio],
   "sample_CouplingBar_over_Matrix_inputform"->B4String[ratiobar],
   "expected_internal_ratio_inputform"->B4String[expected],
   "native_exact_ratio_and_metadata_verified"->ok,
   "registered_invariant_inputform"->B4String[If[Head[inv]===Missing,inv,
      If[kind==="quartic",inv[b4h,b4h2,b4s,b4t],
       If[kind==="mass",inv[b4a,b4b],inv[b4a,b4b,b4c]]]]],
   "unique_arrangements_count"->If[ListQ[unique],Length[unique],Null],
   "physical_Lagrangian_coefficient_verified"->False|>
];
B4Family[name_String]:=Module[{d1,d2,df,alpha,yf,ok,second,rows},
 {d1,d2,df,alpha}=$B4Models[name];yf=(alpha+1)/2;
 ok=CheckAbort[Check[
  ResetModel[];T3RGBetaAddSM[];T3RGBetaAddHeavyFermion[df,yf];
  AddScalar[S1,GaugeRep->T3RGBetaGaugeRep[d1,alpha/2]];
  AddScalar[S2,GaugeRep->T3RGBetaGaugeRep[d2,(alpha+2)/2]];
  T3RGBetaAddT3Yukawas[d1,d2,df,yf];
  T3RGBetaAddMasses[d1,d2,df,yf];
  T3RGBetaAddBasicQuartics[d1,d2];
  T3RGBetaAddMixingQuartic[d1,d2];True,$Failed],$Aborted];
 If[ok=!=True,Return[<|"family"->name,"status"->"BlockedRegistration"|>]];
 second=If[yf===0,F,Fc];
 rows={
  B4CouplingRow[$yukawas,y1,"yukawa",Sqrt[2],{Bar[S1],l,F}],
  B4CouplingRow[$yukawas,y2,"yukawa",Sqrt[2],{S2,l,second}],
  B4CouplingRow[$fermionMasses,MF,"mass",2,{F,second}],
  B4CouplingRow[$quartics,lambdaT3,"quartic",6,{H,H,S1,Bar[S2]}],
  B4CouplingRow[$quartics,lambdaH1,"quartic",6,{Bar[H],H,Bar[S1],S1}],
  B4CouplingRow[$quartics,lambdaH2,"quartic",6,{Bar[H],H,Bar[S2],S2}]
 };
 <|"family"->name,"dimensions_alpha"->{d1,d2,df,alpha},
   "majorana"->(yf===0),"status"->If[
      AllTrue[rows,TrueQ[Lookup[#,"native_exact_ratio_and_metadata_verified",False]]&],
      "AllSixRegistryRatiosMeasured","OneOrMoreRegistryRatiosBlocked"],
   "operators"->rows,"physical_action_certified"->False|>
];
$B4GroupInv=Names["*`GroupInvBar"];
$B4Probe=CheckAbort[Check[
 <|"format"->"b2-rgbeta-internal-normalisation-native-v1",
   "status"->"NativeInternalTensorPrefactorsCapturedPhysicalActionPending",
   "git_head_at_capture"->$B4HeadFromPython,
   "git_source_blobs_at_capture"->$B4BlobsFromPython,
   "git_provenance_capture_method"->"Python wrapper: tree, index, clean-filter source guard",
   "rgbeta_package_file"->B4String[FindFile["RGBeta`"]],
   "GroupInvBar_symbols"->$B4GroupInv,
   "GroupInvBar_symbolic_probe"->If[Length[$B4GroupInv]===1,
     Module[{fn,eps0,s20,cp0,outE,outS,outC},
       fn=ToExpression[First[$B4GroupInv]];
       eps0=eps[SU2L[fund],b4i,b4j];
       s20=delS2[SU2L,b4a,b4i,b4j];
       cp0=I*eps0;
       outE=fn[eps0];outS=fn[s20];outC=fn[cp0];
       <|"eps_unchanged"->TrueQ[outE===eps0],
         "delS2_unchanged"->TrueQ[outS===s20],
         "i_eps_conjugates"->TrueQ[FullSimplify[outC==-I*eps0]],
         "eps_conjugation_inputform"->B4String[outE],
         "delS2_conjugation_inputform"->B4String[outS],
         "complex_epsilon_conjugation_inputform"->B4String[outC],
         "physical_component_orientation_inferred"->False|>],
       <|"status"->"GroupInvBarInternalSymbolNotUniquelyLocated",
         "physical_component_orientation_inferred"->False|>],
   "GroupInvBar_downvalues_inputform"->
      B4GroupInvDownValuesString[$B4GroupInv],
   "classes"->(B4Family/@Keys[$B4Models]),
   "internal_prefactor_is_physical_action_factor"->False,
   "numeric_component_axes_verified"->False,
   "physical_source_certified"->False|>,$Failed],$Aborted];
If[!AssociationQ[$B4Probe],Print["BLOCKED: native evaluation failed"];Exit[5]];

(* RawJSON cannot export symbolic Wolfram atoms such as Missing[...].
   Record exact offending JSON paths, not just a generic Export:: error.
   Never silently elevate an incomplete registry to verified native evidence. *)
B4JSONLeafQ[x_] := StringQ[x] || IntegerQ[x] ||
  (Head[x] === Real && NumberQ[x] && x =!= Indeterminate) ||
  x === True || x === False || x === Null;
B4JSONUnexpected[a_Association,p_String] := Flatten[
  KeyValueMap[Function[{k,v},
    Join[If[StringQ[k],{}, {<|"path"->p,"reason"->"NonStringKey",
       "inputform"->B4String[k]|>}],
      B4JSONUnexpected[v,p<>"."<>B4String[k]]]],a],1];
B4JSONUnexpected[a_List,p_String] := Flatten[
  MapIndexed[B4JSONUnexpected[#1,p<>"["<>ToString[First[#2]]<>"]"]&,a],1];
B4JSONUnexpected[x_,p_String] := If[B4JSONLeafQ[x],{},
  {<|"path"->p,"reason"->"NonJSONValue",
     "inputform"->B4String[x],"head"->B4String[Head[x]]|>}];
B4JSONSafe[a_Association] := Association@KeyValueMap[
  (If[StringQ[#1],#1,B4String[#1]] -> B4JSONSafe[#2])&,a];
B4JSONSafe[a_List] := B4JSONSafe /@ a;
B4JSONSafe[x_] := If[B4JSONLeafQ[x],x,
  "B2_UNSERIALIZABLE_NATIVE_VALUE["<>B4String[x]<>"]"];

$B4JSONIssues=B4JSONUnexpected[$B4Probe,"$"];
If[Length[$B4JSONIssues]>0,
  Print["B2 diagnostic: ",Length[$B4JSONIssues],
        " non-JSON native value(s); refusing verified status."];
  Scan[(Print["B2 JSON issue: ",#["path"]," = ",
        StringTake[#["inputform"],UpTo[240]]])&,$B4JSONIssues];
  $B4Probe=Join[$B4Probe,
    <|"status"->"NativeInternalTensorDataNotJSONSafe",
      "json_preflight_issues"->$B4JSONIssues|>]
];
$B4Safe=B4JSONSafe[$B4Probe];

(* Encode in memory before touching output. An ExportString failure cannot
   leave a zero-byte .json. All values are now explicitly JSON primitive types. *)
$B4Encoded=CheckAbort[
  Check[ExportString[$B4Safe,"RawJSON"],$Failed],$Aborted];
If[!StringQ[$B4Encoded] || StringLength[$B4Encoded]<16,
  Print["BLOCKED: RawJSON in-memory encoding failed; no output created"];
  Exit[6]];
$B4RoundTrip=CheckAbort[Check[ImportString[$B4Encoded,"RawJSON"],$Failed],$Aborted];
If[!AssociationQ[$B4RoundTrip] ||
   Lookup[$B4RoundTrip,"format",None]=!="b2-rgbeta-internal-normalisation-native-v1",
   Print["BLOCKED: RawJSON in-memory round trip failed; no output created"];
   Exit[6]];
If[!DirectoryQ[DirectoryName[$B4Out]],
 CreateDirectory[DirectoryName[$B4Out],CreateIntermediateDirectories->True]];
$B4Tmp=$B4Out<>".partial";
If[FileExistsQ[$B4Tmp],
  Print["BLOCKED: leftover partial file present: ",$B4Tmp];Exit[6]];
$B4Wrote=CheckAbort[
 Check[Export[$B4Tmp,$B4Encoded,"Text",CharacterEncoding->"UTF-8"],$Failed],
 $Aborted];
If[$B4Wrote===$Failed || $B4Wrote===$Aborted ||
   !FileExistsQ[$B4Tmp] || FileByteCount[$B4Tmp]<16,
 Print["BLOCKED: writing temporary JSON text failed: ",$B4Tmp];Exit[6]];
$B4OnDisk=CheckAbort[Check[Import[$B4Tmp,"RawJSON"],$Failed],$Aborted];
If[!AssociationQ[$B4OnDisk] ||
   Lookup[$B4OnDisk,"format",None]=!="b2-rgbeta-internal-normalisation-native-v1",
 Print["BLOCKED: temporary native JSON failed readback: ",$B4Tmp];Exit[6]];
$B4Renamed=CheckAbort[Check[RenameFile[$B4Tmp,$B4Out],$Failed],$Aborted];
If[$B4Renamed===$Failed || $B4Renamed===$Aborted ||
   !FileExistsQ[$B4Out],
 Print["BLOCKED: could not finalize native JSON: ",$B4Tmp];Exit[6]];
If[Length[$B4JSONIssues]>0,
 Print["BLOCKED: diagnostic native JSON written with unresolved values: ",$B4Out];
 Exit[7]];
If[!AllTrue[$B4Probe["classes"],#["status"]==="AllSixRegistryRatiosMeasured"&],
 Print["BLOCKED: native ratios differ; inspect ",$B4Out];Exit[7]];
Print["B2 native internal ratios 5/5, 30/30; physical=0/16: ",$B4Out];
Exit[0];

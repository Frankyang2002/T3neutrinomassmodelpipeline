(* B2: fresh-kernel, READ-ONLY RGBeta extended registry evidence.
   Does NOT claim numerical invariant axes, physical Yukawa Feynman rules,
   scalar quartic action factors, a basis match or physical certification.
   Usage: wolframscript -file scripts/RunPhysicalB2RGBetaRegistryAudit.wl OUTPUT.json
   OUTPUT must not exist. Opt-in general extensions must be off. *)
ClearAll["Global`*"];
If[Length[Rest[$ScriptCommandLine]] =!= 1,
  Print["Usage: OUTPUT.json"]; Exit[2]];
$B2Out=ExpandFileName[Last[$ScriptCommandLine]];
If[FileExistsQ[$B2Out],Print["Refusing existing output: ",$B2Out];Exit[2]];
$B2EnvNames={"T3_EXPERIMENTAL_GENERAL_FERMIONS",
 "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
 "T3_EXPERIMENTAL_GENERAL_QUARTICS",
 "T3_EXPERIMENTAL_FLAVOR_MASSES"};
If[AnyTrue[$B2EnvNames,Environment[#]==="1"&],
 Print["BLOCKED: experimental unrestricted extensions are on"]; Exit[3]];
$B2Root=DirectoryName@DirectoryName@ExpandFileName[$InputFileName];
If[CheckAbort[Check[Needs["RGBeta`"];True,$Failed],$Aborted]=!=True,
 Print["BLOCKED RGBeta load"];Exit[4]];
Get[FileNameJoin[{$B2Root,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
$B2Models=<|"A"->{1,3,2,0},"B"->{2,2,1,-1},
 "C"->{2,2,3,-1},"D"->{3,1,2,0},"E"->{3,3,2,0}|>;
ClearAll[B2KeyName,B2Read,B2String,B2RegistryRow,B2Run];
B2KeyName[x_]:=Last@StringSplit[ToString[x,InputForm],"`"];
B2Read[rec_,name_String]:=Module[{keys},
 If[!AssociationQ[rec],Return[Missing["RecordNotAssociation",name]]];
 keys=Select[Keys[rec],B2KeyName[#]===name&];
 If[Length[keys]!=1,Return[Missing["KeyMissingOrAmbiguous",name]]];
 Lookup[rec,First[keys]]];
B2String[x_]:=ToString[x,InputForm];
B2RegistryRow[registry_,symbol_Symbol,kind_String] := Module[
 {rec,keys,inv,fields,chirality,coup,coupbar},
 rec=If[AssociationQ[registry],Lookup[registry,symbol,Missing["MissingEntry"]],
    Missing["NoRegistry"]];
 keys=If[AssociationQ[rec],(ToString[#,InputForm]& /@ Keys[rec]),{}];
 inv=B2Read[rec,"Invariant"];fields=B2Read[rec,"Fields"];
 chirality=B2Read[rec,"Chirality"];
 coup=B2Read[rec,"Coupling"];coupbar=B2Read[rec,"CouplingBar"];
 <|"kind"->kind,"coupling"->SymbolName[symbol],
   "registry_keys"->keys,"fields_inputform"->B2String[fields],
   "invariant_function_inputform"->B2String[inv],
   "chirality_inputform"->B2String[chirality],
   "coupling_function_inputform"->B2String[coup],
   "coupling_bar_function_inputform"->B2String[coupbar],
   "has_invariant_record"->(Head[inv]=!=Missing),
   "has_fields_record"->(Head[fields]=!=Missing)|>
];
B2Run[family_String]:=Module[{d1,d2,df,alpha,yf,ok},
 {d1,d2,df,alpha}=$B2Models[family];yf=(alpha+1)/2;
 ok=CheckAbort[Check[
   ResetModel[]; T3RGBetaAddSM[]; T3RGBetaAddHeavyFermion[df,yf];
   AddScalar[S1,GaugeRep->T3RGBetaGaugeRep[d1,alpha/2]];
   AddScalar[S2,GaugeRep->T3RGBetaGaugeRep[d2,(alpha+2)/2]];
   T3RGBetaAddT3Yukawas[d1,d2,df,yf];
   T3RGBetaAddMasses[d1,d2,df,yf];
   T3RGBetaAddBasicQuartics[d1,d2];
   T3RGBetaAddMixingQuartic[d1,d2]; True,$Failed],$Aborted];
 If[ok=!=True,Return[<|"family"->family,"dimensions_alpha"->{d1,d2,df,alpha},
     "status"->"BlockedNativeRegistration"|>]];
 <|"family"->family,"dimensions_alpha"->{d1,d2,df,alpha},
   "status"->If[AssociationQ[$quartics] && KeyExistsQ[$quartics,lambdaT3],
       "RegistryRecordsCapturedNoActionCertificate",
       "PartialRegistryCaptureQuarticEntryNotLocated"],
   "y1"->B2RegistryRow[$yukawas,y1,"yukawa"],
   "y2"->B2RegistryRow[$yukawas,y2,"yukawa"],
   "MF"->B2RegistryRow[$fermionMasses,MF,"mass"],
   "lambdaT3"->B2RegistryRow[$quartics,lambdaT3,"quartic"],
   "quartic_registry_is_association"->AssociationQ[$quartics],
   "literal_mix_invariant_inputform"->B2String[T3RGBetaMixInvariant[d1,d2][i,j,s,t]],
   "physical_action_certified"->False|>
];
$B2Data=CheckAbort[Check[
 <|"format"->"b2-rgbeta-extended-native-registry-v1",
   "status"->"ExtendedFormalRegistryProbeCompletedActionPending",
   "rgbeta_eps_rules_inputform"->B2String[DownValues[eps]],
   "rgbeta_del_rules_inputform"->B2String[DownValues[del]],
   "rgbeta_delS2_rules_inputform"->B2String[DownValues[delS2]],
   "five_classes"->(B2Run/@Keys[$B2Models]),
   "missing_registry_entries_are_a_blocker_not_a_physical_witness"->True,
   "physical_source_certified"->False,
   "numeric_RGBeta_components_certified"->False,
   "full_action_feynman_rules_verified"->False,
   "comment"->"Formal registry metadata and expressions only, no physical conversion"|>,
 $Failed],$Aborted];
If[!AssociationQ[$B2Data]||Length[$B2Data["five_classes"]]=!=5,
 Print["BLOCKED native capture"];Exit[5]];
If[!DirectoryQ[DirectoryName[$B2Out]],
 CreateDirectory[DirectoryName[$B2Out],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[$B2Out,$B2Data,"RawJSON"],$Failed]===$Failed,
 Print["BLOCKED: cannot export JSON"];Exit[6]];
Print["B2 extended native formal registry exported (no physical certificate): ",$B2Out];
Exit[0];

(* WS4 read-only RGBeta native Yukawa+mass registry, single fresh kernel.
   Usage: wolframscript -file scripts/RunLODirectUVNativeFermionRegistry.wl A_p0|B_m1|C_m1 OUTPUT.json
   No RGEs/matching/production writes. No physical 4-spinor-to-Weyl claim. *)
ClearAll["Global`*"];
$wsArgs = Rest[$ScriptCommandLine];
If[Length[$wsArgs] =!= 2, Print["Usage: MODEL OUTPUT.json"]; Exit[2]];
$wsModel = $wsArgs[[1]];
$wsOutput = ExpandFileName[$wsArgs[[2]]];
$wsSpecs = <|"A_p0" -> {1,3,2,0}, "B_m1" -> {2,2,1,-1},
  "C_m1" -> {2,2,3,-1}|>;
If[!KeyExistsQ[$wsSpecs,$wsModel],Print["Unsupported model: ",$wsModel];Exit[2]];
If[FileExistsQ[$wsOutput],Print["Refusing to overwrite: ",$wsOutput];Exit[2]];
$wsRoot = DirectoryName@DirectoryName@ExpandFileName[$InputFileName];
$wsStep = "load-RGBeta";
If[CheckAbort[Check[Needs["RGBeta`"];True,$Failed],$Aborted] =!= True,
 Print["RGBeta unavailable"];Exit[3]];
Get[FileNameJoin[{$wsRoot,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
ClearAll[WSFail, WSKey, WSRead, WSYukawa, WSRun];
WSFail[reason_String] := Throw[<|"status" -> "Blocked", "step" -> $wsStep,
 "reason" -> reason|>, "WS4Stop"];
WSKey[k_] := Last@StringSplit[ToString[k,InputForm],"`"];
WSRead[a_Association,n_String] := Module[{keys},
 keys=Select[Keys[a],WSKey[#]===n&];
 If[Length[keys]=!=1, WSFail["Missing or duplicate registry field: "<>n]];
 Lookup[a,First[keys]]];
WSYukawa[y_Symbol,expected_List] := Module[
 {record,fields,chirality,f,fb,e,eb,t,tb,ratio,ratioBar},
 record=Lookup[$yukawas,y,Missing["NoYukawa"]];
 If[!AssociationQ[record],WSFail["Missing Yukawa: "<>SymbolName[y]]];
 fields=WSRead[record,"Fields"];chirality=WSRead[record,"Chirality"];
 If[!SameQ[fields,expected]||chirality=!=Right,
 WSFail["Unexpected Yukawa fields/chirality: "<>SymbolName[y]]];
 f=WSRead[record,"Coupling"]; fb=WSRead[record,"CouplingBar"];
 e=Quiet@Check[f[wsp,wsr],$Failed];eb=Quiet@Check[fb[wsp,wsr],$Failed];
 If[MemberQ[{$Failed,$Aborted},e]||MemberQ[{$Failed,$Aborted},eb],
 WSFail["Yukawa function evaluation failed: "<>SymbolName[y]]];
 t=Cases[e,HoldPattern[_[wsp,wsr]],Infinity];
 tb=Cases[eb,HoldPattern[_[wsp,wsr]],Infinity];
 If[Length[t]!=1||Length[tb]!=1,WSFail["Unexpected Yukawa matrix term count"]];
 If[FreeQ[First[t],Bar[y]]||FreeQ[First[tb],y],
 WSFail["Yukawa conjugation changed"]];
 ratio=FullSimplify[e/First[t]]; ratioBar=FullSimplify[eb/First[tb]];
 If[!TrueQ[FullSimplify[ratio==Sqrt[2]]]||
    !TrueQ[FullSimplify[ratioBar==Sqrt[2]]],
    WSFail["Unexpected Yukawa internal prefactor"]];
 <|"name"->SymbolName[y],"fields_inputform"->ToString[fields,InputForm],
 "chirality"->"Right","coupling_internal_factor"->ToString[ratio,InputForm],
 "couplingbar_internal_factor"->ToString[ratioBar,InputForm],
 "physical_action_coefficient_certified"->False|>
];
WSRun[] := Module[{spec,d1,d2,df,alpha,yf,setup,majorana,f,fc,m,mfields,mchirality,
 inv,invexpr,invkind,mfFn,mfExpr,mfTerms,mfRatio,rows},
 $wsStep="register-model";
 spec=$wsSpecs[$wsModel];{d1,d2,df,alpha}=spec;
 yf=(alpha+1)/2;majorana=(alpha===-1&&OddQ[df]);
 setup=CheckAbort[Check[
 ResetModel[];T3RGBetaAddSM[];
 T3RGBetaAddHeavyFermion[df,yf];
 AddScalar[S1,GaugeRep->T3RGBetaGaugeRep[d1,alpha/2]];
 AddScalar[S2,GaugeRep->T3RGBetaGaugeRep[d2,(alpha+2)/2]];
 T3RGBetaAddT3Yukawas[d1,d2,df,yf];
 T3RGBetaAddMasses[d1,d2,df,yf];True,$Failed],$Aborted];
 If[setup=!=True,WSFail["Native RGBeta registration failed"]];
 $wsStep="check-simultaneous-Yukawas";
 rows={WSYukawa[y1,{Bar[S1],l,F}],
    WSYukawa[y2,{S2,l,If[majorana,F,Fc]}]};
 $wsStep="read-heavy-mass";
 f=Lookup[$fermions,F,Missing["F"]];
 If[!AssociationQ[f],WSFail["F registry missing"]];
 fc=If[majorana,Null,Lookup[$fermions,Fc,Missing["Fc"]]];
 If[!majorana&&!AssociationQ[fc],WSFail["Fc registry missing"]];
 m=Lookup[$fermionMasses,MF,Missing["MF"]];
 If[!AssociationQ[m],WSFail["MF registry missing"]];
 mfields=WSRead[m,"Fields"];mchirality=WSRead[m,"Chirality"];
 If[!SameQ[mfields,If[majorana,{F,F},{F,Fc}]]||mchirality=!=Right,
 WSFail["Unexpected MF field pair/chirality"]];
 inv=WSRead[m,"Invariant"];invexpr=If[df===1,inv[],inv[wsma,wsmb]];
 invkind=If[df===1,"singlet",Which[
   SameQ[invexpr,eps[T3RGBetaRep[df],wsma,wsmb]],"epsilon",
   SameQ[invexpr,del[T3RGBetaRep[df],wsma,wsmb]],"delta",
   True,"unknown"]];
 If[invkind=!=If[df===1,"singlet",If[df===2,"epsilon","delta"]],
 WSFail["Unexpected mass invariant (Dirac doublet must use epsilon)"]];
 mfFn=WSRead[m,"Coupling"];
 mfExpr=Quiet@Check[mfFn[wsmr,wsms],$Failed];
 If[MemberQ[{$Failed,$Aborted},mfExpr],WSFail["MF coupling function failed"]];
 mfTerms=Cases[mfExpr,HoldPattern[_[wsmr,wsms]],Infinity];
 If[Length[mfTerms]!=1||(FreeQ[First[mfTerms],Bar[MF]]&&
    FreeQ[First[mfTerms],MF]),
 WSFail["MF internal coupling matrix orientation unexpected"]];
 mfRatio=FullSimplify[mfExpr/First[mfTerms]];
 If[!TrueQ[NumericQ[mfRatio]],WSFail["MF internal factor not numeric"]];
 $wsStep="complete";
 <|"format"->"t3-ws4-simultaneous-native-rgbeta-yukawa-mass-v1",
 "status"->"NativeRGBetaYukawaAndMassRegistryVerified", "model"->$wsModel,
 "dimensions_alpha"->spec,"majorana"->majorana,"yukawas"->rows,
 "mass"-><|"fields_inputform"->ToString[mfields,InputForm],
 "chirality"->"Right","invariant_kind"->invkind,
 "invariant_inputform"->ToString[invexpr,InputForm],
 "mass_internal_tensor_factor"->ToString[mfRatio,InputForm],
 "mass_matrix_barred_in_internal_coupling"->(!FreeQ[First[mfTerms],Bar[MF]]),
 "factor_is_physical_Lagrangian_coefficient"->False|>,
 "four_spinor_to_weyl_phase_measured"->False,
 "native_su2_numeric_component_axes_measured"->False,
 "physical_Matchete_RGBeta_dictionary_certified"->False,
 "rge_or_matching_performed"->False|>
];
$wsResult=CheckAbort[Catch[WSRun[],"WS4Stop"],$Aborted];
If[!AssociationQ[$wsResult],$wsResult=<|"status"->"Blocked","step"->$wsStep,
 "reason"->"Native evaluation aborted/non-association"|>];
If[!DirectoryQ[DirectoryName[$wsOutput]],
 CreateDirectory[DirectoryName[$wsOutput],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[$wsOutput,$wsResult,"RawJSON"],$Failed]===$Failed,
 Print["Cannot export JSON: ",$wsOutput];Exit[7]];
Print[$wsModel,": ",Lookup[$wsResult,"status","Blocked"]];
If[KeyExistsQ[$wsResult,"reason"],Print["Reason: ",$wsResult["reason"]]];
If[Lookup[$wsResult,"status",""]=!="NativeRGBetaYukawaAndMassRegistryVerified",Exit[6]];
Exit[0];

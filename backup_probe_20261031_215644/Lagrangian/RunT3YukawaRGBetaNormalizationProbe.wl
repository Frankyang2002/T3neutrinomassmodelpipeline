(* Workstream A: isolated, read-only RGBeta native Yukawa convention probe.
   Construct only the gauge/field/Yukawa sector in a fresh kernel. The full
   T3RGBetaBuild also registers heavy masses and quartics; those have no
   bearing on AddYukawa normalization and can fail independently.
   Does not edit production files.
   Command: wolframscript -file Lagrangian/RunT3YukawaRGBetaNormalizationProbe.wl OUTPUT.json A m2
   Supported models: A_m2, A_p0, B_m1, C_m1, D_m2, D_p0, E_m2, E_p0.
   Checks the native AddYukawa Coupling and CouplingBar functions, NOT the
   physical Matchete/RGBeta action normalization or RGBeta Ward identity.
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] < 3,
  Print["Usage: OUTPUT.json FAMILY m2|m1|p0"]; Exit[2]];
$NArgs = Take[$ScriptCommandLine, -3];
$NOut = ExpandFileName[$NArgs[[1]]];
$NModel = ToUpperCase[$NArgs[[2]]] <> "_" <> $NArgs[[3]];
$NSpecs = <|
 "A_m2"->{1,3,2,-2}, "A_p0"->{1,3,2,0},
 "B_m1"->{2,2,1,-1}, "C_m1"->{2,2,3,-1},
 "D_m2"->{3,1,2,-2}, "D_p0"->{3,1,2,0},
 "E_m2"->{3,3,2,-2}, "E_p0"->{3,3,2,0}|>;
If[!KeyExistsQ[$NSpecs,$NModel],Print["Unsupported model: ",$NModel];Exit[2]];
$NRoot = DirectoryName@DirectoryName@ExpandFileName[$InputFileName];
If[CheckAbort[Check[Needs["RGBeta`"];True,$Failed],$Aborted]=!=True,
  Print["RGBeta load failed"];Exit[3]];
Get[FileNameJoin[{$NRoot,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
ClearAll[NFail,NKey,NRead,NOne,NRun];
NFail[why_String]:=Throw[<|"status"->"NativeYukawaNormalizationProbeBlocked",
   "reason"->why,"failed_step"->$NStep,
   "physical_action_equivalence_verified"->False|>,"NStop"];
$NStep="initialization";
NKey[k_] := Last@StringSplit[ToString[k,InputForm],"`"];
NRead[rec_Association,key_String]:=Module[{ks},
  ks=Select[Keys[rec],NKey[#]===key&];
  If[Length[ks]=!=1,NFail["Missing or nonunique native registry key "<>key]];
  Lookup[rec,First[ks]]];
NOne[sym_Symbol,expectedFields_List] := Module[
 {info,fields,orientation,couplingFn,conjugateFn,e1,e2,terms1,terms2,den1,den2,ratio1,ratio2},
 info=Lookup[$yukawas,sym,Missing["NoYukawa"]];
 If[!AssociationQ[info],NFail["Missing Yukawa registry record "<>SymbolName[sym]]];
 fields=NRead[info,"Fields"];
 If[!SameQ[fields,expectedFields],
   NFail[SymbolName[sym]<>" field orientation changed: "<>ToString[fields,InputForm]]];
 orientation=NRead[info,"Chirality"];
 If[orientation=!=Right,NFail["Unexpected chirality for "<>SymbolName[sym]]];
 couplingFn=NRead[info,"Coupling"];
 conjugateFn=NRead[info,"CouplingBar"];
 e1=CheckAbort[Quiet@Check[couplingFn[t3p,t3r],$Failed],$Aborted];
 e2=CheckAbort[Quiet@Check[conjugateFn[t3p,t3r],$Failed],$Aborted];
 If[MemberQ[{$Failed,$Aborted},e1]||MemberQ[{$Failed,$Aborted},e2],
   NFail["Native coupling functions failed to evaluate"]];
 terms1=Cases[e1, HoldPattern[_[t3p,t3r]],Infinity];
 terms2=Cases[e2, HoldPattern[_[t3p,t3r]],Infinity];
 If[Length[terms1]=!=1 || Length[terms2]=!=1,
   NFail["Native coupling contains unexpected number of matrix coefficients"]];
 den1=First[terms1]; den2=First[terms2];
 If[FreeQ[den1,Bar[sym]] || FreeQ[den2,sym],
   NFail["Native Coupling/CouplingBar did not use expected symbol"]];
 ratio1=Quiet@FullSimplify[e1/den1];
 ratio2=Quiet@FullSimplify[e2/den2];
 If[!TrueQ[FullSimplify[ratio1==Sqrt[2]]] ||
    !TrueQ[FullSimplify[ratio2==Sqrt[2]]],
   NFail["Unexpected native Coupling/CouplingBar normalization for "<>SymbolName[sym]<>
     ": "<>ToString[{ratio1,ratio2},InputForm]]];
 <|"coupling"->SymbolName[sym],
   "fields_inputform"->ToString[fields,InputForm],
   "chirality"->"Right",
   "coefficient_of_Matrix_Bar_y"->ToString[ratio1,InputForm],
   "coefficient_of_Matrix_y_in_CouplingBar"->ToString[ratio2,InputForm],
   "native_complex_scalar_yukawa_tensor_prefactor_verified"->True|>
];
NRun[]:=Module[{spec,d1,d2,df,alpha,yF,setup,extra,inv,reg,baseline,expected,extraRow=None},
 $NStep="resolve-model";
 spec=$NSpecs[$NModel];{d1,d2,df,alpha}=spec;
 yF=(alpha+1)/2;
 (* Source-identical field and baseline Yukawa definitions, but deliberately
    exclude T3RGBetaAddMasses and scalar quartics. The RGBeta AddYukawa
    normalization depends on the registered scalar being complex and on
    Chirality -> Right, not on these unrelated interactions. *)
 $NStep="build-isolated-yukawa-sector";
 Print["STEP: ",$NStep," for ",$NModel];
 setup=CheckAbort[Check[
   ResetModel[];
   T3RGBetaAddSM[];
   T3RGBetaAddHeavyFermion[df,yF];
   AddScalar[S1,GaugeRep->T3RGBetaGaugeRep[d1,alpha/2]];
   AddScalar[S2,GaugeRep->T3RGBetaGaugeRep[d2,(alpha+2)/2]];
   T3RGBetaAddT3Yukawas[d1,d2,df,yF];
   True,
   $Failed],$Aborted];
 If[setup=!=True,NFail["Isolated native gauge/field/Yukawa setup failed"]];
 $NStep="verify-baseline-yukawa-registry";
 If[Length[Keys[$yukawas]]=!=5,NFail["Expected five baseline Yukawa records, found "<>
   ToString[Length[Keys[$yukawas]]]]];
 $NStep="evaluate-baseline-native-coupling-functions";
 baseline={
   NOne[y1,{Bar[S1],l,F}],
   NOne[y2,{S2,l,If[alpha===-1,F,Fc]}]
 };
 If[MemberQ[{"A_m2","A_p0","D_m2","D_p0","E_m2","E_p0"},$NModel],
   If[alpha===0,
     extra=y1Conj; expected={S1,l,F}; inv=T3RGBetaYukawaInvariant[d1,df],
     extra=y2Conj; expected={Bar[S2],l,Fc}; inv=T3RGBetaYukawaInvariant[d2,df]];
   (* Intentional: omit native CheckInvariance. It is unresolved for this
      installed RGBeta version; this probe checks registration metadata only. *)
   $NStep="register-optional-extra-yukawa";
   Print["STEP: ",$NStep," for ",$NModel];
   reg=CheckAbort[Check[
     AddYukawa[extra,expected,GroupInvariant->inv,
       CouplingIndices->({gen[#2],heavy[#3]}&),Chirality->Right];True,
     $Failed],$Aborted];
   If[reg=!=True,NFail["Extra AddYukawa registration emitted a message/failed"]];
   If[Length[Keys[$yukawas]]=!=6,NFail["Expected six registered Yukawas"]];
   $NStep="evaluate-extra-native-coupling-functions";
   extraRow=NOne[extra,expected];
 ];
 $NStep="complete";
 <|"status"->"NativeRGBetaYukawaInternalPrefactorsVerified",
  "model"->$NModel,"dimensions_alpha"->spec,
  "baseline_yukawas"->baseline,"extra_yukawa"->extraRow,
  "internal_real_scalar_tensor_prefactor"->"Sqrt[2]",
  "registration_scope"->"IsolatedGaugeFieldsAndYukawas",
  "masses_and_quartics_registered"->False,
  "source_convention"->"Coupling=2/sqrt(2)*Matrix[Bar[y]] for complex scalar, Chirality Right",
  "full_action_normalization_from_internal_prefactor_only"->False,
  "rgbeta_native_ward_invoked"->False,
  "matchete_rgbeta_physical_component_basis_certified"->False,
  "physical_action_equivalence_verified"->False,
  "production_modified"->False|>
];
$NResult=CheckAbort[Catch[NRun[],"NStop"],$Aborted];
If[!AssociationQ[$NResult],$NResult=<|"status"->"NativeYukawaNormalizationProbeBlocked",
 "reason"->"Nonassociation or aborted native result",
 "failed_step"->$NStep,
 "physical_action_equivalence_verified"->False|>];
$NResult=Join[<|"format"->"t3-rgbeta-native-yukawa-prefactor-v1",
 "model"->$NModel|>,$NResult];
If[!DirectoryQ[DirectoryName[$NOut]],
 CreateDirectory[DirectoryName[$NOut],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[$NOut,$NResult,"RawJSON"],$Failed]===$Failed,
 Print["Cannot export report"];Exit[7]];
Print[$NModel,": ",Lookup[$NResult,"status","Missing"]];
If[Lookup[$NResult,"status",""]=!="NativeRGBetaYukawaInternalPrefactorsVerified",
  Print["FAILED STEP: ",Lookup[$NResult,"failed_step","not reported"]];
  Print["REASON: ",Lookup[$NResult,"reason","not reported"]]];
If[Lookup[$NResult,"status",""]=!="NativeRGBetaYukawaInternalPrefactorsVerified",Exit[6]];
Exit[0];

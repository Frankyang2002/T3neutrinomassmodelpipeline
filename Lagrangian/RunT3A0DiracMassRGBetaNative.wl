(* Workstream A — read-only native RGBeta A(0) heavy Dirac mass gauge-tensor probe.
   Usage: wolframscript -file Lagrangian/RunT3A0DiracMassRGBetaNative.wl OUTPUT.json
   Executes the actual T3RGBetaAddHeavyFermion and T3RGBetaAddMasses
   model functions; reads, but never changes, their registrations.
   A successful PROBE does not mean the registered mass is gauge invariant.
*)
ClearAll["Global`*"];
$DMArgs = Rest[$ScriptCommandLine];
If[Length[$DMArgs] =!= 1, Print["Usage: OUTPUT.json"]; Exit[2]];
$DMOut = ExpandFileName[First[$DMArgs]];
$DMDir = DirectoryName@DirectoryName@ExpandFileName[$InputFileName];
$DMStep = "load-rgbeta";
If[CheckAbort[Check[Needs["RGBeta`"] ; True, $Failed], $Aborted] =!= True,
 Print["Could not load RGBeta"]; Exit[3]];
Get[FileNameJoin[{$DMDir,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
ClearAll[DMFail,DMKey,DMRead,DMRun];
DMFail[msg_String] := Throw[<|"status"->"NativeA0DiracMassProbeBlocked",
   "failed_step"->$DMStep,"reason"->msg|>,"DMStop"];
DMKey[k_] := Last@StringSplit[ToString[k,InputForm],"`"];
DMRead[rec_Association,name_String] := Module[{matches},
 matches=Select[Keys[rec],DMKey[#]===name&];
 If[Length[matches]=!=1,DMFail["Missing or nonunique registry key: "<>name]];
 Lookup[rec,First[matches]]
];
DMRun[] := Module[{setup,f,fc,m,fields,inv,expr,ch,repF,repFc,fun,coefficient,
   matrixTerms,matrixTerm,ratio,invariantKind},
 $DMStep="register-native-A_p0-fields-yukawas-and-masses";
 Print["STEP: ",$DMStep];
 setup=CheckAbort[Check[
   ResetModel[];
   T3RGBetaAddSM[];
   T3RGBetaAddHeavyFermion[2,1/2];
   AddScalar[S1,GaugeRep->T3RGBetaGaugeRep[1,0]];
   AddScalar[S2,GaugeRep->T3RGBetaGaugeRep[3,1]];
   T3RGBetaAddT3Yukawas[1,3,2,1/2];
   T3RGBetaAddMasses[1,3,2,1/2];
   True,
   $Failed],$Aborted];
 If[setup=!=True,DMFail["RGBeta A_p0 native setup emitted a message or failed"]];
 $DMStep="read-native-F-Fc-and-MF-registry";
 f=Lookup[$fermions,F,Missing["F"]];
 fc=Lookup[$fermions,Fc,Missing["Fc"]];
 m=Lookup[$fermionMasses,MF,Missing["MF"]];
 If[!And@@(AssociationQ/@{f,fc,m}),DMFail["F/Fc/MF registry entry missing"]];
 repF=DMRead[f,"GaugeRep"]; repFc=DMRead[fc,"GaugeRep"];
 fields=DMRead[m,"Fields"];
 ch=DMRead[m,"Chirality"];
 inv=DMRead[m,"Invariant"];
 expr=inv[dm1,dm2];
 If[!SameQ[repF,{U1Y[1/2],SU2L[fund]}] ||
    !SameQ[repFc,{U1Y[-1/2],SU2L[fund]}],
   DMFail["F/Fc native SU2 or U1 reps changed; current audit requires both SU2 fundamentals"]];
 If[!SameQ[fields,{F,Fc}] || ch=!=Right,
   DMFail["MF field order or chirality changed"]];
 invariantKind=Which[
   SameQ[expr,del[SU2L[fund],dm1,dm2]],"delta",
   SameQ[expr,eps[SU2L[fund],dm1,dm2]],"epsilon",
   True,"other"];
 $DMStep="read-native-MF-tensor-prefactor";
 fun=DMRead[m,"Coupling"];
 coefficient=CheckAbort[Check[fun[dmi,dmj],$Failed],$Aborted];
 If[coefficient===$Failed || coefficient===$Aborted,
  DMFail["MF tensor Coupling function failed to evaluate"]];
 matrixTerms=Cases[coefficient,HoldPattern[_[dmi,dmj]],Infinity];
 If[Length[matrixTerms]=!=1 || FreeQ[First[matrixTerms],Bar[MF]],
  DMFail["MF tensor does not contain one barred matrix coupling"]];
 matrixTerm=First[matrixTerms];
 ratio=Quiet@FullSimplify[coefficient/matrixTerm];
 If[!TrueQ[FullSimplify[ratio==2]],
  DMFail["MF tensor normalization is not 2"]];
 $DMStep="complete";
 <|"format"->"t3-workstream-a-native-a0-mass-gauge-tensor-v1",
  "status"->"NativeA0FandFcMassTensorCaptured",
  "model"->"A_p0",
  "F_gauge_rep_inputform"->ToString[repF,InputForm],
  "Fc_gauge_rep_inputform"->ToString[repFc,InputForm],
  "both_fields_in_SU2_fundamental"->True,
  "F_and_Fc_have_opposite_U1_hypercharge"->True,
  "MF_fields_inputform"->ToString[fields,InputForm],
  "MF_chirality"->"Right",
  "MF_group_invariant_inputform"->ToString[expr,InputForm],
  "MF_group_invariant_kind"->invariantKind,
  "MF_internal_mass_tensor_prefactor"->ToString[ratio,InputForm],
  "native_invariant_checked_by_rgbeta"->False,
  "component_ward_evaluated_by_native_package"->False,
  "mass_beta_evaluated"->False,
  "production_files_modified"->False|>
];
$DMRes=CheckAbort[Catch[DMRun[],"DMStop"],$Aborted];
If[!AssociationQ[$DMRes],$DMRes=<|"status"->"NativeA0DiracMassProbeBlocked",
 "failed_step"->$DMStep,"reason"->"Wolfram aborted or returned no record"|>];
If[!DirectoryQ[DirectoryName[$DMOut]],
 CreateDirectory[DirectoryName[$DMOut],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[$DMOut,$DMRes,"RawJSON"],$Failed]===$Failed,
 Print["Unable to write native JSON: ",$DMOut];Exit[7]];
Print["Native A(0) mass registry: ",Lookup[$DMRes,"status","Missing"]];
Print["Mass invariant: ",Lookup[$DMRes,"MF_group_invariant_inputform","Unavailable"]];
If[Lookup[$DMRes,"status",""]=!="NativeA0FandFcMassTensorCaptured",
 Print["Failed step: ",Lookup[$DMRes,"failed_step","?"]];
 Print["Reason: ",Lookup[$DMRes,"reason","?"]];Exit[6]];
Exit[0];

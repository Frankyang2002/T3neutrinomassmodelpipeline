(* Track A A2: NEW native Matchete Feynman-rule-level witness, three models.
   No Match[], RGE, EFT, source mutation, output overwrite, or invented C matrix.
   The native source mass/Yukawa coefficient tests were completed previously;
   this witness asks the additional, distinct FeynmanRules question.
   Usage: wolframscript -file Numerical/lo_certification/physical_normalisation/parallel_A_spinor/native_rule_probe.wl OUTPUT.json
*)
ClearAll["Global`*"];
taArgs=Rest[$ScriptCommandLine];
If[Length[taArgs]=!=1,Print["Usage: OUTPUT.json"];Exit[2]];
taOut=ExpandFileName[First[taArgs]];
If[FileExistsQ[taOut],Print["BLOCKED: existing output is immutable: ",taOut];Exit[2]];
taDir=DirectoryName[ExpandFileName[$InputFileName]];
taRoot=ExpandFileName[FileNameJoin[{taDir,"..","..","..",".."}]];
taFlags={"T3_EXPERIMENTAL_GENERAL_FERMIONS","T3_EXPERIMENTAL_FLAVOR_MASSES",
 "T3_EXPERIMENTAL_GENERAL_QUARTICS","T3_EXPERIMENTAL_GENERAL_QUADRATICS"};
taPrevious=AssociationMap[Function[n,With[{v=Environment[n]},If[StringQ[v],v,"UNSET"]]],taFlags];
SetEnvironment[Thread[taFlags->ConstantArray["0",Length[taFlags]]]];
If[!AllTrue[taFlags,Environment[#]==="0"&],Print["BLOCKED: restricted flags not zero"];Exit[3]];
(* Needs must execute before expressions mentioning Matchete-exported symbols are parsed. *)
taLoaded=CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]];True,$Failed],$Aborted];
If[taLoaded=!=True,Print["BLOCKED: Matchete not available"];Exit[4]];
Get[FileNameJoin[{taRoot,"Lagrangian","model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{taRoot,"Lagrangian","model","LagrangianBuilder.wl"}]];
ClearAll[TAFail,TAText,TATerms,TAPick,TARule,TARun,TAHashes];
TAFail[reason_String]:=Throw[<|"status"->"BLOCKED","reason"->reason|>,"TAStop"];
TAText[x_]:=ToString[x,InputForm];
TATerms[e_,strings_List]:=Select[
 If[Head[Expand[e]]===Plus,List@@Expand[e],{Expand[e]}],
 Function[t,With[{s=TAText[t]},AllTrue[strings,StringContainsQ[s,#]&]]]];
TAPick[e_,tokens_List]:=Module[{a=TATerms[e,tokens]},
 If[Length[a]=!=1,TAFail["Non-unique native action term: "<>TAText[tokens]]];First[a]];
TAHashes[]:=Module[{paths,p},
 paths={"Lagrangian/model/T3ModelCatalog.wl","Lagrangian/model/LagrangianBuilder.wl",
 "Lagrangian/model/T3Fields.wl","Lagrangian/model/SU2Invariants.wl",
 "Lagrangian/interactions/T3Topology.wl","Lagrangian/interactions/ScalarPotential.wl"};
 Association@Table[p=FileNameJoin[Prepend[FileNameSplit[q],taRoot]];
  If[!FileExistsQ[p],TAFail["Missing required production file: "<>q]];
  q->IntegerString[FileHash[p,"SHA256"],16,64],{q,paths}]];
(* v0.6 Matchete FeynmanRules: a genuinely separate rule derivation.
   If missing/undocumented/unsupported on the local installation, FAIL CLOSED.
   Rephasing of an isolated overall self-energy is NOT measured here. *)
TARule[expr_,name_String]:=Module[{result},
 result=CheckAbort[TimeConstrained[Check[UsingFrontEnd[FeynmanRules[expr]],$Failed],180,$Failed],$Aborted];
 If[MemberQ[{$Failed,$Aborted},result] || Head[result]===FeynmanRules,
  TAFail["FeynmanRules failed or remained unevaluated for "<>name]];
 TAText[result]];
TARun[]:=Module[{all=<||>,specs,entry,id,model,built,lint,free,y1,y1h,y2,y2h,m,
  hcy1,nonhermitian,combined,checks,rulePositive,ruleHostile,per},
 If[!MemberQ[Names["Matchete`*"],"Matchete`FeynmanRules"],
  TAFail["Installed Matchete has no public Matchete`FeynmanRules (requires 0.6 or compatible)"]];
 specs={{"A_p0","A",0},{"B_m1","B",-1},{"C_m1","C",-1}};
 Do[
  {id,model,entry}=spec;
  built=CheckAbort[Check[UsingFrontEnd[BuildT3Lagrangian[T3ModelFromClass[model,entry]]],$Failed],$Aborted];
  If[!AssociationQ[built] || Lookup[built,"Status",""]=!="Success" ||
    !TrueQ[Lookup[built,"FullValidation",False]],TAFail["Restricted native build failed: "<>id]];
  lint=built["LInt"];free=built["LFree"];
  y1=TAPick[lint,{"Coupling[y1,","CG[T3Y1CG,"}];
  y1h=TAPick[lint,{"Bar[Coupling[y1,","CG[Bar[T3Y1CG],"}];
  y2=TAPick[lint,{"Coupling[y2,","CG[T3Y2CG,"}];
  y2h=TAPick[lint,{"Bar[Coupling[y2,","CG[Bar[T3Y2CG],"}];
  m=TAPick[free,{"Coupling[MF,","NewFermion"}];
  combined=y1+y1h+y2+y2h+m;
  (* The bad control is NOT a rephasing: ONLY the y1 H.c. sign is reversed,
     so the action becomes non-Hermitian and native 3-point rules must differ. *)
  hcy1=y1+y1h;nonhermitian=y1-y1h;
  rulePositive=TARule[hcy1,"y1 + H.c. "<>id];
  ruleHostile=TARule[nonhermitian,"y1 - H.c. "<>id];
  checks=<|"native_feynman_rules_return_distinct_HC_sign"->
     (rulePositive=!=ruleHostile),
   "nonempty_HC_rules"->(StringLength[rulePositive]>2),
   "nonempty_hostile_rules"->(StringLength[ruleHostile]>2)|>;
  If[!AllTrue[Values[checks],TrueQ],TAFail["Negative control undetected for "<>id]];
  all[id]=<|
    "y1_hermitian_rule_inputform"->rulePositive,
    "y1_wrong_HC_sign_rule_inputform"->ruleHostile,
    "y2_hermitian_rule_inputform"->TARule[y2+y2h,"y2 + H.c. "<>id],
    "fermion_mass_rule_inputform"->TARule[m,"mass "<>id],
    "combined_yukawa_mass_rule_inputform"->TARule[combined,"combined "<>id],
    "negative_control"->checks,
    "relative_numeric_GammaCC_to_Weyl_phase_measured"->False|>;
  ,{spec,specs}];
 <|"format"->"t3-trackA-native-feynman-rule-probe-v1",
   "status"->"NativeMatcheteRulesCaptured_RelativeWeylPhaseNotEstablished",
   "source_sha256"->TAHashes[],"models"->all,
   "environment_before"->taPrevious,
   "environment_effective"->AssociationMap[Environment,taFlags],
   "native_feynman_rules_invoked"->True,
   "isolated_wrong_HC_sign_negative_control_detected"->True,
   "numerical_CG_epsilon_triplet_basis_map_measured_here"->False,
   "installed_relative_GammaCC_phase_certified"->False,
   "physical_source_certified"->False,
   "matching_or_running_performed"->False|>];
taRes=CheckAbort[Catch[TARun[],"TAStop"],$Aborted];
If[!AssociationQ[taRes]||Lookup[taRes,"status",""]=!=
   "NativeMatcheteRulesCaptured_RelativeWeylPhaseNotEstablished",
 Print["BLOCKED: ",If[AssociationQ[taRes],Lookup[taRes,"reason","No reason"],"evaluation aborted"]];Exit[5]];
If[!DirectoryQ[DirectoryName[taOut]],
 CreateDirectory[DirectoryName[taOut],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[taOut,taRes,"RawJSON"],$Failed]===$Failed,
 Print["BLOCKED: cannot export output"];Exit[6]];
Print["Three native rule-level records, including wrong-HC sign controls, exported: ",taOut];
Print["Installed relative numeric GammaCC/Weyl convention UNRESOLVED; certified 0/16."];
Exit[0];

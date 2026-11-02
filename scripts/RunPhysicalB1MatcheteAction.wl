(* B1 new, restricted, read-only. Fresh wolframscript process only.
   Usage: wolframscript -file scripts/RunPhysicalB1MatcheteAction.wl A_p0|B_m1|C_m1 OUTPUT.json
   Does not call Match, RGE, modify source, or claim a concrete GammaCC phase.
*)
ClearAll["Global`*"];
$b1args=Rest[$ScriptCommandLine];
If[Length[$b1args]=!=2,Print["Usage: MODEL OUTPUT.json"];Exit[2]];
$b1model=$b1args[[1]];$b1out=ExpandFileName[$b1args[[2]]];
$b1cases=<|"A_p0"->{"A",0},"B_m1"->{"B",-1},"C_m1"->{"C",-1}|>;
If[!KeyExistsQ[$b1cases,$b1model],Print["Unsupported model"];Exit[2]];
If[FileExistsQ[$b1out],Print["Refusing to overwrite native witness"];Exit[2]];
$b1flags={"T3_EXPERIMENTAL_GENERAL_FERMIONS",
 "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
 "T3_EXPERIMENTAL_GENERAL_QUARTICS",
 "T3_EXPERIMENTAL_FLAVOR_MASSES"};
$b1oldFlags=AssociationMap[Environment,$b1flags];
SetEnvironment[Thread[$b1flags->ConstantArray["0",Length[$b1flags]]]];
If[!AllTrue[$b1flags,Environment[#]==="0"&],Print["Flag isolation failed"];Exit[2]];
$b1root=DirectoryName@DirectoryName@ExpandFileName[$InputFileName];
(* Native Matchete must be loaded BEFORE parsing definitions containing Bar/NCM etc. *)
If[CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"] ;True],$Failed],$Aborted]=!=True,
 Print["Matchete unavailable"];Exit[3]];
Get[FileNameJoin[{$b1root,"Lagrangian","model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{$b1root,"Lagrangian","model","LagrangianBuilder.wl"}]];
ClearAll[B1Fail,B1TermQ,B1Run,B1Str];
B1Fail[s_String]:=Throw[<|"status"->"Blocked","reason"->s|>,"B1Stop"];
B1Str[x_]:=ToString[x,InputForm];
B1TermQ[t_,y_Symbol,barred_]:=If[TrueQ[barred],
  !FreeQ[t,HoldPattern[Bar[Coupling[y,___]]]],
  !FreeQ[t,HoldPattern[Coupling[y,___]]] &&
    FreeQ[t,HoldPattern[Bar[Coupling[y,___]]]]];
B1Run[]:=Module[{spec,model,built,lint,lf,terms,fterms,rows=<||>,orig,hc,
   massTerms, checks, c, h, nativeProjectors},
 spec=$b1cases[$b1model];model=T3ModelFromClass@@spec;
 If[!AssociationQ[model],B1Fail["Invalid model association"]];
 built=CheckAbort[Check[UsingFrontEnd[BuildT3Lagrangian[model]],$Failed],$Aborted];
 If[!AssociationQ[built]||Lookup[built,"Status",""]=!="Success"||
    !TrueQ[Lookup[built,"FullValidation",False]],B1Fail["Assembled native action failed validation"]];
 lint=built["LInt"];lf=built["LFree"];
 terms=List@@Expand[lint]; fterms=List@@Expand[lf];
 Do[
  orig=Select[terms,B1TermQ[#,y,False]&];hc=Select[terms,B1TermQ[#,y,True]&];
  If[Length[orig]<1||Length[hc]<1,
   B1Fail["Could not isolate original and HC for "<>SymbolName[y]]];
  rows[SymbolName[y]]=<|
   "original_term_count"->Length[orig],"hc_term_count"->Length[hc],
   "original_terms_inputform"->(B1Str/@orig),"hc_terms_inputform"->(B1Str/@hc),
   "original_projectors_inputform"->(B1Str/@Cases[orig,PL|PR,Infinity]),
   "hc_projectors_inputform"->(B1Str/@Cases[hc,PL|PR,Infinity]),
   "original_conjugation_field_inputform"->(B1Str/@Cases[orig,HoldPattern[CConj[_]],Infinity]),
   "symbolic_coupling_only_not_absolute_numerical_coefficient"->True|>;
 ,{y,{y1,y2}}];
 massTerms=Select[fterms,!FreeQ[#,MF]&];
 If[Length[massTerms]<1,B1Fail["Free action lacks MF mass term"]];
 checks=<|
   "action_success"->True,
   "both_yukawas_and_hc_isolated"->True,
   "mass_term_extracted"->True,
   "cg_names_registered_by_build"->MemberQ[Lookup[built,"AllowedInteractions",{}],"Yukawa1_GenericCG"]&&
       MemberQ[Lookup[built,"AllowedInteractions",{}],"Yukawa2_GenericCG"],
   "four_experimental_flags_zero"->AllTrue[$b1flags,Environment[#]==="0"&]|>;
 If[!AllTrue[Values[checks],TrueQ],B1Fail["Executable assembled-action contract failed: "<>B1Str[checks]]];
 nativeProjectors=Cases[{lint,lf},PL|PR,Infinity];
 <|"format"->"t3-b1-native-restricted-action-v1",
   "status"->"NativeActionTermsIsolated_PhysicalDictionaryPending",
   "model"->$b1model,"tests"->checks,
   "source_action_LFree_inputform"->B1Str[lf],
   "source_action_LInt_inputform"->B1Str[lint],
   "source_action_LUV_inputform"->B1Str[built["LUV"]],
   "mass_terms_inputform"->(B1Str/@massTerms),
   "mass_term_count"->Length[massTerms],
   "all_projectors_inputform"->(B1Str/@nativeProjectors),
   "yukawas"->rows,
   "registered_CG_symbol_names"->{"T3Y1CG","T3Y2CG"},
   "cg_component_coefficients_measured_here"->False,
   "source_gammaCC_concrete_matrix_measured"->False,
   "source_mass_flavour_matrix_general_measured"->False,
   "physical_source_certified"->False,
   "matching_or_running_performed"->False,
   "prior_process_flags"->$b1oldFlags,
   "effective_flags"->AssociationMap[Environment,$b1flags]|>
 ];
$b1res=CheckAbort[Catch[B1Run[],"B1Stop"],$Aborted];
If[!AssociationQ[$b1res],$b1res=<|"status"->"Blocked","reason"->"Aborted/non-association"|>];
If[!DirectoryQ[DirectoryName[$b1out]],CreateDirectory[DirectoryName[$b1out],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[$b1out,$b1res,"RawJSON"],$Failed]===$Failed,
 Print["Export failed"];Exit[7]];
Print[$b1model,": ",Lookup[$b1res,"status","Blocked"]];
If[Lookup[$b1res,"status",""]=!="NativeActionTermsIsolated_PhysicalDictionaryPending",Exit[6]];
Exit[0];

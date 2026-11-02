(* WS4 read-only Matchete fermion action witness for A_p0, B_m1, C_m1.
   Usage: wolframscript -file scripts/RunLODirectUVNativeFermionAction.wl MODEL OUTPUT.json
   Restricted opt-ins are forced off in child process. No matching/RGE. *)
ClearAll["Global`*"];
$faArgs=Rest[$ScriptCommandLine];
If[Length[$faArgs]=!=2,Print["Usage: MODEL OUTPUT.json"];Exit[2]];
$faModel=$faArgs[[1]];$faOut=ExpandFileName[$faArgs[[2]]];
$faSpecs=<|"A_p0"->{"A",0},"B_m1"->{"B",-1},"C_m1"->{"C",-1}|>;
If[!KeyExistsQ[$faSpecs,$faModel],Print["Unsupported model"];Exit[2]];
If[FileExistsQ[$faOut],Print["Refusing overwrite"];Exit[2]];
$faFlags={"T3_EXPERIMENTAL_GENERAL_FERMIONS",
 "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
 "T3_EXPERIMENTAL_GENERAL_QUARTICS",
 "T3_EXPERIMENTAL_FLAVOR_MASSES"};
$faPrior=AssociationMap[
 Function[name,With[{v=Environment[name]},If[StringQ[v],v,"UNSET"]]],
 $faFlags];
SetEnvironment[Thread[$faFlags->ConstantArray["0",Length[$faFlags]]]];
If[!AllTrue[$faFlags,Environment[#]==="0"&],Print["Flag isolation failed"];Exit[2]];
$faRoot=DirectoryName@DirectoryName@ExpandFileName[$InputFileName];
If[CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"];True],$Failed],$Aborted]=!=True,
 Print["Matchete unavailable"];Exit[3]];
Get[FileNameJoin[{$faRoot,"Lagrangian","model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{$faRoot,"Lagrangian","model","LagrangianBuilder.wl"}]];
ClearAll[FARun,FAFail];
FAFail[msg_String]:=Throw[<|"status"->"Blocked","reason"->msg|>,"FAStop"];
FARun[]:=Module[{spec,model,built,lint,free,lstr,fstr,spin,expected,checks},
 spec=$faSpecs[$faModel];model=T3ModelFromClass[spec[[1]],spec[[2]]];
 If[!AssociationQ[model],FAFail["Invalid model"]];
 built=CheckAbort[Check[UsingFrontEnd[BuildT3Lagrangian[model]],$Failed],$Aborted];
 If[!AssociationQ[built]||Lookup[built,"Status",""]=!="Success"||
   !TrueQ[Lookup[built,"FullValidation",False]],FAFail["Native Matchete build failed"]];
 lint=built["LInt"];free=built["LFree"];
 lstr=ToString[lint,InputForm];fstr=ToString[free,InputForm];
 (* Use the Matchete-expanded internal expressions, not unexpanded y1[] form. *)
 checks=<|"y1_present"->StringContainsQ[lstr,"Coupling[y1,"],
   "y2_present"->StringContainsQ[lstr,"Coupling[y2,"],
   "both_yukawa_hc_present"->(StringContainsQ[lstr,"Bar[Coupling[y1,"]&&
       StringContainsQ[lstr,"Bar[Coupling[y2,"]),
   "heavy_fermion_free_mass_symbol_present"->StringContainsQ[fstr,"MF"],
   "mass_and_yukawa_use_newfermion"->
     (StringContainsQ[lstr,"NewFermion"]&&StringContainsQ[fstr,"NewFermion"])
 |>;
 If[!AllTrue[Values[checks],TrueQ],FAFail["Assembled native action missing expected Yukawa/mass terms: "<>
    ToString[checks,InputForm]]];
 <|"format"->"t3-ws4-restricted-native-matchete-fermion-action-v1",
 "status"->"NativeRestrictedFermionActionCaptured","model"->$faModel,
 "native_lagrangian_validation"->True,"checks"->checks,
 "mass_free_lagrangian_inputform"->fstr,
 "yukawa_and_mixing_interactions_inputform"->lstr,
 "flag_inheritance"->$faPrior,
 "effective_experimental_flags"->AssociationMap[Environment,$faFlags],
 "physical_charge_conjugation_phase_measured"->False,
 "physical_Matchete_RGBeta_dictionary_certified"->False,
 "matching_or_rge_performed"->False|>
];
$faRes=CheckAbort[Catch[FARun[],"FAStop"],$Aborted];
If[!AssociationQ[$faRes],$faRes=<|"status"->"Blocked","reason"->"Aborted"|>];
If[!DirectoryQ[DirectoryName[$faOut]],
 CreateDirectory[DirectoryName[$faOut],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[$faOut,$faRes,"RawJSON"],$Failed]===$Failed,Exit[7]];
Print[$faModel,": ",Lookup[$faRes,"status","Blocked"]];
If[KeyExistsQ[$faRes,"reason"],Print[$faRes["reason"]]];
If[Lookup[$faRes,"status",""]=!="NativeRestrictedFermionActionCaptured",Exit[6]];
Exit[0];

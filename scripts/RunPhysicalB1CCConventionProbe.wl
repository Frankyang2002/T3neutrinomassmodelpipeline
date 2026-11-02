(* B1 standalone read-only native ABSTRACT CC probe. Fresh kernel only.
   wolframscript -file scripts/RunPhysicalB1CCConventionProbe.wl OUTPUT.json
   It does not equate an abstract GammaCC to a particular 4x4 numeric matrix.
   No matching, RGEs, external file writes, or production changes. *)
ClearAll["Global`*"];
$b1args=Rest[$ScriptCommandLine];
If[Length[$b1args]=!=1,Print["Usage: OUTPUT.json"];Exit[2]];
$b1out=ExpandFileName[First[$b1args]];
If[FileExistsQ[$b1out],Print["Refusing to overwrite: ",$b1out];Exit[2]];
$b1flags={"T3_EXPERIMENTAL_GENERAL_FERMIONS",
 "T3_EXPERIMENTAL_GENERAL_QUADRATICS","T3_EXPERIMENTAL_GENERAL_QUARTICS",
 "T3_EXPERIMENTAL_FLAVOR_MASSES"};
SetEnvironment[Thread[$b1flags->ConstantArray["0",Length[$b1flags]]]];
If[!AllTrue[$b1flags,Environment[#]==="0"&],Print["Flag isolation failed"];Exit[3]];
$b1root=DirectoryName@DirectoryName@ExpandFileName[$InputFileName];
(* Package symbols must be defined before lexical parsing of B1Run below. *)
If[CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"] ;True],$Failed],$Aborted]=!=True,
 Print["Native Matchete unavailable"];Exit[3]];
Get[FileNameJoin[{$b1root,"Lagrangian","model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{$b1root,"Lagrangian","model","LagrangianBuilder.wl"}]];
ClearAll[B1Stop,B1Run,B1String];
B1String[x_]:=ToString[x,InputForm];
B1Stop[msg_String]:=Throw[<|"status"->"Blocked","reason"->msg|>,"B1CCStop"];
B1Run[]:=Module[{model,built,f,lfield,cc,symbolic,checks,manifest},
 model=T3ModelFromClass["A",0];
 If[!AssociationQ[model],B1Stop["A_p0 model unavailable"]];
 built=CheckAbort[Check[UsingFrontEnd[BuildT3Lagrangian[model]],$Failed],$Aborted];
 If[!AssociationQ[built]||Lookup[built,"Status",""]=!="Success"||
    !TrueQ[Lookup[built,"FullValidation",False]],
    B1Stop["Native A_p0 Lagrangian not valid"]];
 f=NewFermion[j,r];lfield=l[i,p];
 cc=DiracProduct[GammaCC];
 checks=<|
  "native_cc_antisymmetry"->TrueQ[FullSimplify[Transp[GammaCC]==-GammaCC]],
  "native_cc_square_minus_identity"->TrueQ[FullSimplify[DiracProduct[GammaCC,GammaCC]==-1]],
  "native_cconj_field_definition"->TrueQ[FullSimplify[
    NCM[Bar[lfield],CConj[f]]==NCM[Bar[lfield],cc,Transp[Bar[f]]]]],
  "native_cconj_projector_flip"->TrueQ[FullSimplify[
    CConj[NCM[PL,f]]==NCM[PR,CConj[f]]]],
  "native_assembled_yukawa_both_present"->(!FreeQ[built["LInt"],y1]&&
                                          !FreeQ[built["LInt"],y2]),
  "native_restricted_flags"->AllTrue[$b1flags,Environment[#]==="0"&]|>;
 If[!AllTrue[Values[checks],TrueQ],B1Stop["Native abstract CC check failed: "<>B1String[checks]]];
 (* DownValues/OwnValues describe the installed symbolic object, not a proof
    of a sign relative to the externally chosen Weyl gamma representation. *)
 manifest=<|
  "symbol_full_name"->Context[GammaCC]<>SymbolName[GammaCC],
  "gamma_cc_ownvalues_count"->Length[OwnValues[GammaCC]],
  "gamma_cc_downvalues_count"->Length[DownValues[GammaCC]],
  "gamma_cc_upvalues_count"->Length[UpValues[GammaCC]],
  "native_cc_squared_inputform"->B1String[DiracProduct[GammaCC,GammaCC]],
  "native_cc_transpose_inputform"->B1String[Transp[GammaCC]],
  "native_y1_bilinear_inputform"->B1String[NCM[Bar[lfield],CConj[f]]],
  "native_y2_bilinear_inputform"->B1String[NCM[Bar[lfield],f]]|>;
 <|
  "format"->"t3-b1-native-cc-abstract-v2",
  "status"->"NativeAbstractCCAlgebraVerified_NumericOrientationPending",
  "model"->"A_p0","checks"->checks,
  "installed_symbol_manifest"->manifest,
  "chosen_numeric_gamma_matrix_test_performed_here"->False,
  "installed_numeric_gamma_cc_sign_measured"->False,
  "cross_package_spinor_phase_certified"->False,
  "physical_source_certified"->False,
  "production_files_modified"->False|>
];
$b1result=CheckAbort[Catch[B1Run[],"B1CCStop"],$Aborted];
If[!AssociationQ[$b1result],$b1result=<|"status"->"Blocked","reason"->"Wolfram aborted"|>];
If[!DirectoryQ[DirectoryName[$b1out]],CreateDirectory[DirectoryName[$b1out],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[$b1out,$b1result,"RawJSON"],$Failed]===$Failed,
 Print["Export failed"];Exit[7]];
Print[Lookup[$b1result,"status","Blocked"]];
If[KeyExistsQ[$b1result,"reason"],Print[$b1result["reason"]]];
Exit[If[Lookup[$b1result,"status",""]==="NativeAbstractCCAlgebraVerified_NumericOrientationPending",0,6]];

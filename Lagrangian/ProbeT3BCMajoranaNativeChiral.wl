(* Read-only B/C(-1) native Matchete chiral-algebra probe. No matching. *)
ClearAll["Global`*"];
args=Rest[$ScriptCommandLine]; If[Length[args]!=1,Print["Expected output JSON path"];Exit[2]];
out=ExpandFileName[First[args]]; dir=DirectoryName@ExpandFileName[$InputFileName];
(* Load Matchete BEFORE parsing any definitions mentioning its exported symbols.
   Otherwise Global`Bar, Global`NCM, etc. are created during parsing and
   shadow Matchete`Bar, Matchete`NCM when Needs loads the package. *)
step="load";
loaded=CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]; True],$Failed],$Aborted];
If[loaded=!=True,
 res=<|"status"->"NativeMajoranaChiralProbeBlocked","step"->"load",
   "reason"->"Matchete failed to load before native symbols were parsed"|>;
 If[!DirectoryQ[DirectoryName[out]],CreateDirectory[DirectoryName[out],CreateIntermediateDirectories->True]];
 Export[out,res,"RawJSON"];Print[res["status"]];Exit[6]
];
fail[msg_]:=Throw[<|"status"->"NativeMajoranaChiralProbeBlocked","step"->step,"reason"->msg|>,"stop"];
run[]:=Module[{results=<||>,model,built,lint,f,spin,explicit,checks},
 Get[FileNameJoin[{dir,"model","T3ModelCatalog.wl"}]];
 Get[FileNameJoin[{dir,"model","LagrangianBuilder.wl"}]];
 Do[
 step="build-"<>cls;
 model=T3ModelFromClass[cls,-1];If[!AssociationQ[model],fail["Model missing: "<>cls]];
 built=CheckAbort[Check[UsingFrontEnd[BuildT3Lagrangian[model]],$Failed],$Aborted];
 If[!AssociationQ[built]||Lookup[built,"Status",""]=!="Success"||!TrueQ[Lookup[built,"FullValidation",False]],fail["Full native build failed: "<>cls]];
 lint=built["LInt"];
 f=If[cls==="B",NewFermion[r],NewFermion[j,r]];
 spin=NCM[Bar[l[i,p]],CConj[f]];
 explicit=NCM[Bar[l[i,p]],DiracProduct[GammaCC],Transp[Bar[f]]];
 checks=<|
 "y1_and_y2_in_native_action"->(!FreeQ[lint,y1]&&!FreeQ[lint,y2]),
 "cconj_definition"->TrueQ[FullSimplify[spin==explicit]],
 "cconj_flips_PL_to_PR"->TrueQ[FullSimplify[CConj[NCM[PL,f]]==NCM[PR,CConj[f]]]],
 "cconj_flips_PR_to_PL"->TrueQ[FullSimplify[CConj[NCM[PR,f]]==NCM[PL,CConj[f]]]],
 "cc_transpose_antisymmetric"->TrueQ[FullSimplify[Transp[GammaCC]==-GammaCC]],
 "cc_square_minus_one"->TrueQ[FullSimplify[DiracProduct[GammaCC,GammaCC]==-1]]|>;
 results[cls<>"_m1"]=<|"checks"->checks,"passed"->AllTrue[Values[checks],TrueQ],"native_y1_bilinear"->ToString[spin,InputForm],"native_y1_explicit_cc"->ToString[explicit,InputForm],"native_y2_bilinear"->ToString[NCM[Bar[l[i,p]],f],InputForm]|>;
 ,{cls,{"B","C"}}];
 <|"format"->"t3-bc-native-majorana-chiral-v1","status"->If[AllTrue[Values[results],TrueQ[Lookup[#,"passed",False]]&],"NativeAbstractChiralChecksPassedPhysicalDictionaryPending","NativeChiralIdentityMismatch"],"models"->results,"concrete_CC_orientation_measured"->False,"rgbeta_cross_package_dictionary_certified"->False,"production_modified"->False|>
];
res=CheckAbort[Catch[run[],"stop"],$Aborted];If[!AssociationQ[res],res=<|"status"->"NativeMajoranaChiralProbeBlocked","step"->step,"reason"->"Aborted"|>];
If[!DirectoryQ[DirectoryName[out]],CreateDirectory[DirectoryName[out],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[out,res,"RawJSON"],$Failed]===$Failed,Exit[7]];
Print[Lookup[res,"status","Unknown"]];If[KeyExistsQ[res,"reason"],Print[res["reason"]]];
Exit[If[Lookup[res,"status",""]==="NativeAbstractChiralChecksPassedPhysicalDictionaryPending",0,6]];

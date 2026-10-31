(* Workstream B: runtime RGBeta mixed scalar mass vertex metadata and native Lam.
   Uses installed RGBeta PackageScope Lam; no invented component basis.
   wolframscript -file scripts/ProbeT3BCMixedBilinearRGBeta.wl NEW_OUTPUT.json
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine]=!=2,
 Print["Usage: ProbeT3BCMixedBilinearRGBeta.wl NEW_OUTPUT.json"];Exit[64]];
out=ExpandFileName[$ScriptCommandLine[[2]]];
If[FileExistsQ[out],Print["Refusing to overwrite ",out];Exit[64]];
root=DirectoryName@DirectoryName@ExpandFileName[$InputFileName];
loaded=CheckAbort[Check[Needs["RGBeta`"];True,$Failed],$Aborted];
If[loaded=!=True,Print["RGBeta unavailable"];Exit[3]];
Get[FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
lamName="RGBeta`PackageScope`Lam";
If[Names[lamName]=!={lamName} ||
    Length[DownValues[Symbol[lamName]]]<1,
 Print["Installed RGBeta Lam unavailable or unsupported"];Exit[4]];
lam=Symbol[lamName];
T3BRFail[model_,why_] := Throw[<|"model"->model,"status"->"Blocked",
  "error"->why,"physical_mapping_certified"->False|>,"T3BRMixed"];
T3BROne[label_,dF_] := Module[
 {built,info,fields,unique,coeff,conjQ,live,arrays,rules,relevant,ell},
 Print["Live RGBeta mixed quadratic for ",label];
 built=CheckAbort[Check[T3RGBetaBuildEFT1[2,2,dF,-1],$Failed],$Aborted];
 If[!AssociationQ[built],T3BRFail[label,"EFT1 model build failed"]];
 If[CheckAbort[Check[
   AddScalarMass[muT3GeneralS1S2,{Bar@S1,Bar@S2},
    GroupInvariant->(eps[SU2L@fund,#1,#2]&),SelfConjugate->False];True,
   $Failed],$Aborted]=!=True,
   T3BRFail[label,"AddScalarMass declaration failed"]];
 If[!KeyExistsQ[$scalarMasses,muT3GeneralS1S2],
   T3BRFail[label,"Native mass coupling missing from live registry"]];
 info=$scalarMasses[muT3GeneralS1S2];
 fields=info[Fields];
 unique=info[UniqueArrangements];
 conjQ=info[SelfConjugate];
 coeff=ToString[info[Coupling][],InputForm];
 If[fields=!={Bar[S1],Bar[S2]} || conjQ=!=False ||
   !ListQ[unique] || Length[unique]=!=12,
   T3BRFail[label,"Unexpected native fields, h.c. flag, or permutations"]];
 (* The ACTUAL installed RGBeta mass-augmented Lam, not a reimplemented
    scalar-potential action. This may be expensive but catches a missing
    registered mass vertex that metadata alone would not detect. *)
 live=CheckAbort[Check[lam[t3BM_a,t3BM_b,t3BM_c,t3BM_d,True],$Failed],$Aborted];
 If[live===$Failed || live===$Aborted ||
   FreeQ[live,muT3GeneralS1S2,{0,Infinity}],
   T3BRFail[label,"Live massive Lam does not contain added mixed mass coupling"]];
 arrays=Cases[live,_SparseArray,Infinity];
 rules=Flatten[ArrayRules/@arrays,1];
 relevant=Select[rules,!FreeQ[#,muT3GeneralS1S2,{0,Infinity}]&];
 If[Length[relevant]<1,
   T3BRFail[label,"Live sparse coupling tensor has no mixed-mass entries"]];
 <|"model"->label,"status"->"LiveRGBetaMixedBilinearRegisteredAndLamDetected",
   "rgbeta_mass_symbol"->"muT3GeneralS1S2",
   "rgbeta_fields"->(ToString[#,InputForm]&/@fields),
   "rgbeta_invariant_inputform"->ToString[info[Invariant][i,j],InputForm],
   "rgbeta_mass_coupling_function_inputform"->coeff,
   "rgbeta_unique_permutations"->Length[unique],
   "rgbeta_self_conjugate"->conjQ,
   "native_lam_symbol"->lamName,
   "native_massive_lam_evaluated"->True,
   "native_lam_mass_sparse_array_rules_count"->Length[relevant],
   "native_lam_mass_rule_samples_inputform"->(ToString[#,InputForm]&/@Take[relevant,UpTo[8]]),
   "rgbeta_potential_convention"->"+muR * eps * Bar[S1] Bar[S2] + h.c.; - sign in Lagrangian",
   "native_rgbeta_epsilon_component_phase_verified"->False,
   "physical_mapping_certified"->False|>
];
report=CheckAbort[Check[<|
 "format"->"t3-bc-live-mixed-rgbeta-bilinear-v1",
 "status"->"LiveRGBetaNativeMassTensorDiagnostic",
 "source_sha256"-><|
  "RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl"->
   IntegerString[FileHash[FileNameJoin[{root,"RGE","running","rgbeta","RunGeneralT3QuadraticRGBeta.wl"}],"SHA256"],16,64],
  "RGE/running/rgbeta/T3RGBetaModel.wl"->
   IntegerString[FileHash[FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaModel.wl"}],"SHA256"],16,64]|>,
 "wolfram_version"->$Version,
 "installed_rgbeta_package_file"->ToString[FindFile["RGBeta`"],InputForm],
 "models"-><|
   "B_-1"->Catch[T3BROne["B_-1",1],"T3BRMixed"],
   "C_-1"->Catch[T3BROne["C_-1",3],"T3BRMixed"]|>,
 "physical_mapping_certified"->False|>,$Failed],$Aborted];
If[!AssociationQ[report],Print["RGBeta mixed bilinear probe blocked"];Exit[5]];
If[!DirectoryQ[DirectoryName[out]],
 CreateDirectory[DirectoryName[out],CreateIntermediateDirectories->True]];
If[Check[Export[out,report,"RawJSON"],$Failed]===$Failed || !FileExistsQ[out],
 Print["Cannot write ",out];Exit[6]];
Print["RGBeta mixed quadratic diagnostic: ",out];
If[!AllTrue[Values[report["models"]],
 Lookup[#,"status",""]==="LiveRGBetaMixedBilinearRegisteredAndLamDetected"&],Exit[7]];

(* Workstream B: runtime RGBeta mixed scalar mass vertex metadata and native Lam.
   Uses installed RGBeta PackageScope Lam; no invented component basis.
   v4: inspect native SparseArray ArrayRules directly instead of searching
   an opaque TStructure/SparseArray with FreeQ; preserve failure diagnostics.
   DownValues uses literal package symbol; indices and native keys are safe.
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
(* DownValues has a holding attribute: DownValues[Symbol[name]] is invalid.
   The known installed v1.2.0 function is a literal, qualified symbol. *)
lamName="RGBeta`PackageScope`Lam";
If[Names[lamName]=!={lamName},
 Print["Expected RGBeta Lam function is not registered: ",lamName];Exit[4]];
lamDefinitionCount=Length[DownValues[RGBeta`PackageScope`Lam]];
If[lamDefinitionCount<1,
 Print["Installed RGBeta Lam has no definitions: ",lamName];Exit[4]];
Print["Resolved native RGBeta Lam: ",lamName,
 " (downvalues=",lamDefinitionCount,")"];

T3BRFail[model_,why_,details_:<||>] := Throw[
 Join[<|"model"->model,"status"->"Blocked", "error"->why,
   "physical_mapping_certified"->False|>,details],"T3BRMixed"];

(* RGBeta`PackageScope`Fields, Invariant, Coupling, ... are NOT the same
   keys as a Global`Fields symbol. This follows scripts/RGBetaNativeAccess.wl.
   Select the unique native key by its unqualified, explicitly checked name. *)
T3BRNativeKeyName[key_] := If[StringQ[key],key,
 Last[StringSplit[ToString[key,InputForm],"`"]]];
T3BRNativeGet[model_String,info_Association,key_String] := Module[
 {matches,value},
 matches=Select[Keys[info],T3BRNativeKeyName[#]===key&];
 If[Length[matches]=!=1,
   T3BRFail[model,"Missing or ambiguous native RGBeta key "<>key<>
     "; available="<>StringRiffle[T3BRNativeKeyName/@Keys[info],","]]];
 value=Lookup[info,First[matches],Missing["MissingNativeKey",key]];
 If[MissingQ[value] || value===$Failed || value===$Aborted,
   T3BRFail[model,"Native RGBeta value unavailable: "<>key]];
 value
];
T3BROne[label_,dF_] := Module[
 {built,added,info,fields,unique,coeff,conjQ,live,arrays,rules,
  relevant,invFn,coupFn,diag,ruleText,lamText},
 Print["Live RGBeta mixed quadratic for ",label];
 built=CheckAbort[Check[T3RGBetaBuildEFT1[2,2,dF,-1],$Failed],$Aborted];
 If[!AssociationQ[built],T3BRFail[label,"EFT1 model build failed"]];
 added=CheckAbort[Check[
   AddScalarMass[muT3GeneralS1S2,{Bar@S1,Bar@S2},
    GroupInvariant->(eps[SU2L@fund,#1,#2]&),SelfConjugate->False],
   $Failed],$Aborted];
 If[added===$Failed || added===$Aborted,
   T3BRFail[label,"AddScalarMass declaration failed"]];
 If[!AssociationQ[$scalarMasses] ||
    !KeyExistsQ[$scalarMasses,muT3GeneralS1S2],
   T3BRFail[label,"Native mass coupling missing from live registry"]];
 info=$scalarMasses[muT3GeneralS1S2];
 If[!AssociationQ[info],T3BRFail[label,"Invalid native mass record"]];
 fields=T3BRNativeGet[label,info,"Fields"];
 unique=T3BRNativeGet[label,info,"UniqueArrangements"];
 conjQ=T3BRNativeGet[label,info,"SelfConjugate"];
 invFn=T3BRNativeGet[label,info,"Invariant"];
 coupFn=T3BRNativeGet[label,info,"Coupling"];
 coeff=ToString[coupFn[],InputForm];
 If[fields=!={Bar[S1],Bar[S2]} || conjQ=!=False ||
   !ListQ[unique] || Length[unique]=!=12,
   T3BRFail[label,"Unexpected native fields, h.c. flag, or permutations"]];
 (* Live mass contribution must appear in the package's own Lam[...,True]
    output. RGBeta returns a TStructure backed by SparseArray. FreeQ over
    the unexpanded top-level object is NOT a reliable scalar-mass witness.
    Inspect ArrayRules on the ACTUAL native arrays, including default rules,
    then require an explicit occurrence of the registered mass symbol in
    a non-default sparse rule. Never reconstruct a hypothetical Lam tensor. *)
 live=CheckAbort[Check[
   RGBeta`PackageScope`Lam[t3BMa,t3BMb,t3BMc,t3BMd,True],
   $Failed],$Aborted];
 If[live===$Failed || live===$Aborted,
   T3BRFail[label,"Native massive Lam evaluation failed or aborted"]];
 arrays=Cases[live,_SparseArray,{0,Infinity}];
 rules=Flatten[ArrayRules/@arrays,1];
 (* String-based search is used only to find the *named symbolic coupling*
    inside the InputForm of actual ArrayRules. The Python verifier separately
    checks the native field/coupling factor and invariant, so this does not
    declare a physical tensor or an independently normalized vertex. *)
 relevant=Select[rules,
   !FreeQ[#,muT3GeneralS1S2,{0,Infinity}] ||
    StringContainsQ[ToString[#,InputForm],"muT3GeneralS1S2"]&];
 ruleText=ToString[#,InputForm]& /@ Take[rules,UpTo[12]];
 lamText=ToString[live,InputForm];
 diag=<|
   "native_lam_head_inputform"->ToString[Head[live],InputForm],
   "native_lam_leaf_count"->LeafCount[live],
   "native_lam_top_level_freeq_mass_symbol"->
     FreeQ[live,muT3GeneralS1S2,{0,Infinity}],
   "native_lam_sparse_array_count"->Length[arrays],
   "native_lam_total_sparse_array_rules_count"->Length[rules],
   "native_lam_mass_sparse_array_rules_count"->Length[relevant],
   "native_lam_sparse_rule_examples_inputform"->ruleText,
   "native_lam_excerpt_inputform"->
     StringTake[lamText,UpTo[1500]],
   "native_mass_registry_keys_inputform"->
     (ToString[#,InputForm]&/@Keys[$scalarMasses])|>;
 Print[label,": massive Lam sparse arrays=",Length[arrays],
   "; native rules=",Length[rules],
   "; mixed-mass rule witnesses=",Length[relevant],
   "; top-level FreeQ=",diag["native_lam_top_level_freeq_mass_symbol"]];
 If[Length[relevant]<1,
   T3BRFail[label,
    "Native massive Lam has no mixed-mass coupling in inspected sparse rules",
    diag]];
 <|"model"->label,"status"->"LiveRGBetaMixedBilinearRegisteredAndLamDetected",
   "rgbeta_mass_symbol"->"muT3GeneralS1S2",
   "rgbeta_fields"->(ToString[#,InputForm]&/@fields),
   "rgbeta_invariant_inputform"->ToString[invFn[i,j],InputForm],
   "rgbeta_mass_coupling_function_inputform"->coeff,
   "rgbeta_unique_permutations"->Length[unique],
   "rgbeta_self_conjugate"->conjQ,
   "native_lam_symbol"->lamName,
   "native_lam_definition_count"->lamDefinitionCount,
   "native_massive_lam_evaluated"->True,
   "native_lam_mass_sparse_array_rules_count"->Length[relevant],
   "native_lam_runtime_diagnostics"->diag,
   "native_lam_mass_rule_samples_inputform"->(ToString[#,InputForm]&/@Take[relevant,UpTo[8]]),
   "rgbeta_potential_convention"->"+muR * eps * Bar[S1] Bar[S2] + h.c.; - sign in Lagrangian",
   "native_rgbeta_epsilon_component_phase_verified"->False,
   "physical_mapping_certified"->False|>
];
(* A Wolfram message inside a native Lam call should block only that model.
   Do not let the outer Check erase BOTH model reports and produce no JSON. *)
T3BRSafeOne[label_String,dF_Integer] := Module[{result},
 result=CheckAbort[
   Check[Catch[T3BROne[label,dF],"T3BRMixed"],$Failed],
   $Aborted];
 If[AssociationQ[result],Return[result]];
 <|"model"->label,"status"->"Blocked",
   "error"->"Live RGBeta evaluation emitted a message, failed, or aborted; inspect Wolfram terminal output",
   "physical_mapping_certified"->False|>
];
report=CheckAbort[<|
 "format"->"t3-bc-live-mixed-rgbeta-bilinear-v1",
 "status"->"LiveRGBetaNativeMassTensorDiagnostic",
 "source_sha256"-><|
  "RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl"->
   IntegerString[FileHash[FileNameJoin[{root,"RGE","running","rgbeta","RunGeneralT3QuadraticRGBeta.wl"}],"SHA256"],16,64],
  "RGE/running/rgbeta/T3RGBetaModel.wl"->
   IntegerString[FileHash[FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaModel.wl"}],"SHA256"],16,64]|>,
 "wolfram_version"->$Version,
 "probe_revision"->"inspect-native-sparse-rules-with-block-diagnostics-v4",
 "installed_rgbeta_package_file"->ToString[FindFile["RGBeta`"],InputForm],
 "models"-><|
   "B_-1"->T3BRSafeOne["B_-1",1],
   "C_-1"->T3BRSafeOne["C_-1",3]|>,
 "physical_mapping_certified"->False|>,$Aborted];
If[!AssociationQ[report],Print["RGBeta mixed bilinear probe blocked"];Exit[5]];
If[!DirectoryQ[DirectoryName[out]],
 CreateDirectory[DirectoryName[out],CreateIntermediateDirectories->True]];
If[Check[Export[out,report,"RawJSON"],$Failed]===$Failed || !FileExistsQ[out],
 Print["Cannot write ",out];Exit[6]];
Print["RGBeta mixed quadratic diagnostic: ",out];
If[!AllTrue[Values[report["models"]],
 Lookup[#,"status",""]==="LiveRGBetaMixedBilinearRegisteredAndLamDetected"&],
 Print["RGBeta mixed bilinear probe blocked; model diagnostic JSON saved"];Exit[7]];
Print["Both native B/C massive Lam witnesses detected; physical map UNCERTIFIED"];


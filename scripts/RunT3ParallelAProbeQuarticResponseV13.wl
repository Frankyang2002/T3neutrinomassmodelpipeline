(* Workstream A v13. Read-only forensic re-probe of a FAILED v10 smoke.
   No RGBeta redefinitions, no component axes, no Wilson evolution.
   Usage: wolframscript -file scripts/RunT3ParallelAProbeQuarticResponseV13.wl FORENSICS.json OUTPUT.json
   The v13 forensic JSON must be source-checked by the Python analyzer first.
   This probe captures LIMITED inputform diagnostics, NOT a pass certificate. *)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 3,
 Print["Usage: native v13 forensic input JSON, new re-probe output JSON"];Exit[64]];
T3A13Input=ExpandFileName[$ScriptCommandLine[[2]]];
T3A13Out=ExpandFileName[$ScriptCommandLine[[3]]];
If[!FileExistsQ[T3A13Input] || FileExistsQ[T3A13Out] ||
   FileExistsQ[T3A13Out<>".partial"],
 Print["v13 BLOCKED: input missing or output already exists"];Exit[6]];
T3A13Root=DirectoryName@DirectoryName@ExpandFileName[$InputFileName];
T3A13Sources={"RGE/running/rgbeta/T3RGBetaModel.wl",
 "scripts/RGBetaNativeAccess.wl",
 "scripts/RunT3ParallelANativeQuarticProjectorsV10.wl"};
T3A13Smoke=FileNameJoin[{T3A13Root,"output","parallel_t3","wilson",
 "rgbeta_quartic_projectors_native_v10_smoke.json"}];
T3A13Flags={"T3_EXPERIMENTAL_GENERAL_FERMIONS",
 "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
 "T3_EXPERIMENTAL_GENERAL_QUARTICS",
 "T3_EXPERIMENTAL_FLAVOR_MASSES"};
If[AnyTrue[T3A13Flags, Environment[#] === "1" &],
 Print["v13 BLOCKED: optional T3 extensions enabled"];Exit[3]];
T3A13Read=Quiet@Check[Import[T3A13Input,"RawJSON"],$Failed];
If[!AssociationQ[T3A13Read] ||
 Lookup[T3A13Read,"format",None]=!="t3-parallel-a-native-rgbeta-quartic-forensics-v13" ||
 Lookup[T3A13Read,"status",None]=!="NativeQuarticSmokeFailureClassified_NotCertified" ||
 Lookup[T3A13Read,"physical_Wilson_evolution_authorized",True]=!=False,
 Print["v13 BLOCKED: invalid forensic input status"];Exit[2]];
If[!FileExistsQ[T3A13Smoke] ||
 IntegerString[FileHash[T3A13Smoke,"SHA256"],16,64]=!=
 Lookup[T3A13Read,"input_smoke_sha256",None],
 Print["v13 BLOCKED: original failed native smoke missing or changed"];Exit[2]];
T3A13Hashes=Lookup[T3A13Read,"input_smoke_source_sha256",<||>];
If[!AssociationQ[T3A13Hashes] ||
 Sort[Keys[T3A13Hashes]]=!=Sort[T3A13Sources] ||
 !AllTrue[T3A13Sources,Function[s,
  FileExistsQ[FileNameJoin[{T3A13Root,s}]] &&
  Lookup[T3A13Hashes,s,""]===IntegerString[
    FileHash[FileNameJoin[{T3A13Root,s}],"SHA256"],16,64]]],
 Print["v13 BLOCKED: native v10 source provenance changed"];Exit[2]];
If[CheckAbort[Check[Needs["RGBeta`"];True,$Failed],$Aborted]=!=True,
 Print["v13 BLOCKED: installed RGBeta unavailable"];Exit[4]];
T3A13Package=Quiet@Check[FindFile["RGBeta`"],$Failed];
If[!StringQ[T3A13Package] || !FileExistsQ[T3A13Package] ||
 IntegerString[FileHash[T3A13Package,"SHA256"],16,64]=!=
 Lookup[T3A13Read,"installed_RGBeta_main_file_sha256",None],
 Print["v13 BLOCKED: native RGBeta package differs from failed smoke"];Exit[2]];
Get[FileNameJoin[{T3A13Root,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
Get[FileNameJoin[{T3A13Root,"scripts","RGBetaNativeAccess.wl"}]];

ClearAll[T3A13Text,T3A13Compact,T3A13Diag,T3A13Replay,T3A13Model];
T3A13Text[x_] := ToString[x,InputForm];
T3A13Compact[x_] := StringTake[T3A13Text[x],UpTo[900]];
T3A13Diag[expr_] := Module[{refined,without,one,two,numerical},
 If[expr===$Failed || expr===$Aborted,
   Return[<|"native_execution_returned_failure"->True|>]];
 refined=TimeConstrained[Quiet@Check[
   FullSimplify[RefineGroupStructures[
     expr /. (Matrix|Tensor)[x_][__] :> x]],$Failed],45,$Failed];
 If[refined===$Failed || refined===$Aborted,
   Return[<|"raw_head"->T3A13Text[Head[expr]],
     "raw_inputform_prefix"->T3A13Compact[expr],
     "refinement_failed_or_timeout"->True|>]];
 without=TimeConstrained[Quiet@Check[FullSimplify[refined/.T3A13Src->0],$Failed],15,$Failed];
 one=TimeConstrained[Quiet@Check[FullSimplify[refined/.T3A13Src->1],$Failed],15,$Failed];
 two=TimeConstrained[Quiet@Check[FullSimplify[refined/.T3A13Src->2],$Failed],15,$Failed];
 numerical=Quiet@Check[N[one,15],$Failed];
 <|"raw_head"->T3A13Text[Head[expr]],
   "raw_inputform_prefix"->T3A13Compact[expr],
   "refined_inputform_prefix"->T3A13Compact[refined],
   "constant_prefix"->T3A13Compact[without],
   "coefficient_prefix"->T3A13Compact[one],
   "double_prefix"->T3A13Compact[two],
   "linear_in_source"->TrueQ[PossibleZeroQ[without]] &&
     TrueQ[PossibleZeroQ[two-2 one]],
   "native_coefficient_numeric"->NumericQ[numerical] &&
     TrueQ[Abs[Im[numerical]]<10^-8],
   "native_coefficient_value"->If[NumericQ[numerical] &&
     TrueQ[Abs[Im[numerical]]<10^-8],N[Re[numerical]],Null]|>
];
(* T3A13Src is a unique placeholder to avoid local symbol scoping in test.
   It is substituted only AFTER native expressions have been built. *)
T3A13Replay[target_Symbol,source_Symbol] := Module[
 {proj,isolated,full,rules,before},
 before=Length[Keys[$quartics]];
 proj=Quiet@Check[T3Stage22NativeGet[T3Stage22NativeInfo[target],"Projector"],$Failed];
 If[proj===$Failed,Return[<|"projector_lookup_failed"->True|>]];
 isolated=TimeConstrained[Quiet@Check[
   Block[{$quartics=KeyTake[$quartics,{source}]},
     Activate[proj[t13a,t13b,t13c,t13d]] *
       Lam[t13a,t13b,t13c,t13d]],$Failed],55,$Failed];
 (* Independent package route: keep full RGBeta quartic registry and set
    all non-source coupling symbols to zero *after* building native Lam.
    This may be slow. The two routes are never substituted for proof. *)
 rules=Thread[DeleteCases[Keys[$quartics],source] -> 0];
 full=TimeConstrained[Quiet@Check[
   Activate[proj[t13a,t13b,t13c,t13d]] *
     Lam[t13a,t13b,t13c,t13d],$Failed],55,$Failed];
 <|"target"->SymbolName[target],"source"->SymbolName[source],
  "isolated_registry"->T3A13Diag[isolated/.source->T3A13Src],
  "full_registry_then_other_couplings_zero"->
    T3A13Diag[If[full===$Failed||full===$Aborted,full,
      (full/.rules)/.source->T3A13Src]],
  "original_registry_restored_after_probe"->(Length[Keys[$quartics]]===before)|>
];
T3A13Model[tag_String,probes_List] := Module[{cls,alpha,dims,build,rows},
 cls=StringTake[tag,1]; alpha=ToExpression[Last[StringSplit[tag,"_"]]];
 dims=Switch[cls,"A",{1,3,2},"B",{2,2,1},"C",{2,2,3},
   "D",{3,1,2},"E",{3,3,2},_,{}];
 If[Length[dims]=!=3 || Length[probes]<1 || Length[probes]>3,
   Return[<|"ok"->False,"reason"->"invalid_model_or_probe_count"|>]];
 build=Quiet@Check[T3RGBetaBuildEFT1@@Append[dims,alpha],$Failed];
 If[!AssociationQ[build],
   Return[<|"ok"->False,"reason"->"native_EFT1_build_failed"|>]];
 rows=Table[Module[{target,source,rr},
   target=Select[Keys[$quartics],SymbolName[#]===Lookup[sel,"target",""]&];
   source=Select[Keys[$quartics],SymbolName[#]===Lookup[sel,"source",""]&];
   If[Length[target]=!=1 || Length[source]=!=1,
     <|"ok"->False,"reason"->"unregistered_target_or_source"|>,
     rr=TimeConstrained[Quiet@Check[
       T3A13Replay[First[target],First[source]],$Failed],180,$Failed];
     If[AssociationQ[rr],rr,
        <|"ok"->False,"reason"->"reprobe_exception_or_timeout"|>]
   ]],{sel,probes}];
 <|"ok"->True,"dimensions"->dims,"alpha"->alpha,"probes"->rows|>
];
T3A13Selection=Lookup[T3A13Read,"minimal_reprobe_targets",{}];
If[!ListQ[T3A13Selection] || Length[T3A13Selection]<1 ||
 Length[T3A13Selection]>3 ||
 !AllTrue[T3A13Selection,AssociationQ[#] &&
   MemberQ[{"A_-4","B_-1","E_-2"},Lookup[#,"model",None]]&],
 Print["v13 BLOCKED: unsafe/empty re-probe selection"];Exit[2]];
T3A13Groups=GroupBy[T3A13Selection,Lookup[#,"model",""]&];
T3A13Models=Association@Table[
 model->T3A13Model[model,Lookup[T3A13Groups,model,{}]],
 {model,Keys[T3A13Groups]}];
T3A13Report=<|
 "format"->"t3-parallel-a-native-rgbeta-quartic-response-reprobe-v13",
 "status"->"NativeResponseForensicsCaptured_NoPhysicsAuthorization",
 "forensic_input_sha256"->IntegerString[FileHash[T3A13Input,"SHA256"],16,64],
 "old_smoke_source_sha256"->T3A13Hashes,
 "installed_RGBeta_main_file_sha256"->Lookup[T3A13Read,"installed_RGBeta_main_file_sha256",None],
 "probe_source_sha256"->IntegerString[FileHash[ExpandFileName[$InputFileName],"SHA256"],16,64],
 "models"->T3A13Models,
 "native_projector_150_diagonals_certified"->False,
 "native_complex_quartic_projectors_certified"->False,
 "physical_Wilson_evolution_authorized"->False,
 "EFT_order"->"UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT"|>;
If[!DirectoryQ[DirectoryName[T3A13Out]],
 CreateDirectory[DirectoryName[T3A13Out],CreateIntermediateDirectories->True]];
T3A13Tmp=T3A13Out<>".partial";
If[Quiet@Check[Export[T3A13Tmp,T3A13Report,"RawJSON"],$Failed]===$Failed ||
 !FileExistsQ[T3A13Tmp],
 Print["v13 BLOCKED: response JSON output failed"];Exit[9]];
If[Quiet@Check[RenameFile[T3A13Tmp,T3A13Out],$Failed]===$Failed,
 Print["v13 BLOCKED: response JSON rename failed"];Exit[9]];
Print["v13 native differential reprobe: ",Length[T3A13Selection],
 " selected source/target couplings -> ",T3A13Out];
Print["Forensic evidence only: native v10 and physical Wilson evolution STILL BLOCKED."];
Exit[0];

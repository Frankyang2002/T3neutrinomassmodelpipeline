(*
 Read-only native RGBeta quartic component *probe*, all 16 T3 models.
 This is intentionally NOT the existing declarations-only exporter.

 Its contract is strict: evaluate native RGBeta `Invariant` functions at
 explicit representation indices, retain coupling and arrangement metadata.
 If a nontrivial tensor is not natively numeric, emit an actionable BLOCKED
 status; NEVER substitute guessed delS2/eps phase or invented prefactors.

 wolframscript -file scripts/RunT3ParallelAExportRGBetaComponents.wl \
   output/parallel_t3/wilson/native_rgbeta_components_probe.json
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
 Print["Usage: RunT3ParallelAExportRGBetaComponents.wl OUTPUT.json"];Exit[64]
];
outFile=ExpandFileName[$ScriptCommandLine[[2]]];
root=DirectoryName@DirectoryName@ExpandFileName[$InputFileName];
registry=Quiet@Check[Import[FileNameJoin[{root,"RGE","general","generated",
                                    "general_rgbeta_quartic_registry.json"}],"RawJSON"],$Failed];
If[!AssociationQ[registry] || Length[Keys[registry["models"]]] =!= 16,
 Print["Blocked_RegistrySourceMissingOrIncomplete"];Exit[2]];
If[Check[Needs["RGBeta`"];True,$Failed] =!= True,
 Print["Blocked_RGBetaNativePackageUnavailable"];Exit[3]];
Get[FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
Get[FileNameJoin[{root,"scripts","RGBetaNativeAccess.wl"}]];
(* Pure native evaluation: NO hand-coded CG tensors or normalization guesses. *)
ClearAll[T3ANativeRecord];
T3ANativeRecord[symbol_Symbol,modelDims_List] := Module[
 {info,fields,inv,coupling,perms,tuples,values,numeric,components,dims,base,names},
 info=T3Stage22NativeInfo[symbol];
 fields=T3Stage22NativeGet[info,"Fields"];
 inv=T3Stage22NativeGet[info,"Invariant"];
 coupling=T3Stage22NativeGet[info,"Coupling"];
 perms=T3Stage22NativeGet[info,"UniqueArrangements"];
 base=fields /. Bar[x_] :> x;
 names=If[Length[base]===4, (If[Head[#]===Symbol, SymbolName[#], "UNKNOWN"] & /@ base),{}];
 dims=names /. {"H"->2,"S1"->modelDims[[1]],"S2"->modelDims[[2]]};
 If[Length[fields]=!=4 || Length[dims]=!=4 || !And@@(IntegerQ[#] && #>0 & /@ dims),
  Return[<|"symbol"->SymbolName[symbol],"status"->"Blocked_QuarticFieldDimensionMismatch"|>]];
 tuples=Tuples[Range/@dims];
 values=Quiet[Check[FullSimplify[inv@@#]&/@tuples,$Failed]];
 numeric=ListQ[values] && Length[values]===Times@@dims && AllTrue[values,NumericQ];
 components=If[numeric,ToString[#,InputForm]&/@values,{}];
 <|"symbol"->SymbolName[symbol],
   "status"->If[numeric,"NativeInvariantComponentsEvaluated_PrefactorNotCertified",
                      "Blocked_NativeInvariantNotNumeric"],
   "native_fields"->ToString[fields,InputForm],
   "dimensions"->dims,
   "components"->components,
   "unresolved_invariant_sample"->If[numeric,"",ToString[Take[values,UpTo[3]],InputForm]],
   "native_coupling_function"->ToString[coupling,InputForm],
   "native_unique_arrangements"->ToString[perms,InputForm],
   "potential_prefactor_status"->"Blocked_RequiresNativeQuarticActionNormalization",
   "same_component_frame_certified"->False|>
];
ClearAll[T3AProbe];
T3AProbe[]:=Module[{all=<||>,settings,meta,build,ext,new,old,expected,
                         symbols,records,dimList,modelName,counts,physicalDim},
 Do[
  Print["T3 A RGBeta native component probe: ",modelName];
  settings=registry["models",modelName];
  meta=Join[settings["dimensions"],{settings["alpha"]}];
  build=Quiet@Check[Apply[T3RGBetaBuild,meta],$Failed];
  If[!AssociationQ[build],Return[$Failed]];
  ext=FileNameJoin[{root,"RGE","general","generated",
                    modelName<>"_extra_quartics.wl"}];
  T3GeneralQuarticAddedSymbols={};
  If[Check[Get[ext];True,$Failed]=!=True,Return[$Failed]];
  new=Lookup[settings["new_rgbeta_symbols"],"name",{}];
  old=Lookup[settings["preexisting_rgbeta_symbols"],"existing_name",{}];
  expected=Join[new,old];
  symbols=Symbol/@expected;
  If[!And@@(KeyExistsQ[$quartics,#]&/@symbols),Return[$Failed]];
  records=Table[
   counts=If[i<=Length[new],settings["new_rgbeta_symbols"][[i,"counts"]],
       settings["preexisting_rgbeta_symbols"][[i-Length[new],"counts"]]];
   (* Ordered four scalar dimensions are model-dependent. Current inventory
      stores counts, while native Fields can permute the slots. Do not deduce
      actual field axes from counts alone; if unavailable, BLOCK. *)
   T3ANativeRecord[symbols[[i]],settings["dimensions"]], {i,Length[expected]}];
  AssociateTo[all,modelName-><|
      "status"->If[AllTrue[records,Lookup[#,"status",""]===
          "NativeInvariantComponentsEvaluated_PrefactorNotCertified"&],
                     "NativeComponentsOnly_NotNormalized",
                     "Blocked_NativeCGNumericsOrBasis"],
      "records"->records|>],
  {modelName,Sort[Keys[registry["models"]]]}];
 <|"format"->"t3-parallel-a-native-rgbeta-component-probe-v1",
   "native_source"->"RGBeta $quartics Association, no symbol-name matching",
   "status"->"NotAPhysicalDictionary_NativePrefactorsAndCGFrameUncertified",
   "model_count"->Length[Keys[all]],"models"->all|>
];
If[FileExistsQ[outFile],Print["Blocked_OutputExists: ",outFile];Exit[5]];
result=Quiet@Check[Catch[T3AProbe[],"T3Stage22Native"],$Failed];
If[!AssociationQ[result],Print["Blocked_NativeTensorProbeFailed"];Exit[6]];
If[!DirectoryQ[DirectoryName[outFile]],CreateDirectory[DirectoryName[outFile],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[outFile,result,"RawJSON"],$Failed]===$Failed,
 Print["Blocked_NativeOutputWrite"];Exit[7]];
Print["Native tensor probe written (NOT normalization certificate): ",outFile];

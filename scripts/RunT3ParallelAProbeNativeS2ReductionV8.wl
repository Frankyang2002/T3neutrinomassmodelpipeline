(* Workstream A v8 -- LIVE RGBeta symbolic symmetric-square reduction.
   This does NOT assign numbers to RGBeta S2 components and does not modify
   the model, produce RGEs, or perform EFT matching.
   Run: wolframscript -file scripts/RunT3ParallelAProbeNativeS2ReductionV8.wl OUTPUT.json
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: RunT3ParallelAProbeNativeS2ReductionV8.wl OUTPUT.json"]; Exit[64]];
out = ExpandFileName[$ScriptCommandLine[[2]]];
If[FileExistsQ[out] || FileExistsQ[out <> ".partial"],
  Print["A v8: refusing to overwrite existing output ",out];Exit[6]];
root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
sources = {"RGE/running/rgbeta/T3RGBetaModel.wl",
  "scripts/RunT3ParallelAProbeNativeS2ReductionV8.wl"};
If[!AllTrue[sources,FileExistsQ[FileNameJoin[{root,#}]] &],
  Print["A v8: source file missing"];Exit[2]];
load = CheckAbort[Check[Needs["RGBeta`"];True,$Failed],$Aborted];
If[load =!= True,Print["A v8: installed RGBeta unavailable"];Exit[3]];
Get[FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
build = CheckAbort[Check[T3RGBetaBuildEFT1[1,3,2,-2],$Failed],$Aborted];
If[!AssociationQ[build],Print["A v8: native A_-2 EFT1 registration failed"];Exit[4]];
fundRep = SU2L[fund];tripRep = SU2L[S2];
representationValid = TrueQ[
  Dim[fundRep]===2 && TraceNormalization[fundRep]===1/2 && Casimir2[fundRep]===3/4];
tripletValid = TrueQ[
  Dim[tripRep]===3 && TraceNormalization[tripRep]===2 && Casimir2[tripRep]===2];
gen = tGen[tripRep,v8A,v8u,v8v];
cg = delS2[SU2L,v8u,v8i,v8j];
genAbstract = TrueQ[Head[gen]===tGen];
cgAbstract = TrueQ[Head[cg]===delS2];
isometry = TrueQ[
  delS2[SU2L,v8r,v8i,v8j] delS2[SU2L,v8s,v8i,v8j] ===
  del[tripRep,v8r,v8s]];
If[!And[representationValid,tripletValid,genAbstract,cgAbstract,isometry],
  Print["A v8: native SU2 S2 preliminary identity changed"];Exit[4]];
(* RefineGroupStructures is PackageScope, not a user-visible numerical matrix API.
   Resolve the ACTUALLY loaded symbol, never create a shadow definition. *)
refinerNames = Names["RGBeta`*`RefineGroupStructures"];
uniqueRefiner = TrueQ[Length[refinerNames]===1 &&
    StringStartsQ[First[refinerNames],"RGBeta`"]];
If[uniqueRefiner,
  withRefiner = CheckAbort[Check[
     With[{f=ToExpression[First[refinerNames]]},
       {f[gen], f[cg]}],$Failed],$Aborted],
  withRefiner = $Failed
];
If[ListQ[withRefiner] && Length[withRefiner]===2,
  {genRefined,cgRefined}=withRefiner,
  genRefined=gen;cgRefined=cg
];
(* Preserve distinct symbolic external labels; RGBeta's integer-index
   probes are NOT component evaluation because repeated labels are summed. *)
genReduced = TrueQ[
  FreeQ[genRefined,tGen[tripRep,___]] &&
  !FreeQ[genRefined,tGen[fundRep,___]] && genRefined=!=gen];
cgReduced = TrueQ[
  FreeQ[cgRefined,delS2[___]] && cgRefined=!=cg];
reduction = genReduced && cgReduced;
report = <|
  "format"->"t3-parallel-a-native-rgbeta-s2-reduction-probe-v8",
  "status"->If[reduction,
     "NativeRGBetaS2SymbolicReductionObserved_NumericalFrameUnspecified",
     "NativeRGBetaS2SymbolicReductionNotExposed_NumericalFrameUnspecified"],
  "native_model"->"A_-2_EFT1_SM_S1_S2",
  "installed_RGBeta_package_file"->ToString[FindFile["RGBeta`"],InputForm],
  "native_refiner_symbol_names"->refinerNames,
  "native_refiner_unique"->uniqueRefiner,
  "fundamental_dimension_and_trace_valid"->representationValid,
  "S2_dimension_trace_and_casimir_valid"->tripletValid,
  "native_s2_generator_still_symbolic"->genAbstract,
  "native_delS2_still_symbolic"->cgAbstract,
  "native_S2_symmetric_isometry_valid"->isometry,
  "native_s2_generator_reduction_observed"->genReduced,
  "native_delS2_reduction_observed"->cgReduced,
  "native_s2_generator_raw_inputform"->ToString[gen,InputForm],
  "native_s2_generator_refined_inputform"->ToString[genRefined,InputForm],
  "native_delS2_raw_inputform"->ToString[cg,InputForm],
  "native_delS2_refined_inputform"->ToString[cgRefined,InputForm],
  "basis_independent_K_eigenvalues_expected"->{"-1","-1","1/2","1/2","1/2","1/2"},
  "native_s2_numeric_components_exposed"->False,
  "native_s2_numeric_frame_certified"->False,
  "native_s2_to_matchete_intertwiner_certified"->False,
  "physical_coupling_dictionary_certified"->False,
  "physical_wilson_evolution_authorized"->False,
  "source_sha256"->Association@Table[
     path->IntegerString[FileHash[FileNameJoin[{root,path}],"SHA256"],16,64],
     {path,sources}]
|>;
If[!DirectoryQ[DirectoryName[out]],
  CreateDirectory[DirectoryName[out],CreateIntermediateDirectories->True]];
tmp=out <> ".partial";
wrote=CheckAbort[Check[Export[tmp,report,"RawJSON"],$Failed],$Aborted];
If[wrote===$Failed || !FileExistsQ[tmp],
  Print["A v8: JSON write failed"];Exit[7]];
check=Quiet@Check[Import[tmp,"RawJSON"],$Failed];
If[!AssociationQ[check] || check["format"] =!= report["format"],
  Print["A v8: JSON readback failed"];Exit[8]];
If[Quiet@Check[RenameFile[tmp,out],$Failed] === $Failed,
  Print["A v8: rename failed"];Exit[9]];
Print["A v8 native S2 symbolic reduction: ",report["status"]];
Print["Generator reduced: ",genReduced,"; delS2 reduced: ",cgReduced];
Print["Numerical RGBeta S2 frame: NOT EXPOSED. Output: ",out];
Exit[0];

(* Live RGBeta Lam diagnostic for the B_-1/C_-1 110101 quartic sector.
   This is source evidence, NOT a physical potential/coupling certificate.
   Run from repo root:
   wolframscript -file scripts/ProbeT3BCNativeRGBetaLam.wl output/audits/t3_bc_native_lam/native_lam_v1.json
   Existing reports are never overwritten. *)

ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: ProbeT3BCNativeRGBetaLam.wl NEW_OUTPUT.json"];
  Exit[64]
];
outFile = ExpandFileName[$ScriptCommandLine[[2]]];
If[FileExistsQ[outFile],
  Print["Output already exists; choose a new file: ",outFile]; Exit[64]];
root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];

loaded = CheckAbort[Check[Needs["RGBeta`"]; True,$Failed],$Aborted];
If[loaded =!= True, Print["Cannot load installed RGBeta"]; Exit[3]];
Get[FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
Get[FileNameJoin[{root,"scripts","RGBetaNativeAccess.wl"}]];

ClearAll[T3BCLamFail,T3BCLamProbeOne];
T3BCLamFail[msg_] := (Print["ERROR: ",msg];Throw[$Failed,"T3BCLamProbe"]);

(* RGBeta declares PackageScope["Lam"] in FieldsAndCouplings.m.
   Wolfram's structured-package convention creates that symbol in the
   RGBeta`PackageScope` context, NOT RGBeta` or RGBeta`Private`.
   The previous narrow Names[] search therefore returned {} and exited.
   Resolve the known declared symbol only after RGBeta is loaded; fail
   closed if it is absent or has no implementation. No dummy Lam is made. *)
lamName = "RGBeta`PackageScope`Lam";
If[Names[lamName] =!= {lamName},
  Print["Installed RGBeta does not expose expected package-scoped Lam: ",lamName];
  Print["Available symbols ending in `Lam: ",Names["*`Lam"]];
  Exit[4]
];
lamFunction = Symbol[lamName];
lamDefinitionCount = With[{f=lamFunction},Length[DownValues[f]]];
If[lamDefinitionCount < 1,
  Print["Installed RGBeta Lam has no downvalues: ",lamName]; Exit[4]
];
Print["Resolved installed RGBeta Lam: ",lamName,
      " (downvalues=",lamDefinitionCount,")"];

T3BCLamProbeOne[label_String,dF_Integer] := Module[
  {built, extension, names, records, expr, repr, sparse, dims, samples,
   relevant, verdict, symbols, group, active},
  Print["Evaluating native RGBeta Lam for ",label, " ..."];
  built = CheckAbort[Check[T3RGBetaBuild[2,2,dF,-1],$Failed],$Aborted];
  If[!AssociationQ[built],T3BCLamFail[label<>": RGBeta model build failed"]];
  extension=FileNameJoin[{root,"RGE","general","generated",label<>"_extra_quartics.wl"}];
  If[!FileExistsQ[extension],T3BCLamFail[label<>": missing generated quartic source"]];
  loaded=CheckAbort[Check[Get[extension];True,$Failed],$Aborted];
  If[loaded =!= True,T3BCLamFail[label<>": generated extension failed"]];
  names={"lambdaT3GeneralQ110101I1","lambdaHHdagS1barS2barCross"};
  symbols=Symbol/@ names;
  If[!And@@(KeyExistsQ[$quartics,#]&/@symbols),
    T3BCLamFail[label<>": missing registered 110101 quartics"]];
  records=T3Stage22NativeRecord/@symbols;
  If[!And@@(StringReplace[Lookup[#,"fields_inputform",""],WhitespaceCharacter..->""]===
                   "{H,Bar[H],Bar[S1],Bar[S2]}" &/@ records),
    T3BCLamFail[label<>": unexpected native field order"]];
  If[!And@@(StringContainsQ[Lookup[#,"coupling_factor_inputform",""],
                                     "6*Matrix["] &/@ records),
    T3BCLamFail[label<>": unexpected registered coupling normalization"]];

  (* ACTUAL runtime call into the installed RGBeta Lam implementation.
     Do not replace this by reconstructed native invariant expressions. *)
  expr=CheckAbort[Check[lamFunction[t3a,t3b,t3c,t3d],$Failed],$Aborted];
  If[expr===$Failed || expr===$Aborted,
    T3BCLamFail[label<>": native Lam evaluation failed"]];
  repr=ToString[expr,InputForm];
  sparse=Cases[expr,_SparseArray,Infinity];
  dims=Dimensions/@sparse;
  samples=If[Length[sparse]>0,
    (ToString[#,InputForm]&/@ Take[ArrayRules[First[sparse]],UpTo[12]]),
    {}];
  relevant=AssociationThread[names,
    (StringContainsQ[repr,#]&/@names)];
  active=And@@Values[relevant];
  verdict=If[active && Length[sparse]>0,
    "NativeLamExpressionAndSparseTensorCaptured",
    "NativeLamEvaluationRequiresInspection"];
  Print[label,": Lam tensor captured; sparse arrays=",Length[sparse],
        "; 110101 couplings visible=",active,"; status=",verdict];
  <|
    "model"->label,"status"->verdict,"native_registered_symbols"->names,
    "actual_installed_rgbeta_lam_called"->True,
    "lam_symbol"->lamName,
    "lam_expression_head"->ToString[Head[expr],InputForm],
    "lam_expression_leaf_count"->LeafCount[expr],
    "lam_expression_characters"->StringLength[repr],
    "lam_expression_sha256"->IntegerString[Hash[repr,"SHA256"],16,64],
    "lam_expression_excerpt_inputform"->StringTake[repr,Min[4000,StringLength[repr]]],
    "lam_contains_couplings"->relevant,
    "sparse_array_count"->Length[sparse],
    "sparse_array_dimensions"->dims,
    "first_sparse_array_rules_excerpt"->samples,
    "native_coupling_factors"->Lookup[records,"coupling_factor_inputform"],
    "native_unique_arrangements"->Lookup[records,"unique_arrangements_inputform"],
    "native_self_conjugate"->Lookup[records,"self_conjugate"],
    "lam_action_potential_reconstructed"->False,
    "physical_coupling_mapping_certified"->False,
    "matchete_runtime_evaluated"->False,
    "physical_one_loop_matching_certified"->False
  |>
];

result=CheckAbort[Check[Catch[
  <|
    "format"->"t3-bc-native-rgbeta-lam-diagnostic-v1",
    "status"->"LiveRGBetaLamDiagnosticOnly",
    "source"->"installed RGBeta runtime Lam",
    "wolfram_version"->$Version,
    "package_file"->ToString[FindFile["RGBeta`"],InputForm],
    "models"->Association @ Table[
      spec[[1]] -> T3BCLamProbeOne[spec[[1]],spec[[2]]],
      {spec,{{"B_-1",1},{"C_-1",3}}}
    ],
    "physical_coupling_mapping_certified"->False,
    "physical_matching_certified"->False
  |>,
"T3BCLamProbe"],$Failed],$Aborted];
If[!AssociationQ[result] || Length[Keys[result["models"]]] =!= 2,
  Print["Native RGBeta Lam diagnostic FAILED"];Exit[5]];
If[!DirectoryQ[DirectoryName[outFile]],
  CreateDirectory[DirectoryName[outFile],CreateIntermediateDirectories->True]];
written=Quiet@Check[Export[outFile,result,"RawJSON"],$Failed];
If[written===$Failed || !FileExistsQ[outFile],
  Print["Failed to write native Lam output"];Exit[6]];
Print["Native Lam diagnostic: ",outFile];
Print["Physical potential/coupling normalization and matching: NOT CERTIFIED"];

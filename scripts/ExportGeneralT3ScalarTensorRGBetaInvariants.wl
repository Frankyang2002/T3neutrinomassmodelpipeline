(* Stage 22 native RGBeta raw SU(2) invariant-expression exporter.

   Unlike the generated candidate registry, expressions here are EVALUATED
   FROM RGBeta's actually stored $quartics[symbol, Invariant] after building
   each model. Field order, native unique arrangements and native coupling
   normalisation function are exported with each expression.

   This is a DIAGNOSTIC expression export, not common-basis components,
   not an independent physical transformation proof. In particular native
   delS2 and Matchete's triplet CGs require a verified SU(2) intertwiner.

   wolframscript -file scripts/ExportGeneralT3ScalarTensorRGBetaInvariants.wl \
     output/general_t3_stage22/rgbeta_native_invariants.json
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: ExportGeneralT3ScalarTensorRGBetaInvariants.wl OUTPUT.json"];
  Exit[64]
];
outPath = ExpandFileName[$ScriptCommandLine[[2]]];
root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
registryPath = FileNameJoin[{root,"RGE","general","generated","general_rgbeta_quartic_registry.json"}];
reg = Quiet @ Check[Import[registryPath,"RawJSON"],$Failed];
If[!AssociationQ[reg] || Length[Keys[Lookup[reg,"models",<||>]]] =!= 16,
  Print["Incomplete general quartic registry"]; Exit[2]
];
load = CheckAbort[Check[Needs["RGBeta`"]; True,$Failed],$Aborted];
If[load =!= True, Print["RGBeta unavailable"]; Exit[3]];
Get[FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
Get[FileNameJoin[{root,"scripts","RGBetaNativeAccess.wl"}]];
T3S22NativeFail[msg_] := (Print[msg];Throw[$Failed,"T3Stage22Native"]);

T3S22RGBetaInvariantEvidence[] := Module[
  {models=<||>,entry,meta,build,ext,new,old,ordered,syms,info,inv,
   expr,rows,all=0, registered, symbols, canonical, dims},
  Do[
    Print["Stage 22 native RGBeta invariant expression: ", label];
    entry = reg["models",label];
    meta = Join[entry["dimensions"],{entry["alpha"]}];
    build = CheckAbort[Check[Apply[T3RGBetaBuild,meta],$Failed],$Aborted];
    If[!AssociationQ[build], T3S22NativeFail["RGBeta model build failed: "<>label]];
    ext = FileNameJoin[{root,"RGE","general","generated",label<>"_extra_quartics.wl"}];
    T3GeneralQuarticAddedSymbols = {};
    registered = CheckAbort[Check[Get[ext];True,$Failed],$Aborted];
    If[registered =!= True,
      T3S22NativeFail["RGBeta extension load failed: "<>label]
    ];
    new = Lookup[entry["new_rgbeta_symbols"],"name",{}];
    old = Lookup[entry["preexisting_rgbeta_symbols"],"existing_name",{}];
    ordered = Join[new,old];
    If[Length[ordered] =!= entry["total_additional_quartics"],
      T3S22NativeFail["Native symbol inventory mismatch: "<>label]
    ];
    syms = Symbol /@ ordered;
    If[!And@@(KeyExistsQ[$quartics,#]& /@ syms),
      T3S22NativeFail["Missing native registered quartic: "<>label]
    ];
    rows = Table[
      info = T3Stage22NativeRecord[syms[[i]]];
      inv = info["_native_invariant"];
      (* Use symbolic components, not numeric substitutions. This preserves
         native del / eps / delS2 and avoids manufacturing a triplet basis. *)
      expr = CheckAbort[Check[inv[t3i1,t3i2,t3i3,t3i4],$Failed],$Aborted];
      If[expr === $Failed || expr === $Aborted ||
         !FreeQ[expr,$Failed|$Aborted|Indeterminate|_Missing],
        T3S22NativeFail["RGBeta native invariant evaluation failed: "<>
          label<>"/"<>ordered[[i]]]
      ];
      <|
        "symbol" -> ordered[[i]],
        "native_fields_inputform" -> info["fields_inputform"],
        "native_invariant_function_inputform" -> info["invariant_inputform"],
        "native_invariant_expression_inputform" -> T3Stage22NativeString[expr],
        "native_coupling_function_inputform" -> info["coupling_factor_inputform"],
        "native_unique_arrangements_inputform" -> info["unique_arrangements_inputform"],
        "self_conjugate" -> info["self_conjugate"],
        "component_tensor_exported" -> False,
        "verified_potential_prefactor" -> Null
      |>,
      {i,Length[ordered]}
    ];
    all += Length[rows];
    AssociateTo[models,label-><|"invariants"->rows|>],
    {label,Sort[Keys[reg["models"]]]}
  ];
  If[Length[Keys[models]] =!= 16 || all =!= 50,
    T3S22NativeFail["RGBeta native expression export incomplete"]
  ];
  <|
    "format"->"stage22-native-rgbeta-invariants-v1",
    "source"->"RGBetaRegisteredQuarticInvariantFunctions",
    "status"->"NativeInvariantExpressionsNotCommonBasisTensors",
    "native_expressions_verified"->True,
    "native_invariant_component_alignment_verified"->False,
    "native_potential_prefactors_verified"->False,
    "total_slots"->all,
    "models"->models
  |>
];

If[FileExistsQ[outPath],DeleteFile[outPath]];
result=CheckAbort[Check[Catch[T3S22RGBetaInvariantEvidence[],"T3Stage22Native"],$Failed],$Aborted];
If[!AssociationQ[result] || Lookup[result,"total_slots",0] =!= 50 ||
   Length[Keys[Lookup[result,"models",<||>]]] =!= 16,
  Print["Stage 22 RGBeta native invariant expression export FAILED"]; Exit[5]
];
If[!DirectoryQ[DirectoryName[outPath]],
  CreateDirectory[DirectoryName[outPath],CreateIntermediateDirectories->True]
];
write=Quiet@Check[Export[outPath,result,"RawJSON"],$Failed];
If[write === $Failed || !FileExistsQ[outPath],
  Print["RGBeta native invariant expression JSON write failure"];Exit[6]
];
Print["Native RGBeta invariant expressions written (50 slots): ",outPath];

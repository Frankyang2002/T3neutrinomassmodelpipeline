(* Stage 22: verify *native* RGBeta quartic declaration internals.
   This extracts stored field lists, invariant pure functions, coefficient
   normalization functions and permutation metadata from RGBeta v1.2.0.
   IT IS NOT AN EXPLICIT COMPONENT TENSOR EXPORT. It cannot certify a
   Matchete/RGBeta physical mapping on its own.

   wolframscript -file scripts/ExportGeneralT3ScalarTensorRGBetaDeclarations.wl \
     output/general_t3_stage22/rgbeta_native_declarations.json
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: wolframscript -file scripts/ExportGeneralT3ScalarTensorRGBetaDeclarations.wl OUTPUT.json"];
  Exit[64]
];
outPath = ExpandFileName[$ScriptCommandLine[[2]]];
root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
registryPath = FileNameJoin[{root,"RGE","general","generated","general_rgbeta_quartic_registry.json"}];
reg = Quiet@Check[Import[registryPath,"RawJSON"],$Failed];
If[!AssociationQ[reg] || Length[Keys[Lookup[reg,"models",<||>]]] =!= 16,
  Print["RGBeta candidate registry missing/incomplete"];Exit[2]
];
load = CheckAbort[Check[Needs["RGBeta`"]; True,$Failed],$Aborted];
If[load =!= True,Print["RGBeta package unavailable"];Exit[3]];
Get[FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
Get[FileNameJoin[{root,"scripts","RGBetaNativeAccess.wl"}]];

T3Stage22InspectRGBeta[] := Module[{rows=<||>, label, settings, meta,
  old, new, expected, symbols, info, record, build, ext, registered},
  Do[
    Print["Stage 22 RGBeta declaration probe: ",label];
    settings = reg["models",label];
    meta = Join[settings["dimensions"],{settings["alpha"]}];
    build = CheckAbort[Check[Apply[T3RGBetaBuild,meta],$Failed],$Aborted];
    If[!AssociationQ[build],Print["Failed RGBeta model build: ",label]; Return[$Failed]];
    ext = FileNameJoin[{root,"RGE","general","generated",label<>"_extra_quartics.wl"}];
    T3GeneralQuarticAddedSymbols = {};
    registered = CheckAbort[Check[Get[ext];True,$Failed],$Aborted];
    If[registered =!= True,Print["Failed RGBeta extra vertices: ",label];Return[$Failed]];
    new = Lookup[settings["new_rgbeta_symbols"],"name",{}];
    old = Lookup[settings["preexisting_rgbeta_symbols"],"existing_name",{}];
    expected = Join[new,old];
    If[Length[expected] =!= settings["total_additional_quartics"],
      Print["Incomplete slot inventory: ",label];Return[$Failed]
    ];
    symbols = Symbol /@ expected;
    If[!And@@(KeyExistsQ[$quartics,#]& /@ symbols),
      Print["Missing native RGBeta stored quartic: ",label];Return[$Failed]
    ];
    record = Table[
      info = T3Stage22NativeRecord[symbols[[i]]];
      <|
        "symbol"->expected[[i]],
        "fields_inputform"->info["fields_inputform"],
        "invariant_inputform"->info["invariant_inputform"],
        "coupling_factor_inputform"->info["coupling_factor_inputform"],
        "unique_arrangements_inputform"->info["unique_arrangements_inputform"],
        "self_conjugate"->info["self_conjugate"]
      |>,
      {i,Length[expected]}
    ];
    AssociateTo[rows,label-><|"tensors_not_exported"->True,
       "declarations"->record|>],
    {label,Sort[Keys[reg["models"]]]}
  ];
  <|
    "format"->"stage22-native-rgbeta-declarations-v1",
    "source"->"RGBetaRegisteredQuarticAssociation",
    "status"->"NativeDeclarationMetadataOnly",
    "matchete_rgbeta_normalisation_verified"->False,
    "models"->rows
  |>
];

If[FileExistsQ[outPath], DeleteFile[outPath]];
result = CheckAbort[Check[Catch[T3Stage22InspectRGBeta[], "T3Stage22Native"], $Failed], $Aborted];
If[!AssociationQ[result],Print["Native RGBeta declaration probe failed"];Exit[5]];
If[!DirectoryQ[DirectoryName[outPath]],
  CreateDirectory[DirectoryName[outPath],CreateIntermediateDirectories->True]
];
write = Quiet @ Check[Export[outPath,result,"RawJSON"], $Failed];
If[write === $Failed || !FileExistsQ[outPath],
  Print["RGBeta declaration JSON write failed"]; Exit[6]
];
Print["Native RGBeta declaration evidence written (resolved keys, 50 slots): ",outPath];

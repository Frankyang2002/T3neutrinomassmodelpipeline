(* Stage 22 Matchete exact SU(2) tensor exporter, revision 2.
   Native Matchete tensors have axes ONLY for non-singlet fields. JSON
   retains four field dimensions (with a singleton axis for each SU(2) singlet).
   Fail closed on incomplete records. No RGBeta normalisation is inferred.

   wolframscript -file scripts/ExportGeneralT3Stage22Matchete.wl \
     output/general_t3_stage22/matchete_tensors.json
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: ExportGeneralT3Stage22Matchete.wl OUTPUT.json"]; Exit[64]
];
outPath = ExpandFileName[$ScriptCommandLine[[2]]];
root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
registryPath = FileNameJoin[{root, "RGE", "general", "generated", "general_rgbeta_quartic_registry.json"}];
reg = Quiet @ Check[Import[registryPath, "RawJSON"], $Failed];
If[!AssociationQ[reg] || Length[Keys[Lookup[reg, "models", <||>]]] =!= 16,
  Print["Matchete registry missing or incomplete"]; Exit[2]
];
loadStatus = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]; True], $Failed], $Aborted];
If[loadStatus =!= True, Print["Matchete / FrontEnd unavailable"]; Exit[3]];
Get[FileNameJoin[{root, "Lagrangian", "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{root, "Lagrangian", "model", "LagrangianBuilder.wl"}]];

T3S22Fail[message_] := (Print[message]; Throw[$Failed, "T3Stage22Failure"]);

T3S22TensorRecord[model_Association, counts_List, index_Integer, tensor_] := Module[
  {names, dimensions, active, raw, rawDimensions, components, simplified, expectedLength},
  names = T3ExtendedSpeciesList[counts];
  dimensions = T3ExtendedFieldDimensions[model, names];
  active = Select[dimensions, # > 1 &];
  raw = Normal[tensor];
  rawDimensions = Dimensions[raw];
  expectedLength = Times @@ dimensions;
  (* Matchete's Keep index list omits SU(2) singlets. Its tensor rank is
     usually Length[active], rather than the four scalar field slots. *)
  If[rawDimensions =!= active && rawDimensions =!= dimensions &&
     !(active === {} && (rawDimensions === {} || rawDimensions === {1})),
    Print["Unexpected native tensor axes for ", counts,
      ": got ", rawDimensions, ", active=", active,
      ", complete=", dimensions];
    Return[$Failed]
  ];
  components = Flatten[{raw}];
  If[Length[components] =!= expectedLength,
    Print["Wrong native tensor component count for ", counts,
      ": got ", Length[components], ", expected ", expectedLength];
    Return[$Failed]
  ];
  simplified = FullSimplify /@ components;
  If[!AllTrue[simplified, NumericQ],
    Print["Non-numeric Matchete CG component(s) for ", counts,
      ": ", Select[simplified, !NumericQ[#] &]];
    Return[$Failed]
  ];
  <|
    "symbol" -> ("lambdaT3ExtraQ" <> StringJoin[ToString /@ counts] <>
      "I" <> ToString[index]),
    "species" -> names,
    "counts" -> counts,
    "dimensions" -> dimensions,
    "components" -> (ToString[#, InputForm] & /@ simplified),
    "native_non_singlet_dimensions" -> active,
    "native_tensor_dimensions" -> rawDimensions,
    "potential_prefactor" -> "1"
  |>
];

T3S22ExportAll[] := Module[
  {all=<||>, cls, alpha, model, settings, grouped, countsList, counts,
   expected, tensors, records, entries, rows, sm, native, total=0},
  Do[
    Print["Stage 22 Matchete tensor export: ", label];
    cls = First[StringSplit[label, "_"]];
    alpha = ToExpression[Last[StringSplit[label, "_"]]];
    settings = reg["models", label];
    model = T3ModelFromClass[cls, alpha];
    If[!AssociationQ[model], T3S22Fail["Unsupported model " <> label]];
    ResetAll[];
    sm = CheckAbort[Check[LoadModel["SM"], $Failed], $Aborted];
    If[sm === $Failed || sm === $Aborted,
      T3S22Fail["Matchete SM group registration failed for " <> label]
    ];
    grouped = Join[
      ({#["counts"], #["name"]} &) /@ settings["new_rgbeta_symbols"],
      ({#["counts"], #["existing_name"]} &) /@ settings["preexisting_rgbeta_symbols"]
    ];
    countsList = DeleteDuplicates[First /@ grouped];
    records = {};
    Do[
      counts = countsList[[q]];
      expected = Count[grouped, {counts, _}];
      native = CheckAbort[
        Check[T3ExtendedBosonicTensorBasis[model,T3ExtendedSpeciesList[counts]], $Failed],
        $Aborted
      ];
      If[!AssociationQ[native] || native["PhysicalCount"] =!= expected,
        T3S22Fail["Matchete physical invariant rank mismatch: " <> label <>
          " counts=" <> ToString[counts, InputForm] <> " expected=" <> ToString[expected]]
      ];
      tensors = native["Basis"];
      rows = Table[T3S22TensorRecord[model,counts,i,tensors[[i]]], {i,expected}];
      If[!AllTrue[rows, AssociationQ],
        T3S22Fail["Matchete tensor serialization failed: " <> label <>
          " counts=" <> ToString[counts, InputForm]]
      ];
      records = Join[records, rows],
      {q, Length[countsList]}
    ];
    If[Length[records] =!= settings["total_additional_quartics"],
      T3S22Fail["Incomplete tensor list for " <> label]
    ];
    total += Length[records];
    AssociateTo[all, label -> <|
      "tensors" -> records, "holomorphic_orientation" -> "direct"|>],
    {label, Sort[Keys[reg["models"]]]}
  ];
  If[Length[Keys[all]] =!= 16 || total =!= 50,
    T3S22Fail["Matchete export did not include all 16 models / 50 slots"]
  ];
  <|
    "format" -> "stage22-exact-tensors-v1",
    "source" -> "MatcheteInvariantTensors",
    "basis" -> "Matchete",
    "model_count" -> 16,
    "quartic_slot_count" -> total,
    "tensor_axes_convention" -> "Four species positions, singleton axes inserted for SU(2) singlets",
    "potential_prefactor_convention" -> "+lambda[] times contracted fields in PlusHc[contribution]",
    "models" -> all
  |>
];

(* Remove stale evidence before running: a failed export must never leave an
   old file to be mistaken for a successful current export. *)
If[FileExistsQ[outPath], DeleteFile[outPath]];
result = CheckAbort[
  Check[UsingFrontEnd[Catch[T3S22ExportAll[], "T3Stage22Failure"]], $Failed],
  $Aborted
];
If[!AssociationQ[result] || Lookup[result,"model_count",0] =!= 16 ||
   Lookup[result,"quartic_slot_count",0] =!= 50 ||
   Length[Keys[Lookup[result,"models",<||>]]] =!= 16,
  Print["Stage 22 Matchete export FAILED: no valid tensor file written"];
  Exit[5]
];
If[!DirectoryQ[DirectoryName[outPath]],
  CreateDirectory[DirectoryName[outPath], CreateIntermediateDirectories -> True]
];
writeStatus = Quiet @ Check[Export[outPath, result, "RawJSON"], $Failed];
If[writeStatus === $Failed || !FileExistsQ[outPath],
  Print["Failed to write Matchete tensor JSON"]; Exit[6]
];
Print["Native Matchete tensors written (16 models, 50 slots): ", outPath];

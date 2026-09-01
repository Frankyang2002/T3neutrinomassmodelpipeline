(*
  TestT3RGEExchangeExport.wl

  Build the scotogenic-like benchmark T3-B, alpha=-1, export the RGE tensor
  exchange JSON, and fail if the export contains no quartic or Yukawa data.

  This test mirrors the production RunModel.wl initialization:
    Matchete -> T3ModelCatalog -> LagrangianBuilder -> BuildT3Lagrangian
    -> T3RGETensorExport -> JSON
*)

ClearAll[
  TestFail,
  scriptDirectory,
  projectRoot,
  catalogPath,
  builderPath,
  exporterPath,
  outputDirectory,
  outputPath,
  matcheteLoaded,
  model,
  build,
  data,
  exportResult
];

TestFail[message_] := (
  Print["FAIL: ", message];
  Exit[1]
);

scriptDirectory = DirectoryName @ ExpandFileName[$InputFileName];

projectRoot = ExpandFileName @ FileNameJoin[{
  scriptDirectory,
  "..",
  ".."
}];

catalogPath = FileNameJoin[{
  projectRoot,
  "wolfram",
  "t3",
  "T3ModelCatalog.wl"
}];

builderPath = FileNameJoin[{
  projectRoot,
  "wolfram",
  "t3",
  "LagrangianBuilder.wl"
}];

exporterPath = FileNameJoin[{
  projectRoot,
  "wolfram",
  "rge",
  "T3RGETensorExport.wl"
}];

outputDirectory = FileNameJoin[{
  projectRoot,
  "wolfram",
  "output",
  "T3_B_alpha_m1"
}];

outputPath = FileNameJoin[{
  outputDirectory,
  "rge_tensor_exchange.json"
}];

If[!FileExistsQ[catalogPath],
  TestFail["T3ModelCatalog.wl not found at " <> catalogPath]
];

If[!FileExistsQ[builderPath],
  TestFail["LagrangianBuilder.wl not found at " <> builderPath]
];

If[!FileExistsQ[exporterPath],
  TestFail["T3RGETensorExport.wl not found at " <> exporterPath]
];

If[!DirectoryQ[outputDirectory],
  CreateDirectory[
    outputDirectory,
    CreateIntermediateDirectories -> True
  ]
];


(* ------------------------------------------------------------------------- *)
(* Matchete                                                                   *)
(* ------------------------------------------------------------------------- *)

Print["Loading Matchete..."];

matcheteLoaded = UsingFrontEnd[
  Needs["Matchete`"];
  True
];

If[!TrueQ[matcheteLoaded],
  TestFail["Matchete failed to load."]
];

Print["Matchete loaded successfully."];


(* ------------------------------------------------------------------------- *)
(* Project modules                                                            *)
(* ------------------------------------------------------------------------- *)

Get[catalogPath];
Get[builderPath];
Get[exporterPath];


(* ------------------------------------------------------------------------- *)
(* Scotogenic-like benchmark: T3-B, alpha=-1                                 *)
(* ------------------------------------------------------------------------- *)

model = T3ModelFromClass["B", -1];

If[model === $Failed,
  TestFail["Could not construct T3-B alpha=-1 model."]
];

Print[
  "Model: ",
  model["Class"],
  ", alpha=-1",
  ", Scalar1=",
  {model["Scalar1", "SU2"], model["Scalar1", "Y"]},
  ", Scalar2=",
  {model["Scalar2", "SU2"], model["Scalar2", "Y"]},
  ", Fermion=",
  {model["Fermion", "SU2"], model["Fermion", "Y"]}
];


(* ------------------------------------------------------------------------- *)
(* Build UV model first                                                       *)
(* ------------------------------------------------------------------------- *)

Print["Building UV model..."];

build = CheckAbort[
  UsingFrontEnd[
    BuildT3Lagrangian[model]
  ],
  $Aborted
];

If[
  !AssociationQ[build] ||
  build === $Aborted ||
  build === $Failed,
  TestFail["BuildT3Lagrangian failed."]
];

If[build["Status"] =!= "Success",
  TestFail["UV model validation failed."]
];

Print["Accepted interactions: ", build["AllowedInteractions"]];


(* ------------------------------------------------------------------------- *)
(* Build and validate exchange                                                *)
(* ------------------------------------------------------------------------- *)

Print["Building RGE tensor exchange..."];

data = CheckAbort[
  UsingFrontEnd[
    T3RGEExportAssociation[model]
  ],
  $Aborted
];

If[
  !AssociationQ[data] ||
  data === $Aborted ||
  data === $Failed,
  TestFail["T3RGEExportAssociation failed."]
];

Print[
  "Quartic component terms: ",
  Length[data["quartic_components"]]
];

Print[
  "Raw Yukawa component terms: ",
  Length[data["raw_yukawa_components"]]
];

If[
  Length[data["quartic_components"]] <= 0,
  TestFail["No scalar quartic components were exported."]
];

If[
  Length[data["raw_yukawa_components"]] <= 0,
  TestFail["No raw Yukawa components were exported."]
];


(* ------------------------------------------------------------------------- *)
(* Write JSON                                                                 *)
(* ------------------------------------------------------------------------- *)

exportResult = ExportT3RGETensors[
  model,
  outputPath
];

If[exportResult === $Failed || !FileExistsQ[outputPath],
  TestFail["RGE tensor JSON export failed."]
];

Print["Exchange JSON: ", outputPath];
Print["PASS"];
Exit[0];

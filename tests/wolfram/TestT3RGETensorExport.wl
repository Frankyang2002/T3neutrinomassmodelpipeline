(*
  TestT3RGETensorExport.wl

  Smoke test for the T3 -> RGE tensor exporter.
*)

ClearAll[
  TestFail,
  scriptDirectory,
  projectRoot,
  builderPath,
  catalogPath,
  exporterPath,
  model,
  data,
  matcheteLoaded
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

builderPath = FileNameJoin[{
  projectRoot,
  "wolfram",
  "t3",
  "LagrangianBuilder.wl"
}];

catalogPath = FileNameJoin[{
  projectRoot,
  "wolfram",
  "t3",
  "T3ModelCatalog.wl"
}];

exporterPath = FileNameJoin[{
  projectRoot,
  "wolfram",
  "rge",
  "T3RGETensorExport.wl"
}];

If[
  !FileExistsQ[builderPath],
  TestFail["LagrangianBuilder.wl not found at " <> builderPath]
];

If[
  !FileExistsQ[catalogPath],
  TestFail["T3ModelCatalog.wl not found at " <> catalogPath]
];

If[
  !FileExistsQ[exporterPath],
  TestFail["T3RGETensorExport.wl not found at " <> exporterPath]
];


(* ------------------------------------------------------------------------- *)
(* Load Matchete exactly as the production runner does                        *)
(* ------------------------------------------------------------------------- *)

Print["Loading Matchete..."];

matcheteLoaded = UsingFrontEnd[
  Needs["Matchete`"];
  True
];

If[
  !TrueQ[matcheteLoaded],
  TestFail["Matchete failed to load."]
];

Print["Matchete loaded successfully."];


(* ------------------------------------------------------------------------- *)
(* Load project modules                                                       *)
(* ------------------------------------------------------------------------- *)

Print["Loading T3 model catalog..."];
Get[catalogPath];

Print["Loading T3 builder..."];
Get[builderPath];

Print["Loading RGE tensor exporter..."];
Get[exporterPath];


(* ------------------------------------------------------------------------- *)
(* Construct benchmark model                                                  *)
(* ------------------------------------------------------------------------- *)

model = T3ModelFromClass["E", 0];

If[
  model === $Failed,
  TestFail["Could not construct T3-E alpha=0 model."]
];

Print[
  "Model: ",
  model["Class"],
  ", Scalar1 = ",
  {model["Scalar1", "SU2"], model["Scalar1", "Y"]},
  ", Scalar2 = ",
  {model["Scalar2", "SU2"], model["Scalar2", "Y"]},
  ", Fermion = ",
  {model["Fermion", "SU2"], model["Fermion", "Y"]}
];


(* ------------------------------------------------------------------------- *)
(* Build the UV model first                                                   *)
(* ------------------------------------------------------------------------- *)

Print["Building benchmark UV model..."];

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

If[
  build["Status"] =!= "Success",
  TestFail["UV model validation failed."]
];

Print["UV model built successfully."];
Print["Accepted interactions: ", build["AllowedInteractions"]];


(* ------------------------------------------------------------------------- *)
(* Export RGE component tensors                                               *)
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

Print[];
Print[
  "Quartic component terms: ",
  Length[data["quartic_components"]]
];

Print[
  "Raw Yukawa component terms: ",
  Length[data["raw_yukawa_components"]]
];

Print[
  "Status: ",
  InputForm[data["status"]]
];


(* ------------------------------------------------------------------------- *)
(* Validation                                                                 *)
(* ------------------------------------------------------------------------- *)

If[
  Length[data["quartic_components"]] <= 0,
  TestFail["No scalar quartic components were exported."]
];

If[
  Length[data["raw_yukawa_components"]] <= 0,
  TestFail["No Yukawa components were exported."]
];

Print["PASS"];
Exit[0];

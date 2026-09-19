(* ExportUVCGRegistry.wl
   Export the three UV T3 Clebsch-Gordan tensors directly from the same
   Matchete registry populated by BuildT3Lagrangian.

   Usage:
     wolframscript -file Lagrangian/interactions/ExportUVCGRegistry.wl dS1 dS2 dF alpha output.json
*)

ClearAll["Global`*"];

If[Length[$ScriptCommandLine] < 6,
  Print[
    "Usage: wolframscript -file ExportUVCGRegistry.wl ",
    "dS1 dS2 dF alpha output.json"
  ];
  Exit[2];
];

scriptDirectory = DirectoryName @ ExpandFileName[$InputFileName];
projectRoot = DirectoryName[scriptDirectory];

dS1 = ToExpression[$ScriptCommandLine[[2]]];
dS2 = ToExpression[$ScriptCommandLine[[3]]];
dF  = ToExpression[$ScriptCommandLine[[4]]];
alpha = ToExpression[$ScriptCommandLine[[5]]];
outputPath = ExpandFileName[$ScriptCommandLine[[6]]];

If[!DirectoryQ[DirectoryName[outputPath]],
  CreateDirectory[DirectoryName[outputPath], CreateIntermediateDirectories -> True]
];

Print["Loading Matchete..."];
matcheteLoaded = UsingFrontEnd[Needs["Matchete`"]; True];
If[!TrueQ[matcheteLoaded],
  Print["ERROR: Matchete failed to load."];
  Exit[3];
];

Get[FileNameJoin[{projectRoot, "Lagrangian", "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{projectRoot, "Lagrangian", "model", "LagrangianBuilder.wl"}]];

model = T3ModelFromDimensions[dS1, dS2, dF, alpha];
If[model === $Failed,
  Print["ERROR: invalid T3 dimensions."];
  Exit[4];
];

build = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
If[
  build === $Failed || build === $Aborted ||
  !AssociationQ[build] || Lookup[build, "Status", ""] =!= "Success",
  Print["ERROR: UV T3 build failed."];
  Exit[5];
];

cgNames = {T3Y1CG, T3Y2CG, T3MixCG};

cgRegistry = DeleteCases[
  Map[
    Function[name,
      Module[{reps, tensor},
        reps = Quiet @ Check[
          Matchete`PackageScope`$CGproperties[name, Indices],
          $Failed
        ];
        tensor = Quiet @ Check[
          Matchete`CGManipulations`PackagePrivate`$CGtensors[name],
          $Failed
        ];

        If[
          reps === $Failed || tensor === $Failed ||
          MissingQ[reps] || MissingQ[tensor],
          Nothing,
          <|
            "Name" -> ToString[Unevaluated[name], InputForm],
            "RepsInputForm" -> ToString[InputForm[reps]],
            "TensorInputForm" -> ToString[InputForm[tensor]]
          |>
        ]
      ]
    ],
    cgNames
  ],
  Nothing
];

payload = <|
  "status" -> If[Length[cgRegistry] === 3, "Success", "Incomplete"],
  "metadata" -> <|
    "dS1" -> dS1,
    "dS2" -> dS2,
    "dF" -> dF,
    "alpha" -> alpha,
    "YS1" -> ToString[InputForm[model["Scalar1", "Y"]]],
    "YS2" -> ToString[InputForm[model["Scalar2", "Y"]]],
    "YF" -> ToString[InputForm[model["Fermion", "Y"]]]
  |>,
  "CGRegistry" -> cgRegistry
|>;

Export[outputPath, payload, "RawJSON"];

Print["UV CG registry export: ", payload["status"]];
Print["Output: ", outputPath];
Print["CG names: ", Lookup[cgRegistry, "Name", {}]];

If[payload["status"] =!= "Success", Exit[6]];

(*
  GeneralRGE.wl

  Generalized RGBeta entry point for the T3 pipeline.

  RGE.wl is retained unchanged as legacy/scotogenic-era code.  This file is
  intentionally separate so the generalized implementation can be built and
  validated without changing the old reference.

  Arguments:
    1. Input model JSON path
    2. Output JSON path

  Current scope:
    - load and validate a generalized model description;
    - reset RGBeta;
    - establish a clean entry point for model definition and beta-function
      export.

  The RGBeta model-definition and beta-function extraction are still to be
  implemented.  Until those steps exist, this script exports a status JSON
  rather than pretending that an RGE calculation has been completed.
*)

<< RGBeta

ClearAll[
  Fail,
  args,
  modelPath,
  outputPath,
  model,
  modelName,
  requiredKeys,
  missingKeys,
  result
];

Fail[message_] := (Print["ERROR: ", message]; Exit[1]);

args = Rest[$ScriptCommandLine];

If[
  Length[args] < 2,
  Fail["Expected input-model and output paths."]
];

{modelPath, outputPath} = Take[args, 2];

If[
  !FileExistsQ[modelPath],
  Fail["Model file does not exist: " <> modelPath]
];

model = Quiet@Check[
  Import[modelPath, "RawJSON"],
  $Failed
];

If[
  model === $Failed || !AssociationQ[model],
  Fail["Could not import model JSON."]
];

requiredKeys = {"model_name"};
missingKeys = Select[
  requiredKeys,
  !KeyExistsQ[model, #] &
];

If[
  missingKeys =!= {},
  Fail[
    "Input JSON is missing required key(s): " <>
    StringRiffle[missingKeys, ", "]
  ]
];

modelName = model["model_name"];

Print["Loaded generalized RGE model: ", modelName];

ResetModel[];

result = <|
  "model_name" -> modelName,
  "status" -> "FoundationReady",
  "rgbeta_model_defined" -> False,
  "beta_functions_generated" -> False,
  "message" ->
    "Generalized RGE entry point is ready; RGBeta model definition and beta extraction remain to be implemented."
|>;

Quiet@Check[
  Export[outputPath, result, "RawJSON"],
  Fail["Could not export RGE status JSON: " <> outputPath]
];

Print["RGE status written to: ", outputPath];

(*
  RGE.wl

  Not Done, not used yet as well
  Command-line entry point for RGBeta-based RGE generation.

  Arguments:
    1. Input model JSON path
    2. Output JSON path

  At present this script only loads the model and resets RGBeta. The actual
  model-definition, beta-function calculation and export steps remain to be
  implemented.
*)

<< RGBeta

ClearAll[Fail, args, modelPath, outputPath, model, modelName];

(* Print a clear error and terminate with a non-zero process status. *)
Fail[message_] := (Print["ERROR: ", message]; Exit[1]);

args = Rest[$ScriptCommandLine];
If[Length[args] < 2, Fail["Expected input-model and output paths."]];

{modelPath, outputPath} = Take[args, 2];
If[!FileExistsQ[modelPath], Fail["Model file does not exist: " <> modelPath]];

model = Quiet@Check[Import[modelPath, "RawJSON"], $Failed];
If[model === $Failed || !AssociationQ[model], Fail["Could not import model JSON."]];

modelName = Lookup[model, "model_name", Missing["model_name"]];
If[MissingQ[modelName], Fail["Input JSON has no \"model_name\" key."]];

Print["Loaded model: ", modelName];
ResetModel[];

(* outputPath is intentionally retained for the future export step. *)

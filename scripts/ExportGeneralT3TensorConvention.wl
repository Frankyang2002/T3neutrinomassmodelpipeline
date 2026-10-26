(* Experimental SU(2) tensor diagnostics for the two C5-active B/C alpha=-1
   Higgs quartics. This does NOT identify an RGBeta coupling value. It reports
   the exact Matchete tensor relative to canonical fundamental-index tensors,
   plus a residual. Matchete convention may be different; in that case the
   script fails closed rather than assigning an unjustified normalisation.

   Initialise Matchete SU2L by LoadModel["SM"] AFTER each ResetAll[],
   before the BSM SU(2) tensors are generated. This is mandatory for the
   standalone exporter, just as in the production Lagrangian builder.

   wolframscript -file scripts/ExportGeneralT3TensorConvention.wl \
      output/general_t3_trial/matchete_active_quartic_tensors.json
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 1,
  Print["Usage: ExportGeneralT3TensorConvention.wl OUTPUT.json"];
  Exit[64]
];
outputPath = ExpandFileName[First[args]];
projectRoot = DirectoryName[DirectoryName @ ExpandFileName[$InputFileName]];

loadResult = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]; True], $Failed], $Aborted];
If[loadResult =!= True, Print["Cannot load Matchete with FrontEnd."]; Exit[3]];
Get[FileNameJoin[{projectRoot, "Lagrangian", "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{projectRoot, "Lagrangian", "model", "LagrangianBuilder.wl"}]];

(* Tensor positions are Hdag[a], Hdag[b], S[c], S[d].
   The first structure has conjugate and fundamental SU(2) indices: deltas.
   The second has four conjugate fundamental SU(2) indices: epsilons.
   The factors of 1/2 are a canonical convention, not Matchete's convention.
*)
delta[a_Integer, b_Integer] := KroneckerDelta[a, b];
eps[a_Integer, b_Integer] := Which[a == 1 && b == 2, 1,
  a == 2 && b == 1, -1, True, 0];
refDelta = Table[(delta[a, c] delta[b, d] + delta[a, d] delta[b, c])/2,
  {a, 2}, {b, 2}, {c, 2}, {d, 2}];
refEps = Table[(eps[a, c] eps[b, d] + eps[a, d] eps[b, c])/2,
  {a, 2}, {b, 2}, {c, 2}, {d, 2}];

T3TensorDescribe[model_Association, counts_List, reference_, name_String] := Module[
  {names, data, basis, tensor, refVector, rawVector, factor, residual,
   consistent, dims},
  names = T3ExtendedSpeciesList[counts];
  dims = T3ExtendedFieldDimensions[model, names];
  data = CheckAbort[Check[T3ExtendedBosonicTensorBasis[model, names], $Failed], $Aborted];
  If[!AssociationQ[data], Return[<|"status" -> "GenerationFailed",
    "coupling" -> name|>]];
  If[data["PhysicalCount"] =!= 1, Return[<|"status" -> "UnexpectedMultiplicity",
    "coupling" -> name, "count" -> data["PhysicalCount"]|>]];
  tensor = Normal[First[data["Basis"]]];
  If[Dimensions[tensor] =!= {2, 2, 2, 2},
    Return[<|"status" -> "UnexpectedTensorDimensions", "coupling" -> name,
      "dimensions" -> Dimensions[tensor]|>]
  ];
  refVector = Flatten[reference];
  rawVector = Flatten[tensor];
  factor = FullSimplify[Total[Conjugate[refVector] rawVector]/
    Total[Conjugate[refVector] refVector]];
  residual = FullSimplify[rawVector - factor refVector];
  consistent = TrueQ[And @@ (PossibleZeroQ /@ residual)] &&
    !TrueQ[PossibleZeroQ[factor]];
  <|
    "status" -> If[consistent, "CanonicalTensorProportional", "BasisMismatchReview"],
    "coupling" -> name,
    "species" -> names,
    "field_counts" -> counts,
    "representation_dimensions" -> dims,
    "matchete_tensor_dimensions" -> Dimensions[tensor],
    "canonical_reference" -> If[reference === refDelta,
      "(delta_ac delta_bd+delta_ad delta_bc)/2",
      "(eps_ac eps_bd+eps_ad eps_bc)/2"],
    "matchete_over_canonical_factor_inputform" -> ToString[factor, InputForm],
    "all_component_residuals_zero" -> consistent,
    "raw_tensor_inputform" -> ToString[tensor, InputForm],
    "canonical_tensor_inputform" -> ToString[reference, InputForm]
  |>
];

T3TensorRun[] := Module[{result = {}, model, elements, fail, smBaseline},
  Do[
    (* ResetAll[] removes Matchete's registered SM groups.  The BSM tensor
       machinery calls EnsureBSMRepresentation[], which in turn calls
       DefineRepresentation[rep, SU2L, ...].  SU2L MUST be restored first.
       This mirrors BuildT3Lagrangian in LagrangianBuilder.wl. *)
    ResetAll[];
    smBaseline = CheckAbort[Check[LoadModel["SM"], $Failed], $Aborted];
    If[smBaseline === $Failed || smBaseline === $Aborted,
      Print["ERROR: LoadModel[\"SM\"] failed before tensor generation for ", class];
      Return[<|
        "status" -> "SMModelLoadFailed",
        "failed_model" -> class <> "_m1",
        "results" -> result,
        "rgbeta_mapping_verified" -> False
      |>]
    ];
    Print["Loaded SM gauge registry before BSM tensor generation for ", class];
    model = T3ModelFromClass[class, -1];
    If[model === $Failed, Return[$Failed]];
    (* Models use complex S1,S2 and independent scalar fields. *)
    elements = {
      T3TensorDescribe[model, {0, 2, 0, 0, 2, 0}, refDelta,
        "lambdaT3ExtraQ020020I1"],
      T3TensorDescribe[model, {0, 2, 0, 2, 0, 0}, refEps,
        "lambdaT3ExtraQ020200I1"]
    };
    AppendTo[result, <|"model" -> class <> "_m1", "tensors" -> elements|>],
    {class, {"B", "C"}}
  ];
  fail = AnyTrue[Flatten[Lookup[result, "tensors", {}], 1],
    Lookup[#, "status", ""] =!= "CanonicalTensorProportional" &];
  <|
    "status" -> If[fail, "ReviewRequired", "MatcheteTensorConventionVerified"],
    "scope" -> "two C5-active structures for B/C(-1)",
    "comparison" -> "Matchete tensor relative to explicit canonical doublet tensor",
    "rgbeta_mapping_verified" -> False,
    "note" -> "The ratio shown is a Matchete-to-canonical SU2 factor, NOT a Matchete-to-RGBeta coupling factor. RGBeta delS2 normalisation and Lagrangian sign remain to be checked.",
    "results" -> result
  |>
];

result = CheckAbort[UsingFrontEnd[T3TensorRun[]], $Aborted];
If[!AssociationQ[result], Print["Tensor extraction aborted/failed."]; Exit[5]];
If[!DirectoryQ[DirectoryName[outputPath]],
  CreateDirectory[DirectoryName[outputPath], CreateIntermediateDirectories -> True]
];
Export[outputPath, result, "RawJSON"];
Print["Status: ", result["status"]];
Scan[Function[item, Scan[Function[t,
  Print[item["model"], " ", t["coupling"], " ", t["status"],
    " Matchete/canonical=", Lookup[t,
      "matchete_over_canonical_factor_inputform", "unknown"]]],
  item["tensors"]]], result["results"]];
Print["Report: ", outputPath];
If[result["status"] =!= "MatcheteTensorConventionVerified", Exit[6]];
Exit[0];

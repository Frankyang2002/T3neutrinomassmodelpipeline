(* Native Matchete v0.6.0 A_-4 Feynman-rule extraction.
   Deliberately does NOT claim to evaluate a connected UV four-point amplitude:
   FeynmanRules returns vertices, not a contracted heavy propagator exchange.
   Writes new evidence ONLY below output/general_observable/lo_fermion_boundary/.
   Launch from the repository with wolframscript -file PATH.
*)
ClearAll["Global`*"];
$T3NativeProbeVersion = "t3-native-matchete-fr-evidence-v1";
probeFile = ExpandFileName[$InputFileName];
repoRoot = Nest[DirectoryName, probeFile, 4];
probeDir = DirectoryName[probeFile];
sourceFile = FileNameJoin[{repoRoot, "output", "general_observable",
   "general_wilson_boundary_v1", "A_m4", "wilson_boundary.json"}];
probeOut = FileNameJoin[{probeDir, "normalization_probes"}];
If[!DirectoryQ[probeOut], CreateDirectory[probeOut, CreateIntermediateDirectories -> True]];
reportPath = FileNameJoin[{probeOut, "native_feynman_a_m4_" <>
   StringReplace[CreateUUID[], "-" -> ""] <> ".json"}];
ClearAll[writeProbe, failProbe, termsOf, description, inputFormOf];
writeProbe[rec_Association] := Module[{payload},
  payload = Join[<|
    "format" -> $T3NativeProbeVersion,
    "model" -> "A_-4",
    "matchete_gamma_cc_to_weyl_phase_certified" -> False,
    "independent_connected_uv_four_point_amplitude_computed" -> False,
    "physical_master_wilson_normalization_certified" -> False,
    "original_16_model_release_changed" -> False
    |>, rec];
  Export[reportPath, payload, "RawJSON"];
  Print["NATIVE_FR_REPORT=", reportPath];
  payload
];
failProbe[phase_String, why_String] := (
  writeProbe[<|"status" -> "Blocked", "blocked_phase" -> phase,
    "reason" -> why|>];
  Print["NATIVE_FR_STATUS=Blocked; PHASE=", phase, "; REASON=", why];
  Exit[2]
);
termsOf[e_] := Module[{x = Expand[e]}, If[Head[x] === Plus, List @@ x, {x}]];
inputFormOf[e_] := ToString[InputForm[e]];
(* Retain complete vertex expressions when reasonably sized. Longer results
   carry a digest/size and a diagnostic excerpt, never a false 'complete' flag. *)
description[e_] := Module[{s = inputFormOf[e], limit = 100000, h},
  h = Hash[s, "SHA256", "HexString"];
  <|"head" -> ToString[Head[e], InputForm],
    "inputform_sha256" -> h,
    "inputform_characters" -> StringLength[s],
    "inputform_complete" -> (StringLength[s] <= limit),
    "inputform" -> If[StringLength[s] <= limit, s, Null],
    "inputform_prefix" -> StringTake[s, UpTo[Min[6000, StringLength[s]]]]|>
];

Print["NATIVE_FR_PROBE_VERSION=", $T3NativeProbeVersion];
Print["SOURCE_ROOT=", repoRoot];
If[!FileExistsQ[sourceFile], failProbe["source", "Native A_-4 source JSON unavailable"]];
native = Quiet[Check[Import[sourceFile, "RawJSON"], $Failed]];
If[!AssociationQ[native] || Lookup[native, "model", ""] =!= "A_-4",
  failProbe["source", "Native source JSON invalid/model mismatch"]];
sources = Lookup[native, "operator_terms", {}];
source = Select[sources, Lookup[#, "source_index", None] === 2 &];
If[Length[source] =!= 1,
  failProbe["source", "Expected exactly one native A_-4 source index 2"]];
sourceText = Lookup[First[source], "source_term_inputform", ""];
If[!StringQ[sourceText] ||
   !StringContainsQ[sourceText, "GammaCC, Proj[-1]"] ||
   !StringContainsQ[sourceText, "Coupling[MF,"] ||
   !StringContainsQ[sourceText, "CG[T3MixCG"] ||
   StringContainsQ[sourceText, "GammaCC, Proj[1]"],
  failProbe["source", "Source term fails holomorphic GammaCC or mass/CG checks"]];

(* The checked repository loads Matchete through the frontend for proper
   package initialization, as independently established by the V2 probe. *)
load = CheckAbort[UsingFrontEnd[Needs["Matchete`"]; True], $Aborted];
If[load =!= True || Names["*`GammaCC"] === {} ||
    Names["*`FeynmanRules"] === {},
  failProbe["load", "Matchete GammaCC/FeynmanRules unavailable"]];
Print["NATIVE_FR_MATCHeTE_LOADED=True"];

sourceDir = FileNameJoin[{repoRoot, "Lagrangian"}];
Get[FileNameJoin[{sourceDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{sourceDir, "model", "LagrangianBuilder.wl"}]];
If[!And @@ (Environment[#] === "1" & /@ {
  "T3_EXPERIMENTAL_GENERAL_FERMIONS",
  "T3_EXPERIMENTAL_GENERAL_QUARTICS",
  "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
  "T3_EXPERIMENTAL_FLAVOR_MASSES"}),
  failProbe["environment", "Full-general flags must all be 1"]];
model = T3ModelFromClass["A", -4, 3];
If[!AssociationQ[model] || SetT3HeavyFields[{"F"}] =!= True,
  failProbe["model", "Cannot prepare A_-4 with F heavy and scalars light"]];
Print["NATIVE_FR_STAGE=BuildT3Lagrangian"];
built = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
If[!AssociationQ[built] || Lookup[built, "Status", ""] =!= "Success" ||
    Lookup[built, "FullValidation", False] =!= True ||
    Lookup[built, "FermionYukawaInventoryAssembled", False] =!= True,
  failProbe["model", "Complete A_-4 UV Lagrangian failed validation"]];
If[Lookup[built, "AllScalarExtensionsRequested", False] =!= True ||
    Lookup[built, "AdditionalQuadraticsIncluded", -1] =!=
      Lookup[built, "AdditionalQuadraticsExpected", -2] ||
    Lookup[built, "AdditionalQuarticsIncluded", -1] =!=
      Lookup[built, "AdditionalQuarticsExpected", -2],
  failProbe["model", "Full-general scalar inventory inconsistent"]];

(* Restrict the UV Feynman-rule calculation to the TWO topology Yukawa
   vertices plus their Hermitian conjugates. The other couplings do not enter
   the requested tree exchange; the full model has already been validated. *)
uvYukawaTerms = Select[termsOf[built["LInt"]], Function[t,
  Module[{s = inputFormOf[t]},
    (StringContainsQ[s, "Coupling[y1,"] ||
     StringContainsQ[s, "Coupling[y2,"]) &&
    StringContainsQ[s, "NewFermion"] &&
    StringContainsQ[s, "Field[l,"]
  ]
]];
uvY1 = Count[uvYukawaTerms, t_ /;
   StringContainsQ[inputFormOf[t], "Coupling[y1,"]];
uvY2 = Count[uvYukawaTerms, t_ /;
   StringContainsQ[inputFormOf[t], "Coupling[y2,"]];
If[uvY1 =!= 2 || uvY2 =!= 2 || Length[uvYukawaTerms] =!= 4,
  failProbe["UV vertices", "Expected two topology Yukawa terms and each HC"]];
uvYukawa = RelabelIndices[Total[uvYukawaTerms]];

(* Import only the previously verified Matchete expression, not a manually
   reconstructed spinor operator. The Python preflight verifies its SHA chain. *)
contact = Quiet[Check[ToExpression[sourceText, InputForm], $Failed]];
If[contact === $Failed || FreeQ[contact, GammaCC, Infinity] ||
    FreeQ[contact, Proj[-1], Infinity],
  failProbe["EFT1 vertex", "Could not parse native PL contact operator"]];
Print["NATIVE_FR_STAGE=FeynmanRules_UV_3point"];
frOptions = inputFormOf[Options[FeynmanRules]];
uvRules = CheckAbort[UsingFrontEnd[FeynmanRules[uvYukawa]], $Aborted];
If[uvRules === $Aborted || uvRules === $Failed ||
   Head[uvRules] === FeynmanRules || uvRules === Null || uvRules === {},
  failProbe["UV Feynman rules", "No nonempty native UV three-point vertices"]];
Print["NATIVE_FR_STAGE=FeynmanRules_EFT1_4point"];
eftRules = CheckAbort[UsingFrontEnd[FeynmanRules[contact]], $Aborted];
If[eftRules === $Aborted || eftRules === $Failed ||
   Head[eftRules] === FeynmanRules || eftRules === Null || eftRules === {},
  failProbe["EFT1 Feynman rules", "No nonempty native EFT1 LLSS contact vertex"]];

report = <|
  "status" -> "NativeUV3pointAndEFT1ContactVerticesExtracted_AmplitudeOpen",
  "uv_bsm_yukawa_terms" -> Length[uvYukawaTerms],
  "uv_y1_with_hc" -> uvY1,
  "uv_y2_with_hc" -> uvY2,
  "native_source_index" -> 2,
  "native_source_inputform" -> sourceText,
  "feynman_rule_options" -> frOptions,
  "uv_yukawa_lagrangian" -> description[uvYukawa],
  "uv_three_point_rules" -> description[uvRules],
  "eft1_four_point_rules" -> description[eftRules],
  "uv_mass_lagrangian_source" -> description[FreeLag[NewFermion]],
  "uv_three_point_rules_extracted" -> True,
  "eft1_four_point_rules_extracted" -> True,
  "heavy_propagator_and_both_exchange_orderings_contracted" -> False,
  "why_no_absolute_sign" ->
    "FeynmanRules returned vertices, not the connected two-vertex UV amplitude; spinor-flow and heavy-propagator contraction must be evaluated independently"
|>;
writeProbe[report];
Print["NATIVE_FR_STATUS=", report["status"]];
Print["NATIVE_FR_UV_HEAD=", Head[uvRules]];
Print["NATIVE_FR_EFT_HEAD=", Head[eftRules]];
Exit[0];
